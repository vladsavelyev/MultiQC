from __future__ import division
import json
import logging
import operator
import os
from collections import OrderedDict, defaultdict
from itertools import islice

from multiqc import config
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

log = logging.getLogger(__name__)

GENOME_SIZE = 3200000000
EXOME_SIZE  = 320000000

read_format = '{:,.1f} ' + config.read_count_prefix
if config.read_count_multiplier == 1:
    read_format = '{:,.0f}'


class Metrics:
    def __init__(self):
        self.reads = 0
        self.bases = 0
        self.perfect_index = 0
        self.yield_q30 = 0
        self.qscore_sum = 0
        self.unknown_barcodes = 0
        self.trimmed_bases = 0
        self.read_stats = defaultdict(int)

    def generate_metics_dict(self, total_stats):
        data = {
            'total': self.reads,
            'total_yield': self.bases,
            'clusters_pct': self.reads / total_stats.reads * 100.0,
            'yield_pct': self.bases / total_stats.bases * 100.0,
        }

        if self.yield_q30:
            data['yieldQ30'] = self.yield_q30
            if self.bases > 0:
                data['percent_Q30'] = float(self.yield_q30) / float(self.bases) * 100.0

            data['depth'] = self.yield_q30 / GENOME_SIZE

        if self.perfect_index and self.reads > 0:
            data['percent_perfectIndex'] = '{0:.1f}'.format(float(100.0 * self.perfect_index / self.reads))

        if self.qscore_sum and self.bases > 0:
            data['mean_qscore'] = float(self.qscore_sum) / float(self.bases) * 100.0

        return data


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        self.lanes = defaultdict(Metrics)
        self.samples = defaultdict(Metrics)
        self.sample_by_lane = defaultdict(lambda: defaultdict(Metrics))
        self.total_stats = Metrics()
        self.undetermined = Metrics()
        self.undetermined_by_lane = defaultdict(Metrics)

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='bcl2fastq',
            anchor='bcl2fastq',
            href="https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html",
            info="can be used to both demultiplex data and convert BCL files"
                 " to FASTQ file formats for downstream analysis."
        )

        # Gather data from all json files
        found_files = [f for f in self.find_log_files('bcl2fastq')]
        if len(found_files) == 0:
            raise UserWarning

        source_path_per_sample = self._merge_bcl2fastq_runs(found_files)

        if len(found_files) > 1:
            log.warning("Detected multiple bcl2fastq runs from the same sequencer output. "
                        "They will be merged, undetermined stats will recalculated, "
                        "unknown barcodes stats will be dropped.")
        self._recalculate_undetermined()

        # Filter to strip out ignored sample names
        self.samples = self.ignore_samples(self.samples)
        self.sample_by_lane = {lid: self.ignore_samples(samples) for lid, samples in self.sample_by_lane.items()}

        # Return with Warning if no files are found
        if len(self.samples) == 0:
            raise UserWarning

        # Print source files
        for s in source_path_per_sample.keys():
            self.add_data_source(
                s_name=s,
                source=",".join(list(set(source_path_per_sample[s]))),
                module='bcl2fastq',
                section='bcl2fastq-bysample'
            )

        # Add sample counts to general stats table
        samples_table = self._make_samples_table()
        self.add_general_stats(samples_table)
        self.write_data_file(samples_table, 'multiqc_bcl2fastq_bysample')

        # Add section for summary stats per flow cell
        lanes_table = self._make_lanes_table()
        self.add_section (
            name='Lane Statistics',
            anchor='bcl2fastq-lanestats',
            description='Statistics about each lane for each flowcell',
            plot=self.lane_stats_table(lanes_table)
        )
        self.write_data_file(lanes_table, 'multiqc_bcl2fastq_bylane')

        # Add section for counts by lane
        categories_readtype = OrderedDict()
        categories_readtype['perfect'] = {'name': 'Perfect Index Reads'}
        categories_readtype['imperfect'] = {'name': 'Mismatched Index Reads'}
        categories_readtype['undetermined'] = {'name': 'Undetermined Reads'}
        self.add_section (
            name='Read clusters by lane',
            anchor='bcl2fastq-bylane',
            description='Number of read clusters (e.g. pairs) per lane, '
                        'categorised as "perfect index", "mismatched index" or "undetermined".',
            helptext="""Perfect index reads are those that do not have a single mismatch.
                All samples of a lane are combined. Undetermined reads are treated as a third category.""",
            plot=bargraph.plot(
                self.get_counts_bar_data_per_category(self.lanes, self.undetermined_by_lane),
                categories_readtype,
                {
                    'id': 'bcl2fastq_lane_counts',
                    'title': 'bcl2fastq: Clusters by lane',
                    'ylab': 'Number of clusters',
                    'hide_zero_cats': False
                }
            )
        )

        # Add section for counts by sample
        # lcats = set()  # get cats for per-lane tab
        # for lid, samples in self.sample_by_lane.items():
        #     lcats.update(samples.keys())
        # lcats = list(sorted(lcats))
        categories_lanes = list(sorted(set(self.sample_by_lane.keys())))
        self.add_section(
            name='Clusters by sample',
            anchor='bcl2fastq-bysample',
            description='Number of reads per sample.',
            helptext="""Perfect index reads are those that do not have a single mismatch.
                All samples are aggregated across lanes combinned. Undetermined reads are ignored.
                Undetermined reads are treated as a separate sample.""",
            plot=bargraph.plot(
                [
                    MultiqcModule.get_counts_bar_data_per_category(self.samples),
                    MultiqcModule.get_perfect_counts_bar_data_per_lane(self.sample_by_lane)
                ],
                [categories_readtype, categories_lanes],
                {
                    'id': 'bcl2fastq_sample_counts',
                    'title': 'bcl2fastq: Clusters by sample',
                    'hide_zero_cats': False,
                    'ylab': 'Number of clusters',
                    'data_labels': ['Index mismatches', 'Counts per lane']
                }
            )
        )

        if len(found_files) == 1:
            # Add section with undetermined barcodes
            self.add_section(
                name="Undetermined barcodes by lane",
                anchor="undetermine_by_lane",
                description="Count of the top twenty most abundant undetermined barcodes by lanes",
                plot=bargraph.plot(
                    MultiqcModule.get_bar_data_from_undetermined(self.lanes),
                    None,
                    {
                        'id': 'bcl2fastq_undetermined',
                        'title': 'bcl2fastq: Undetermined barcodes by lane',
                        'ylab': 'Count',
                        'tt_percentages': False,
                        'use_legend': True,
                        'tt_suffix': 'reads'
                    }
                )
            )

    def _merge_bcl2fastq_runs(self, found_files):
        source_path_per_sample = defaultdict(set)

        # 225957331   = sum(dR['NumberReads'] for dR in conversionResult.get("DemuxResults", []))
        # 2699078367  = conversionResult['Undetermined']['NumberReads']
        # 2925035698  = conversionResult['Undetermined']['NumberReads'] + sum([dR['NumberReads'] for dR in conversionResult.get("DemuxResults", [])])
        # 2925035698  = conversionResult['TotalClustersPF']
        # 3830022144  = conversionResult['TotalClustersRaw']
        # 302.0       = conversionResult['Yield'] / conversionResult['TotalClustersPF']
        # 230.641168  = conversionResult['Yield'] / conversionResult['TotalClustersRaw']
        #
        # Yield = 302*TotalClustersPF, TotalClustersPF = sum(NumberReads for all Sample and Undetermined)

        contents_by_sequencer_run = defaultdict(list)

        for file_i, myfile in enumerate(found_files):
            try:
                content = json.loads(myfile["f"])
            except ValueError:
                log.warning('Could not parse file as json: {}'.format(myfile["fn"]))
                return

            sequencer_run_id = content["RunId"]
            filename = os.path.join(myfile['root'], myfile["fn"])
            contents_by_sequencer_run[sequencer_run_id].append((content, filename))

        if len(contents_by_sequencer_run) > 1:
            log.error("Input data belongs to multiple sequencer runs. It's not supported by the bcl2fastq module. " +
                      "Please run MultiQC on each run separately to generate separate reports. "
                      "Files: " + str(found_files), " runs: " + str(contents_by_sequencer_run.keys()))
            raise UserWarning

        sequencer_run_id, contents = list(contents_by_sequencer_run.items())[0]

        for file_i, (content, filename) in enumerate(contents):
            for lane_data in content.get("ConversionResults", []):
                lane_number = lane_data["LaneNumber"]
                lane_id = MultiqcModule.prepend_runid(sequencer_run_id, 'L{}'.format(lane_number))
                if lane_id in self.lanes:
                    log.debug("Duplicate runId/lane combination found! Overwriting: {}".format(lane_id))
                lane_d = self.lanes[lane_id]
                lane_d.reads = lane_data['TotalClustersPF']
                lane_d.bases = lane_data['Yield']
                if file_i == 0:
                    self.total_stats.reads += lane_d.reads
                    self.total_stats.bases += lane_d.bases
                else:
                    assert lane_d.reads == lane_data['TotalClustersPF'], (lane_d.reads, lane_data['TotalClustersPF'])
                    assert lane_d.bases == lane_data['Yield'], (lane_d.bases, lane_data['Yield'])

                for sample_data in lane_data.get("DemuxResults", []):
                    if sample_data["SampleId"] == sample_data["SampleName"]:
                        sample_id = sample_data["SampleName"]
                    else:
                        sample_id = "{}-{}".format(sample_data["SampleId"], sample_data["SampleName"])
                    if sample_id in self.sample_by_lane[lane_id]:
                        log.debug("Duplicate runId/lane/sample combination found! "
                                  "Overwriting: {}, {}".format(lane_id, sample_id))

                    source_path_per_sample[sample_id].add(filename)

                    sample_d      = self.samples[sample_id]
                    lane_sample_d = self.sample_by_lane[lane_id][sample_id]

                    sample_d.reads      += sample_data["NumberReads"]
                    sample_d.bases      += sample_data["Yield"]
                    lane_sample_d.reads += sample_d.reads
                    lane_sample_d.bases += sample_d.bases

                    for index_metric in sample_data.get("IndexMetrics", []):
                        lane_d.perfect_index        += index_metric["MismatchCounts"]["0"]
                        sample_d.perfect_index      += index_metric["MismatchCounts"]["0"]
                        lane_sample_d.perfect_index += index_metric["MismatchCounts"]["0"]

                    for read_metric in sample_data.get("ReadMetrics", []):
                        r = read_metric["ReadNumber"]
                        lane_d.yield_q30                                   += read_metric["YieldQ30"]
                        lane_d.qscore_sum                                  += read_metric["QualityScoreSum"]
                        sample_d.yield_q30                                 += read_metric["YieldQ30"]
                        sample_d.qscore_sum                                += read_metric["QualityScoreSum"]
                        sample_d.read_stats["R{}_yield".format(r)]         += read_metric["Yield"]
                        sample_d.read_stats["R{}_Q30".format(r)]           += read_metric["YieldQ30"]
                        sample_d.read_stats["R{}_trimmed_bases".format(r)] += read_metric["TrimmedBases"]
                    # Remove unpopulated read keys
                    for r in range(1, 5):
                        if not sample_d.read_stats["R{}_yield".format(r)] and \
                           not sample_d.read_stats["R{}_Q30".format(r)] and \
                           not sample_d.read_stats["R{}_trimmed_bases".format(r)]:
                            sample_d.read_stats.pop("R{}_yield".format(r))
                            sample_d.read_stats.pop("R{}_Q30".format(r))
                            sample_d.read_stats.pop("R{}_trimmed_bases".format(r))

                # Add undetermined barcodes and reads, if we only process on bcl2fastq run.
                # if we have multiple runs however, we need to recalculate undetermined from total stats
                if len(found_files) == 1:
                    try:
                        unknown_barcodes = content['UnknownBarcodes'][lane_number - 1]['Barcodes']
                    except IndexError:
                        unknown_barcodes = next(
                            (item['Barcodes'] for item in content['UnknownBarcodes'] if item['Lane'] == 8),
                            None
                        )
                    lane_d.unknown_barcodes = unknown_barcodes

                    if "Undetermined" in lane_data:
                        undet = self.undetermined_by_lane[lane_id]
                        undet.reads       = lane_data["Undetermined"]["NumberReads"]
                        undet.bases = lane_data["Undetermined"]["Yield"]
                        for read_metric in lane_data["Undetermined"]["ReadMetrics"]:
                            undet.yield_q30     += read_metric["YieldQ30"]
                            undet.qscore_sum    += read_metric["QualityScoreSum"]
                            undet.trimmed_bases += read_metric["TrimmedBases"]

                        self.undetermined.reads         += undet.reads
                        self.undetermined.bases         += undet.bases
                        self.undetermined.yield_q30     += undet.yield_q30
                        self.undetermined.qscore_sum    += undet.qscore_sum
                        self.undetermined.trimmed_bases += undet.trimmed_bases
                else:
                    self.undetermined_by_lane[lane_id].reads = None
                    self.undetermined_by_lane[lane_id].bases = None

        return source_path_per_sample

    def _recalculate_undetermined(self):
        for lane_id in self.lanes:
            determined = Metrics()
            determined.reads = sum([s.reads for s in self.sample_by_lane[lane_id].values()])
            determined.bases = sum([s.bases for s in self.sample_by_lane[lane_id].values()])
            lane = self.lanes[lane_id]
            undet_reads = lane.reads - determined.reads
            undet_bases = lane.bases - determined.bases
            undet = self.undetermined_by_lane[lane_id]
            if undet.reads is not None:
                assert undet.reads == undet_reads
                assert undet.bases == undet_bases
            else:
                undet.reads = undet_reads
                undet.bases = undet_bases
                self.undetermined.reads += undet_reads
                self.undetermined.bases += undet_bases

    @staticmethod
    def get_unknown_barcodes(lane_unknown_barcode):
        """ Python 2.* dictionaries are not sorted.
        This function return an `OrderedDict` sorted by barcode count.
        """
        try:
            sorted_barcodes = OrderedDict(
                sorted(
                    lane_unknown_barcode.items(),
                    key=operator.itemgetter(1),
                    reverse=True
                )
            )
        except AttributeError:
            sorted_barcodes = None
        return sorted_barcodes

    def _make_samples_table(self):
        data = dict()
        for sample_id, d in self.samples.items():

            data[sample_id] = d.generate_metics_dict(self.total_stats)

            pct_q30_per_read = dict()
            for r in range(1, 5):
                # Zero division is possible
                try:
                    pct_q30_per_read[r] = '{0:.1f}'.format(float(
                        100.0 * d.read_stats["R{}_Q30".format(r)] / d.read_stats["R{}_yield".format(r)]))
                except ZeroDivisionError:
                    pct_q30_per_read[r] = '0.0'
                except KeyError:
                    pass

                try:
                    data[sample_id]["percent_R{}_Q30".format(r)] = pct_q30_per_read[r]
                    data[sample_id]["R{}_trimmed_bases".format(r)] = d.read_stats["R{}_trimmed_bases".format(r)]
                except KeyError:
                    pass

        return data

    def _make_lanes_table(self):
        return {
            lane_id: d.generate_metics_dict(self.total_stats)
            for lane_id, d in self.lanes.items()
        }

    def add_general_stats(self, data):
        headers = OrderedDict()
        headers['depth'] = {
            'title': 'Est. depth'.format(config.read_count_prefix),
            'description': 'Estimated depth based on the number of bases with quality score greater or equal to Q30, '
                           'assuming human whole genome',
                           # ' (or exome, if the sample id contains "rna" substring)',
            'min': 0,
            'suffix': 'X',
            'scale': 'BuPu'
        }
        headers['total'] = {
            'title': 'Clusters'.format(config.read_count_prefix),
            'description': 'Total number of clusters (read pairs) for this sample as determined by bcl2fastq demultiplexing ({})'.format(config.read_count_desc),
            'scale': 'Blues',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count',
            'format': read_format,
        }
        headers['clusters_pct'] = {
            'title': 'Clst %',
            'description': 'Percentage of clusters (read pairs) for this sample in the run, '
                           'as determined by bcl2fastq demultiplexing',
            'scale': 'RdYlGn',
            'max': 100,
            'min': 0,
            'suffix': '%'
        }
        headers['total_yield'] = {
            'title': 'Yield ({})'.format(config.base_count_prefix),
            'description': 'Total number of sequenced bases for this sample ({})'.format(config.base_count_desc),
            'scale': 'Greens',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count',
            'format': '{:,.1f}',
        }
        headers['yield_pct'] = {
            'title': 'Yld %',
            'description': 'Percentage of sequenced bases for this sample in this run',
            'scale': 'RdYlGn',
            'max': 100,
            'min': 0,
            'suffix': '%'
        }
        headers['yieldQ30'] = {
            'title': 'Yld &ge;Q30 (Mb)'.format(config.base_count_prefix),
            'description': 'Number of bases with a Phred score of 30 or higher ({})'.format(config.base_count_desc),
            'scale': 'Greens',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count',
            'format': '{:,.1f}',
            'hidden': True,
        }
        headers['percent_Q30'] = {
            'title': 'Yld &ge;Q30',
            'description': 'Percentage of bases with a Phred score of 30 or higher',
            'scale': 'RdYlGn',
            'max': 100,
            'min': 0,
            'suffix': '%',
        }
        # If no data for a column, header will be automatically ignored
        for r in range(1, 5):
            hide_col = True
            for s in data:
                try:
                    if float(data[s]["percent_R{}_Q30".format(r)]) > 0:
                        hide_col = False
                except KeyError:
                    pass
            headers['percent_R{}_Q30'.format(r)] = {
                'title': 'R{} &ge;Q30'.format(r),
                'description': 'Percent of bases in R{} with a Phred score of 30 or higher'.format(r),
                'scale': 'RdYlGn',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'hidden': hide_col
            }
        headers['percent_perfectIndex'] = {
            'title': 'Perf. idx',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%'
        }
        # If no data for a column, header will be automatically ignored
        for r in range(1, 5):
            hide_col = True
            for s in data:
                try:
                    if data[s]["R{}_trimmed_bases".format(r)] > 0:
                        hide_col = False
                except KeyError:
                    pass
            try:
                headers['R{}_trimmed_bases'.format(r)] = {
                    'title': 'R{} trimmed'.format(r),
                    'description': 'Number of bases trimmed ({})'.format(config.base_count_desc),
                    'scale': 'RdYlBu',
                    'modify': lambda x: x * 0.000001,
                    'hidden': hide_col,
                    'format': '{:,.1f} ' + config.base_count_prefix,
                }
            except KeyError:
                pass
        self.general_stats_addcols(data, headers)
        return data

    @staticmethod
    def lane_stats_table(data):
        """ Return a table with overview stats for each bcl2fastq lane for a single flow cell """
        headers = OrderedDict()
        headers['depth'] = {
            'title': 'Est. depth'.format(config.read_count_prefix),
            'description': 'Estimated depth based on the number of bases with quality score greater or equal to Q30',
            'min': 0,
            'suffix': 'X',
            'scale': 'BuPu'
        }
        headers['total'] = {
            'title': 'Clusters'.format(config.read_count_prefix),
            'description': 'Total number of clusters (read pairs) for this lane ({})'.format(config.read_count_desc),
            'scale': 'Blues',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count',
            'format': read_format,
        }
        headers['total_yield'] = {
            'title': 'Total yield'.format(config.base_count_prefix),
            'description': 'Tota number of bases yielded for this lane ({})'.format(config.base_count_desc),
            'scale': 'Greens',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count',
            'format': '{:,.1f}&nbsp;' + config.base_count_prefix,
        }
        headers['percent_Q30'] = {
            'title': 'Yield &ge;Q30',
            'description': 'Percentage of bases with greater than or equal to Q30 quality score',
            'suffix': '%',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn'
        }
        headers['mean_qscore'] = {
            'title': 'Mean quality',
            'description': 'Average phred qualty score',
            'min': 0,
            'scale': 'Spectral'
        }
        headers['percent_perfectIndex'] = {
            'title': 'Reads with perfect idx',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%'
        }
        table_config = {
            'namespace': 'bcl2fastq',
            'id': 'bcl2fastq-lane-stats-table',
            'table_title': 'bcl2fastq Lane Statistics',
            'col1_header': 'Run ID - Lane',
            'no_beeswarm': True
        }
        return table.plot(data, headers, table_config)

    @staticmethod
    def prepend_runid(run_id, rest):
        return str(run_id) + " - " + str(rest)

    @staticmethod
    def get_counts_bar_data_per_category(lanes_or_samples, undetermined_by_each=None):
        bar_data = {}
        for id, stats in lanes_or_samples.items():
            if stats.perfect_index is not None:
                bar_data[id] = {
                    "perfect": stats.perfect_index,
                    "imperfect": stats.reads - stats.perfect_index,
                }
                if undetermined_by_each:
                    bar_data[id]["undetermined"] = undetermined_by_each[id].reads
        return bar_data

    @staticmethod
    def get_perfect_counts_bar_data_per_lane(sample_by_lane):
        bar_data = defaultdict(dict)
        for lid, samples in sample_by_lane.items():
            for sid, stats in samples.items():
                if stats.perfect_index is not None:
                    bar_data[sid][lid] = stats.perfect_index
        return bar_data

    @staticmethod
    def get_bar_data_from_undetermined(lanes):
        """ Get data to plot for undetermined barcodes.
        """
        bar_data = defaultdict(dict)
        # get undetermined barcodes for each lanes
        for lid, lane in lanes.items():
            try:
                for barcode, count in islice(lane.unknown_barcodes.items(), 20):
                    bar_data[barcode][lane] = count
            except AttributeError:
                pass

        # sort results
        bar_data = OrderedDict(sorted(
            bar_data.items(),
            key=lambda x: sum(x[1].values()),
            reverse=True
        ))
        return OrderedDict(
            (key, value) for key, value in islice(bar_data.items(), 20)
        )

