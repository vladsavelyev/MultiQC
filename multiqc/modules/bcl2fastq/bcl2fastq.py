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


class MultiqcModule(BaseMultiqcModule):
    sample_metrics = [
        "total",
        "total_yield",
        "perfectIndex",
        "yieldQ30",
        "qscore_sum",
    ]

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='bcl2fastq',
            anchor='bcl2fastq',
            href="https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html",
            info="can be used to both demultiplex data and convert BCL files"
                 " to FASTQ file formats for downstream analysis."
        )

        # Gather data from all json files

        bcl2fastq_datas = []
        found_files = [f for f in self.find_log_files('bcl2fastq')]
        for myfile in found_files:
            bcl2fastq_datas.append(self.parse_file_as_json(myfile))

        if len(bcl2fastq_datas) == 0:
            raise UserWarning

        # Make sure all data belong to one sequencer run
        datas_by_sequencerrun = defaultdict(list)
        for bcl2fastq_data in bcl2fastq_datas:
            for sequencer_run_id, sequencer_run_data in bcl2fastq_data.items():
                datas_by_sequencerrun[sequencer_run_id].append(sequencer_run_data)
        if len(datas_by_sequencerrun) > 1:
            log.error("Multiple runs found in files", found_files)
            raise UserWarning
        sequencer_run_id, bcl2fastq_runs = list(datas_by_sequencerrun.items())[0]

        if len(bcl2fastq_datas) > 1:
            log.warning("Warning: multiple bcl2fastq runs detected. "
                        "They will be merged, undetermined information will be dropped.")

        # Collect counts by lane and sample (+source_files)
        bcl2fastq_bylane, bcl2fastq_bysample, bcl2fastq_bysample_lane, source_files \
            = MultiqcModule.split_data_by_lane_and_sample(sequencer_run_id, bcl2fastq_runs)

        # Filter to strip out ignored sample names
        bcl2fastq_bylane = self.ignore_samples(bcl2fastq_bylane)
        bcl2fastq_bysample = self.ignore_samples(bcl2fastq_bysample)
        bcl2fastq_bysample_lane = self.ignore_samples(bcl2fastq_bysample_lane)

        # Return with Warning if no files are found
        if len(bcl2fastq_bylane) == 0 and len(bcl2fastq_bysample) == 0:
            raise UserWarning

        # Print source files
        for s in source_files.keys():
            self.add_data_source(
                s_name=s,
                source=",".join(list(set(source_files[s]))),
                module='bcl2fastq',
                section='bcl2fastq-bysample'
            )

        # Add sample counts to general stats table
        self.add_general_stats(bcl2fastq_bysample)
        self.write_data_file(
            {str(k): bcl2fastq_bylane[k] for k in bcl2fastq_bylane.keys()},
            'multiqc_bcl2fastq_bylane'
        )
        self.write_data_file(bcl2fastq_bysample, 'multiqc_bcl2fastq_bysample')

        # Add section for summary stats per flow cell
        self.add_section (
            name = 'Lane Statistics',
            anchor = 'bcl2fastq-lanestats',
            description = 'Statistics about each lane for each flowcell',
            plot = self.lane_stats_table(bcl2fastq_bylane)
        )

        # Add section for counts by lane
        cats = OrderedDict()
        cats["perfect"] = {'name': 'Perfect Index Reads'}
        cats["imperfect"] = {'name': 'Mismatched Index Reads'}
        cats["undetermined"] = {'name': 'Undetermined Reads'}
        self.add_section (
            name = 'Clusters by lane',
            anchor = 'bcl2fastq-bylane',
            description = 'Number of reads per lane (with number of perfect index reads).',
            helptext = """Perfect index reads are those that do not have a single mismatch.
                All samples of a lane are combined. Undetermined reads are treated as a third category.""",
            plot = bargraph.plot(
                self.get_counts_bar_data_per_category(bcl2fastq_bylane),
                cats,
                {
                    'id': 'bcl2fastq_lane_counts',
                    'title': 'bcl2fastq: Clusters by lane',
                    'ylab': 'Number of clusters',
                    'hide_zero_cats': False
                }
            )
        )

        # Add section for counts by sample
        # get cats for per-lane tab
        lcats = set()
        for s_name in bcl2fastq_bysample_lane:
            lcats.update(bcl2fastq_bysample_lane[s_name].keys())
        lcats = sorted(list(lcats))
        self.add_section (
            name = 'Clusters by sample',
            anchor = 'bcl2fastq-bysample',
            description = 'Number of reads per sample.',
            helptext = """Perfect index reads are those that do not have a single mismatch.
                All samples are aggregated across lanes combinned. Undetermined reads are ignored.
                Undetermined reads are treated as a separate sample.""",
            plot = bargraph.plot(
                [
                    MultiqcModule.get_counts_bar_data_per_category(bcl2fastq_bysample),
                    MultiqcModule.get_perfect_counts_bar_data_per_lane(bcl2fastq_bysample_lane)
                ],
                [cats, lcats],
                {
                    'id': 'bcl2fastq_sample_counts',
                    'title': 'bcl2fastq: Clusters by sample',
                    'hide_zero_cats': False,
                    'ylab': 'Number of clusters',
                    'data_labels': ['Index mismatches', 'Counts per lane']
                }
            )
        )

        if len(bcl2fastq_datas) == 1:
            # Add section with undetermined barcodes
            self.add_section(
                name = "Undetermined barcodes by lane",
                anchor = "undetermine_by_lane",
                description = "Count of the top twenty most abundant undetermined barcodes by lanes",
                plot = bargraph.plot(
                    MultiqcModule.get_bar_data_from_undetermined(bcl2fastq_bylane),
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

    @staticmethod
    def parse_file_as_json(myfile):
        bcl2fastq_data = dict()

        try:
            content = json.loads(myfile["f"])
        except ValueError:
            log.warning('Could not parse file as json: {}'.format(myfile["fn"]))
            return
        runId = content["RunId"]
        if runId not in bcl2fastq_data:
            bcl2fastq_data[runId] = dict()
        run_data = bcl2fastq_data[runId]
        for conversionResult in content.get("ConversionResults", []):
            l = conversionResult["LaneNumber"]
            lane = 'L{}'.format(conversionResult["LaneNumber"])
            if lane in run_data:
                log.debug("Duplicate runId/lane combination found! Overwriting: {}".format(MultiqcModule.prepend_runid(runId, lane)))
            run_data[lane] = {
                "total": 0,
                "total_yield": 0,
                "perfectIndex": 0,
                "samples": dict(),
                "yieldQ30": 0,
                "qscore_sum": 0
            }
            # Add undetermined barcodes
            unknown_barcode = dict()
            for lane_data in content.get("UnknownBarcodes", list()):
                if lane_data["Lane"] == l:
                    unknown_barcode = lane_data["Barcodes"]
                    break
            run_data[lane]["unknown_barcodes"] = unknown_barcode

            for demuxResult in conversionResult.get("DemuxResults", []):
                if demuxResult["SampleName"] == demuxResult["SampleId"]:
                    sample = demuxResult["SampleName"]
                else:
                    sample = "{}-{}".format(demuxResult["SampleId"], demuxResult["SampleName"])
                if sample in run_data[lane]["samples"]:
                    log.debug("Duplicate runId/lane/sample combination found! Overwriting: {}, {}".format(MultiqcModule.prepend_runid(runId, lane), sample))
                run_data[lane]["samples"][sample] = {
                    "total": 0,
                    "total_yield": 0,
                    "perfectIndex": 0,
                    "filename": os.path.join(myfile['root'], myfile["fn"]),
                    "yieldQ30": 0,
                    "qscore_sum": 0,
                }
                for r in range(1,5):
                    run_data[lane]["samples"][sample]["R{}_yield".format(r)] = 0
                    run_data[lane]["samples"][sample]["R{}_Q30".format(r)] = 0
                    run_data[lane]["samples"][sample]["R{}_trimmed_bases".format(r)] = 0
                run_data[lane]["total"] += demuxResult["NumberReads"]
                run_data[lane]["total_yield"] += demuxResult["Yield"]
                run_data[lane]["samples"][sample]["total"] += demuxResult["NumberReads"]
                run_data[lane]["samples"][sample]["total_yield"] += demuxResult["Yield"]
                for indexMetric in demuxResult.get("IndexMetrics", []):
                    run_data[lane]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                    run_data[lane]["samples"][sample]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                for readMetric in demuxResult.get("ReadMetrics", []):
                    r = readMetric["ReadNumber"]
                    run_data[lane]["yieldQ30"] += readMetric["YieldQ30"]
                    # run_data[lane]["depth"] += run_data[lane]["yieldQ30"] / GENOME_SIZE
                    run_data[lane]["qscore_sum"] += readMetric["QualityScoreSum"]
                    run_data[lane]["samples"][sample]["yieldQ30"] += readMetric["YieldQ30"]
                    # run_data[lane]["samples"][sample]["depth"] += readMetric["YieldQ30"] / GENOME_SIZE
                    run_data[lane]["samples"][sample]["qscore_sum"] += readMetric["QualityScoreSum"]
                    run_data[lane]["samples"][sample]["R{}_yield".format(r)] += readMetric["Yield"]
                    run_data[lane]["samples"][sample]["R{}_Q30".format(r)] += readMetric["YieldQ30"]
                    run_data[lane]["samples"][sample]["R{}_trimmed_bases".format(r)] += readMetric["TrimmedBases"]
                # Remove unpopulated read keys
                for r in range(1,5):
                    if not run_data[lane]["samples"][sample]["R{}_yield".format(r)] and \
                       not run_data[lane]["samples"][sample]["R{}_Q30".format(r)] and \
                       not run_data[lane]["samples"][sample]["R{}_trimmed_bases".format(r)]:
                        run_data[lane]["samples"][sample].pop("R{}_yield".format(r))
                        run_data[lane]["samples"][sample].pop("R{}_Q30".format(r))
                        run_data[lane]["samples"][sample].pop("R{}_trimmed_bases".format(r))

            undeterminedYieldQ30 = 0
            undeterminedQscoreSum = 0
            undeterminedTrimmedBases = 0
            if "Undetermined" in conversionResult:
                for readMetric in conversionResult["Undetermined"]["ReadMetrics"]:
                    undeterminedYieldQ30 += readMetric["YieldQ30"]
                    undeterminedQscoreSum += readMetric["QualityScoreSum"]
                    undeterminedTrimmedBases += readMetric["TrimmedBases"]
                run_data[lane]["samples"]["undetermined"] = {
                    "total": conversionResult["Undetermined"]["NumberReads"],
                    "total_yield": conversionResult["Undetermined"]["Yield"],
                    "perfectIndex": 0,
                    "yieldQ30": undeterminedYieldQ30,
                    "qscore_sum": undeterminedQscoreSum,
                    "trimmed_bases": undeterminedTrimmedBases
                }

        # Calculate Percents and averages
        for lane_id, lane in run_data.items():
            try:
                lane["percent_Q30"] = (float(lane["yieldQ30"])
                    / float(lane["total_yield"])) * 100.0
            except ZeroDivisionError:
                lane["percent_Q30"] = "NA"
            try:
                lane["percent_perfectIndex"] = (float(lane["perfectIndex"])
                    / float(lane["total"])) * 100.0
            except ZeroDivisionError:
                lane["percent_perfectIndex"] = "NA"
            try:
                lane["mean_qscore"] = float(lane["qscore_sum"]) / float(lane["total_yield"])
            except ZeroDivisionError:
                lane["mean_qscore"] = "NA"
            for sample_id, sample in lane["samples"].items():
                try:
                    sample["percent_Q30"] = (float(sample["yieldQ30"]) / float(sample["total_yield"])) * 100.0
                except ZeroDivisionError:
                    sample["percent_Q30"] = "NA"
                try:
                    sample["percent_perfectIndex"] = (float(sample["perfectIndex"])
                        / float(sample["total"])) * 100.0
                except ZeroDivisionError:
                    sample["percent_perfectIndex"] = "NA"
                try:
                    sample["mean_qscore"] = float(sample["qscore_sum"]) / float(sample["total_yield"])
                except ZeroDivisionError:
                    sample["mean_qscore"] = "NA"

        return bcl2fastq_data

    @staticmethod
    def _process_run_lane_data(lane_data, uniq_lane_id, is_single_run=True):
        processed_lane_data = defaultdict(int)
        bysample_lane_data = defaultdict(dict)
        source_files = dict()
        processed_run_data = {  # sums up everything, including undetermined
            m: sum([raw_s[m] for raw_s in lane_data["samples"].values()])
            for m in [
                "total",
                "total_yield",
            ]
        }

        # all_lane_stats = bcl2fastq_bylane[uniq_lane_id]
        for m in [
                "total",
                "total_yield",
                "perfectIndex",
                "yieldQ30",
                "qscore_sum",
                "percent_Q30",
                "percent_perfectIndex",
                "mean_qscore",
                "unknown_barcodes",
            ]:
            val = lane_data[m]
            if isinstance(val, int) or isinstance(val, float):
                processed_lane_data[m] += val

        if is_single_run == 1:
            processed_lane_data["undetermined"] = lane_data["samples"].get("undetermined", {}).get("total", "NA")
            processed_lane_data["unknown_barcodes"] = MultiqcModule.get_unknown_barcodes(lane_data['unknown_barcodes'])

        ################
        ### Building per-sample-lane stats
        for sample_id, raw_s in lane_data["samples"].items():
            s = dict()
            for m in MultiqcModule.sample_metrics:
                s[m] = raw_s[m]

            # Undetermined samples did not have R1 and R2 information
            for r in range(1, 5):
                for m in ["yield", "Q30", "trimmed_bases"]:
                    try:
                        s["R{}_{}".format(r, m)] = raw_s["R{}_{}".format(r, m)]
                    except KeyError:
                        pass

            if sample_id != "undetermined":
                total_size = GENOME_SIZE
                # if 'rna' in sample_id:
                #     total_size = EXOME_SIZE
                s["depth"] = s["yieldQ30"] / total_size

            if sample_id != "undetermined":
                if sample_id not in source_files:
                    source_files[sample_id] = []
                source_files[sample_id].append(raw_s["filename"])

            bysample_lane_data[sample_id][uniq_lane_id] = s

        return processed_run_data, processed_lane_data, bysample_lane_data, source_files

    @staticmethod
    def _process_sample_data(sample_data_by_lane):
        sample_data = defaultdict(int)
        for lane_id, sample_lane_data in sample_data_by_lane.items():
            for m, val in sample_lane_data.items():
                if isinstance(val, int) or isinstance(val, float):
                    sample_data[m] += val

            s = sample_data  # to populate
            try:
                s["percent_Q30"] = (float(sample_lane_data["yieldQ30"]) / float(s["total_yield"])) * 100.0
            except ZeroDivisionError:
                s["percent_Q30"] = "NA"
            try:
                s["percent_perfectIndex"] = (float(s["perfectIndex"]) / float(s["total"])) * 100.0
            except ZeroDivisionError:
                s["percent_perfectIndex"] = "NA"
            try:
                s["mean_qscore"] = (float(s["qscore_sum"]) / float(s["total_yield"]))
            except ZeroDivisionError:
                s["mean_qscore"] = "NA"

            # Remove unpopulated read keys
            for r in range(1, 5):
                try:
                    if not s["R{}_yield".format(r)] and \
                       not s["R{}_Q30".format(r)] and \
                       not s["R{}_trimmed_bases".format(r)]:
                        s.pop("R{}_yield".format(r))
                        s.pop("R{}_Q30".format(r))
                        s.pop("R{}_trimmed_bases".format(r))
                except KeyError:
                    pass
        return sample_data

    @staticmethod
    def split_data_by_lane_and_sample(sequencer_run_id, bcl2fastq_runs):
        data_bylane = defaultdict(lambda: defaultdict(int))
        data_bysample_bylane = defaultdict(dict)
        data_bysample = dict()
        source_files = dict()
        total_run_data = defaultdict(int)  # based on single bcl2fastq run (we assume all runs belong to the same sequencer run)

        for bcl2fastqw_run_i, bcl2fastq_run_data in enumerate(bcl2fastq_runs):
            for lane_id, lane_data in bcl2fastq_run_data.items():
                unique_lane_id = MultiqcModule.prepend_runid(sequencer_run_id, lane_id)
                processed_run_lane_data, processed_lane_data, new_bysample_lane_data, new_source_files = \
                    MultiqcModule._process_run_lane_data(lane_data, unique_lane_id,
                                                         is_single_run=len(bcl2fastq_runs) == 1)

                if bcl2fastqw_run_i == 0:
                    for m, val in processed_run_lane_data.items():
                        total_run_data[m] += val

                data_bylane[unique_lane_id] = processed_lane_data
                for s, ld in new_bysample_lane_data.items():
                    data_bysample_bylane[s].update(ld)
                source_files.update(new_source_files)

            ################
            ### Adding up to the sample-level stats
            for sample_id, sample_data_bylane in data_bysample_bylane.items():
                if sample_id == 'undetermined' and len(bcl2fastq_runs) > 1:
                    # Undetermined stats for multiple bcl2fastq runs are off as they don't consider other runs
                    # recalculating them, assuming that all bcl2fastq runs make up to a full sequencer run
                    continue

                data_bysample[sample_id] = MultiqcModule._process_sample_data(sample_data_bylane)

        # #############
        # ### Summarising total run data to properly calculate undetermined and percentages
        determined = {
            m: sum(data[m] for sid, data in data_bysample.items() if sid != 'undetermined')
            for m in total_run_data.keys()
        }
        undetermined = {m: total_run_data[m] - determined[m] for m in determined.keys()}
        if len(bcl2fastq_runs) > 1:
            data_bysample['undetermined'] = undetermined
        # else:
        #     assert data_bysample['undetermined']['total'] == undetermined['total'], (data_bysample['undetermined'], undetermined)

        for sample_id, sample_data in data_bysample.items():
            sample_data['clusters_pct'] = sample_data['total'] / total_run_data['total'] * 100.0
            sample_data['yield_pct'] = sample_data['total_yield'] / total_run_data['total_yield'] * 100.0

        return data_bylane, data_bysample, data_bysample_bylane, source_files

    @staticmethod
    def get_unknown_barcodes(lane_unknown_barcode):
        """ Python 2.* dictionaries are not sorted.
        This function return an `OrderedDict` sorted by barcode count.
        """
        try:
            sorted_barcodes = OrderedDict(
                sorted(
                    lane_unknown_barcode.items(),
                    key = operator.itemgetter(1),
                    reverse = True
                )
            )
        except AttributeError:
            sorted_barcodes = None
        return sorted_barcodes

    def add_general_stats(self, bcl2fastq_bysample):
        data = dict()
        for sample_id, sample in bcl2fastq_bysample.items():
            percent_R_Q30 = dict()
            for r in range(1,5):
                # Zero division is possible
                try:
                    percent_R_Q30[r] = '{0:.1f}'.format(float(100.0 * sample["R{}_Q30".format(r)] / sample["R{}_yield".format(r)]))
                except ZeroDivisionError:
                    percent_R_Q30[r] = '0.0'
                except KeyError:
                    pass

            data[sample_id] = {
                'total': sample['total'],
                'total_yield': sample['total_yield'],
                'clusters_pct': sample['clusters_pct'],
                'yield_pct': sample['yield_pct'],
            }

            if 'perfectIndex' in sample:
                try:
                    perfect_percent = '{0:.1f}'.format(float(100.0 * sample["perfectIndex"] / sample["total"]))
                except ZeroDivisionError:
                    perfect_percent = '0.0'
                data[sample_id]['perfectIndex'] = perfect_percent

            if 'yieldQ30' in sample:
                data[sample_id]['yieldQ30'] = sample['yieldQ30']

            if sample_id != 'undetermined':
                data[sample_id]['depth'] = sample['depth']

            for r in range(1,5):
                try:
                    data[sample_id]["percent_R{}_Q30".format(r)] = percent_R_Q30[r]
                    data[sample_id]["R{}_trimmed_bases".format(r)] = sample["R{}_trimmed_bases".format(r)]
                except KeyError:
                    pass

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
            'title': 'Clust. %',
            'description': 'Percentage of clusters (read pairs) for this sample in the run, '
                           'as determined by bcl2fastq demultiplexing',
            'scale': 'RdYlGn',
            'max': 100,
            'min': 0,
            'suffix': '%'
        }
        headers['total_yield'] = {
            'title': 'Total Yield ({})'.format(config.base_count_prefix),
            'description': 'Total number of sequenced bases for this sample ({})'.format(config.base_count_desc),
            'scale': 'Greens',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count',
            'format': '{:,.1f}',
        }
        headers['yield_pct'] = {
            'title': 'Yield %',
            'description': 'Percentage of sequenced bases for this sample in this run',
            'scale': 'RdYlGn',
            'max': 100,
            'min': 0,
            'suffix': '%'
        }
        headers['yieldQ30'] = {
            'title': 'Yield &ge;Q30 (Mb)'.format(config.base_count_prefix),
            'description': 'Number of bases with a Phred score of 30 or higher ({})'.format(config.base_count_desc),
            'scale': 'Greens',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count',
            'format': '{:,.1f}',
        }
        # If no data for a column, header will be automatically ignored
        for r in range(1,5):
            hideCol = True
            for s in data:
                try:
                    if float(data[s]["percent_R{}_Q30".format(r)]) > 0:
                        hideCol = False
                except KeyError:
                    pass
            headers['percent_R{}_Q30'.format(r)] = {
                'title': 'R{} Yield &ge;Q30'.format(r),
                'description': 'Percent of bases in R{} with a Phred score of 30 or higher'.format(r),
                'scale': 'RdYlGn',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'hidden': hideCol
            }
        headers['perfectPercent'] = {
            'title': 'Perf Idx',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%'
        }
        # If no data for a column, header will be automatically ignored
        for r in range(1,5):
            hideCol = True
            for s in data:
                try:
                    if data[s]["R{}_trimmed_bases".format(r)] > 0:
                        hideCol = False
                except KeyError:
                    pass
            try:
                headers['R{}_trimmed_bases'.format(r)] = {
                    'title': 'R{} trimmed'.format(r),
                    'description': 'Number of bases trimmed ({})'.format(config.base_count_desc),
                    'scale': 'RdYlBu',
                    'modify': lambda x: x * 0.000001,
                    'hidden': hideCol,
                    'format': '{:,.1f} ' + config.base_count_prefix,
                }
            except KeyError:
                pass
        self.general_stats_addcols(data, headers)

    def lane_stats_table(self, bcl2fastq_bylane):
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
            'title': 'Total Yield'.format(config.base_count_prefix),
            'description': 'Tota number of bases yielded for this lane ({})'.format(config.base_count_desc),
            'scale': 'Greens',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count',
            'format': '{:,.1f}&nbsp;' + config.base_count_prefix,
        }
        headers['percent_Q30'] = {
            'title': 'Bases &ge; Q30',
            'description': 'Percentage of bases with greater than or equal to Q30 quality score',
            'suffix': '%',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn'
        }
        headers['mean_qscore'] = {
            'title': 'Mean Quality',
            'description': 'Average phred qualty score',
            'min': 0,
            'scale': 'Spectral'
        }
        headers['percent_perfectIndex'] = {
            'title': 'Perfect Index',
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
        return table.plot(bcl2fastq_bylane, headers, table_config)

    @staticmethod
    def prepend_runid(runId, rest):
        return str(runId)+" - "+str(rest)

    @staticmethod
    def get_counts_bar_data_per_category(bcl2fastq_bysample):
        bar_data = {}
        for sample_id, sample_data in bcl2fastq_bysample.items():
            if 'perfectIndex' in sample_data:
                bar_data[sample_id] = {
                    "perfect": sample_data["perfectIndex"],
                    "imperfect": sample_data["total"] - sample_data["perfectIndex"],
                }
                if "undetermined" in sample_data:
                    bar_data[sample_id]["undetermined"] = sample_data["undetermined"]
        return bar_data

    @staticmethod
    def get_perfect_counts_bar_data_per_lane(bcl2fastq_bysample_lane):
        bar_data = {}
        for sname, per_lane in bcl2fastq_bysample_lane.items():
            bar_data[sname] = {}
            for lane, stats in per_lane.items():
                if 'perfectIndex' in stats:
                    bar_data[sname][lane] = stats["perfectIndex"]
        return bar_data

    @staticmethod
    def get_bar_data_from_undetermined(bcl2fastq_bylane):
        """ Get data to plot for undetermined barcodes.
        """
        bar_data = defaultdict(dict)
        # get undetermined barcodes for each lanes
        for lane_id, lane in bcl2fastq_bylane.items():
            try:
                for barcode, count in islice(lane['unknown_barcodes'].items(), 20):
                    bar_data[barcode][lane_id] = count
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
