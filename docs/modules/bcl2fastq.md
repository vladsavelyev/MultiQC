---
Name: bcl2fastq
URL: https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html
Description: >
  bcl2fastq can be used to both demultiplex data and convert BCL files to
  FASTQ file formats for downstream analysis.
---

There are two versions of this software: bcl2fastq for MiSeq and HiSeq
sequencing systems running RTA versions earlier than 1.8, and bcl2fastq2 for
Illumina sequencing systems running RTA version 1.18.54 and above. This module
currently only covers output from the latter.

#### Calculate estimated depth

You can specify a genome size in config

It's often useful to talk about sequencing yield in terms of estimated depth of coverage.
In order to make MultiQC show the estimated depth for each sample, specify the reference genome/target size in your [MultiQC configuration](http://multiqc.info/docs/#configuring-multiqc):

```yaml
bcl2fastq:
    genome_size: 3049315783
```

The coverage depth will be estimated as the yield Q30 dvivided by the genome size.

MultiQC comes with effective genome size presets for Human and Mouse, so you can
provide the genome build name instead, like this: `genome_size: hg38_genome`. The
following values are supported: `hg19_genome`, `hg38_genome`, `mm10_genome`.
