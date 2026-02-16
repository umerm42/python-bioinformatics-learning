# RNA-seq QC Pipeline (FastQC + fastp + MultiQC)

Config-driven, parallel per-sample QC pipeline for paired-end FASTQ (.fq.gz).

## What it does
Per sample:
- FastQC on raw reads
- fastp trimming (HTML/JSON)
- FastQC on trimmed reads

Once per run:
- MultiQC summary report
- QC_REPORT.md summary

## Requirements
- Python 3.9+
- PyYAML: `pip install pyyaml` or `conda install pyyaml`
- Tools in PATH:
  - fastqc
  - fastp
  - multiqc

## Input file: samples.tsv (tab-separated)
Required columns:

```tsv
sample	fq1	fq2
R001	/path/to/R001_R1.fq.gz	/path/to/R001_R2.fq.gz
