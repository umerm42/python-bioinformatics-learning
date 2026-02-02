# python-bioinformatics-learning
Daily Python + Bioinformatics learning projects


## Week 1 Summary

- Python basics and scripting
- FASTQ parsing and QC
- Metadata cleaning and validation
- RNA-seq count table processing with Pandas
- Mini RNA-seq preprocessing pipeline
- Pipeline logic and architecture

### Key concepts I can explain
- FASTQ QC validates raw data; analysis starts after counts/quantification
- Clean metadata is required to avoid wrong sample grouping and silent errors
- Pipelines are modular for debugging, reuse, and scalability
- Inconsistent sample IDs can silently break merges and comparisons
- Pandas enables fast, reliable table operations for omics matrices

**Goal:** Transition from script runner to pipeline designer

## Week 2 Day 8 – Pandas 

Learned to:
- Compute sample-wise library sizes for QC
- Filter genes based on expression consistency
- Reshape RNA-seq matrices from wide to long format
- Merge count tables with metadata
- Aggregate gene expression by biological condition