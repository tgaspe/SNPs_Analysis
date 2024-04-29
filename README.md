# SNP Analysis Workflow
## Overview
This Snakemake workflow is designed to automate the process of SNP analysis from raw sequencing data. It covers several key steps: FastQC for quality control, BWA for alignment of reads to chromosome 21, SNP calling using bcftools, SNP normalization & filtering, SNP annotation using SNPeff, and extraction of specific SNPs associating with particular genes.

## Workflow Steps
- FastQC: Quality control checks on raw sequence data.
- BWA MEM: Aligns reads to a reference genome.
- bcftools: Calls SNPs from the alignment.
- SNP normalization & filtering: Normalizes and filters SNP calls.
- SNP annotation (SNPeff): Annotates SNPs with genomic information.
- Extracting specific SNPs: Extracts SNPs associated with the genes APP, SOD1, and DYRK1A.
## Requirements
- BWA
- Samtools
- bcftools
- vt (variant tools for normalization)
- SNPeff
- Python 3.x
- Matplotlib (for plotting graphs)
