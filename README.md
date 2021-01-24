# ChIP-Seq analysis pipeline (ChAP)

## Introduction

This ChIP-Seq analsysis pipeline will allow for step-by-step analysis of your paired-end ChIP-Seq NGS data. 

## Software dependencies

* hisat2
* samtools
* bedtools
* Picard
* deeptools
* FastQC
* MultiQC
* MACS2

## Usage

Usage: `./chip-seq-pipeline.sh [ fastqc ] [ align <species> ] [ dedup ] [ pca ] [ bigwig ]`

Run `chip-seq-pipeline.sh` from a folder that contains the subfolder `raw-data`. This folder must contain all the fastq.gz files.
It is recommended that the analysis is performed in a step-wise manner.
For example:
1. `./chip-seq-pipeline.sh fastqc`: check quality of fastq.gz files
2. `./chip-seq-pipeline.sh align human`: align files to human genome (no deduplication performed)
3. One can use SeqMonk to map the BAM files and decide whether deduplication is required
4. `./chip-seq-pipeline.sh dedup`: perform deduplication
5. `./chip-seq-pipeline.sh pca`: perform principle component analysis (PCA)
6. `./chip-seq-pipeline.sh bigwig`: create BigWig files

