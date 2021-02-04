# ChIP-Seq analysis pipeline (ChAP)

## Introduction

This ChIP-Seq analsysis pipeline (ChAP) will allow for step-by-step analysis of your paired-end ChIP-Seq NGS data. 

## Software dependencies

* [hisat2](http://daehwankimlab.github.io/hisat2/)
* [bwa](http://bio-bwa.sourceforge.net/)
* [samtools](http://www.htslib.org/)
* [bedtools](bedtools.readthedocs.io)
* [Picard](https://broadinstitute.github.io/picard/)
* [deeptools](https://deeptools.readthedocs.io)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info/)
* [MACS2](https://pypi.org/project/MACS2/)
* [HOMER](http://homer.ucsd.edu/homer/ngs/index.html)
* [ngs.plot](https://github.com/shenlab-sinai/ngsplot)
* [R](https://www.r-project.org/)
	* [tidyverse](https://www.tidyverse.org/)

## Installation

To download the repository (anywhere in /home):
> `git clone https://github.com/niekwit/ChIP-Seq-pipeline.git`

All the dependencies are expected to be in the environment variables, except for Picard (this can be installed anywhere in /home)

## Usage

Create a folder with any name that contains the sub-folder `raw-data` at any location. This sub-folder contains all the fastq.gz files. From the main folder run:
> `path/to/chip-seq-pipeline.sh [ fastqc ] [ align <aligner-pe/se> <species> ] [ dedup ] [ downsample ] [ qc ] [ bigwig ] [ peaks ] [ ngsplot ]`

It is recommended that the analysis is performed in a step-wise manner.
For example:
1. `./chip-seq-pipeline.sh fastqc`: check quality of fastq.gz files
2. `./chip-seq-pipeline.sh align bwa-pe human`: align files to human genome (no deduplication performed)
3. One can use a genome browser to map the BAM files and decide whether deduplication is required
4. `./chip-seq-pipeline.sh dedup`: perform deduplication
5. `./chip-seq-pipeline.sh downsampling`: perform downsampling on input file if read number is higher than corresponding ChIP sample
6. `./chip-seq-pipeline.sh qc`: perform quality control of bam files and replicates
7. `./chip-seq-pipeline.sh bigwig`: create BigWig files
8. `./chip-seq-pipeline.sh peaks`: call/annotate peaks with MACS2/HOMER
9. `./chip-seq-pipeline.sh ngsplot`: generate metagene plots and heatmaps
