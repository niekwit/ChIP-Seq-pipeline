#!/bin/bash

#This pipeline performs ChIP-Seq analsyis for paired-end data (Niek Wit, University of Cambridge, 2021)
PICARD="/home/niek/Documents/scripts/Picard/picard.jar"
max_threads=$(nproc --all) #determines CPU thread count
SCRIPT_DIR="/home/niek/Documents/scripts/ChIP-Seq-scripts/"

usage() {                                    
	echo "Usage: $0 [ rename ] [ fastqc ] [ align <species> ] [ dedup ] [ qc ] [ bigwig ] [ peaks ] [ ngsplot ]"
	echo "rename: renames fastq files from configuration file (rename.config)"
	echo "fastqc: performs FastQC and MultiQC on fq.gz files"
	echo "align: aligns fq.gz files to selected genome (no deduplication) using HISAT2"
	echo -e "Available genomes for alignment:\n\t<human>: hg19\n\t<mouse>: mm9"
	echo "dedup: removes duplicates using PICARD"
	echo "qc: performs quality control on alignment files"
	echo "bigwig: creates bigWig files"
	echo "peaks: calls peaks with MACS2 using selected BAM files (enter samples in macs2-input.csv)"
	echo "ngsplot: generates metagene plots and heatmaps with ngsplot"
	exit 2
}

while getopts '?h' c
do
	case $c in
		h|?) usage;;
	esac
done

if [[ $@ == *"rename"* ]];
then
	source "${SCRIPT_DIR}rename.sh"
	rename
fi

if [[ $@ == *"fastqc"* ]];
then
	source "${SCRIPT_DIR}fastqc.sh"
	fastqc
fi

if [[ $@ == *"align"* ]];
then
    	source "${SCRIPT_DIR}align.sh"
	align
fi

if [[ $@ == *"dedup"* ]] || [[ $@ == *"deduplication"* ]];
then
    	source "${SCRIPT_DIR}dedup.sh"
	dedup
fi

if [[ $@ == *"bigwig"* ]];
then
    	source "${SCRIPT_DIR}bigwig.sh"
	bigwig
fi

if [[ $@ == *"peak"* ]];
then
    	source "${SCRIPT_DIR}peaks.sh"
	peaks
fi

if [[ $@ == *"ngsplot"* ]];
then
    	source "${SCRIPT_DIR}ngsplot.sh"
	ngsplot
fi

if [[ $@ == *"qc"* ]] || [[ $@ == *"QC"* ]];
then
    	source "${SCRIPT_DIR}qc.sh"
	qc
fi
