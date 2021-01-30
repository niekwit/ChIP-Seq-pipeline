#!/bin/bash

#This pipeline performs ChIP-Seq analsyis for paired-end data (Niek Wit, University of Cambridge, 2021)
SCRIPT_DIR="/home/niek/Documents/scripts/ChIP-Seq-pipeline/"
PICARD="/home/niek/Documents/scripts/Picard/picard.jar"
max_threads=40 #$(nproc --all) #determines CPU thread count

usage() {                                    
	echo "Usage: $0 [ rename ] [ fastqc ] [ align <species> ] [ dedup ] [ pca ] [ bigwig ] [ peaks ] [ ngsplot ]"
	echo "rename: renames fastq files from configuration file (rename.config)"
	echo "fastqc: performs FastQC and MultiQC on fq.gz files"
	echo "align: aligns fq.gz files to selected genome (no deduplication) using HISAT2"
	echo -e "Available genomes for alignment:\n\t<human>: hg38\n\t<mouse>: mm9"
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
	source rename.sh
	rename
fi

if [[ $@ == *"fastqc"* ]];
then
	source fastqc.sh
	fastqc
fi

if [[ $@ == *"align"* ]];
then
    	source align.sh
	align
fi

if [[ $@ == *"dedup"* ]] || [[ $@ == *"deduplication"* ]];
then
    	source dedup.sh
	dedup
fi

if [[ $@ == *"bigwig"* ]];
then
    	source bigwig.sh
	bigwig
fi

if [[ $@ == *"peak"* ]];
then
    	source peaks.sh
	peaks
fi

if [[ $@ == *"ngsplot"* ]];
then
    	source ngsplot.sh
	ngsplot
fi

if [[ $@ == *"qc"* ]] || [[ $@ == *"QC"* ]];
then
    	source qc.sh
	qc
fi
