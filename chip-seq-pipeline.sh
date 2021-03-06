#!/bin/bash

#This pipeline performs ChIP-Seq analsyis for paired-end data (Niek Wit, University of Cambridge, 2021)

PICARD=$(find $HOME -name picard.jar)
SCRIPT_DIR=$(find $HOME -type d -name "ChIP-Seq-pipeline")
max_threads=$(nproc --all) #determines CPU thread count
#source "${SCRIPT_DIR}/genome-index.conf" #loads locations of genome indeces, etc
align_prog=""
dedup=""
bigwig=""
fastqc=""
downsample=""
qc=""

usage() {
	echo "Usage: $0 [rename] [fastqc] [align <species>] [dedup] [qc] [bigwig] [peaks] [ngsplot]"
	echo "rename: renames fastq files from configuration file (rename.config)"
	echo "fastqc: performs FastQC and MultiQC on fq.gz files"
	echo "align: aligns fq.gz files to selected genome (no deduplication) using bwa or HISAT2 in single-end or paired-end mode:"
	echo -e "\tAvailable alignment algorithms:\n\t\t<bwa-se>: bwa single-end mode\n\t\t<bwa-pe>: bwa paired-end mode\n\t\t<hisat-se>: HISAT2 single-end mode\n\t\t<hisat-pe>: HISAT2 paired-end mode"
	echo -e "\tAvailable genomes for alignment:\n\t\t<human>: hg19\n\t\t<mouse>: mm9"
	echo "dedup: removes duplicates using PICARD"
	echo "downsample: downsamples input files to corresponding ChIP file"
	echo "qc: performs quality control on alignment files"
	echo "bigwig: creates bigWig files"
	echo "peaks: calls peaks with MACS2 using selected BAM files (enter samples in macs2-input.csv)"
	echo "ngsplot: generates metagene plots and heatmaps with ngsplot"
	exit 2
}


while getopts 'qsfbda:g:?h' c
do
	case $c in
		a) align_prog=$OPTARG ;;
		b) bigwig="TRUE" ;;
		g) genome=$OPTARG ;;
		d) dedup="TRUE" ;;
		f) fastqc="TRUE";;
		s) down_sample="TRUE";;
		q) qc="TRUE";;
		h|?) usage ;;
	esac
done

if [[ $fastqc == "TRUE" ]]; then
	source "${SCRIPT_DIR}/fastqc.sh"
fi

if [[ -n $align_prog ]]; then
	source "${SCRIPT_DIR}/align.sh"
fi

if [[ $dedup == "TRUE" ]]; then
	source "${SCRIPT_DIR}/dedup.sh"
fi

if [[ $bigwig == "TRUE" ]]; then
	source "${SCRIPT_DIR}/bigwig.sh"
fi

if [[ $downsample == "TRUE" ]]; then
	source "${SCRIPT_DIR}/downsample.sh"
fi


: <<'END'
if [[ $@ == *"rename"* ]];
then
	source "${SCRIPT_DIR}/rename.sh"
fi

if [[ $@ == *"fastqc"* ]];
then
	source "${SCRIPT_DIR}/fastqc.sh"
fi

if [[ $@ == *"align"* ]];
then
    	source "${SCRIPT_DIR}/align.sh"
fi

if [[ $@ == *"dedup"* ]] || [[ $@ == *"deduplication"* ]];
then
    	source "${SCRIPT_DIR}/dedup.sh"
fi

if [[ $@ == *"downsample"* ]] || [[ $@ == *"downsampling"* ]];
then
    	source "${SCRIPT_DIR}/downsample.sh"
fi

if [[ $@ == *"bigwig"* ]];
then
    	source "${SCRIPT_DIR}/bigwig.sh"
fi

if [[ $@ == *"peak"* ]];
then
    	source "${SCRIPT_DIR}/peaks.sh"
fi

if [[ $@ == *"ngsplot"* ]];
then
    	source "${SCRIPT_DIR}/ngsplot.sh"
fi

if [[ $@ == *"qc"* ]] || [[ $@ == *"QC"* ]];
then
    	source "${SCRIPT_DIR}/qc.sh"
fi
END
