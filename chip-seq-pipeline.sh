#!/bin/bash

#This pipeline performs ChIP-Seq analsyis for paired-end data (Niek Wit, University of Cambridge, 2021)

SCRIPT_DIR="/home/niek/Documents/scripts/ChIP-Seq-pipeline/"
PICARD="/home/niek/Documents/scripts/Picard/picard.jar"
#genome=""

usage() {                                    
	echo "Usage: $0 [fastqc] [ align <species> ] [ dedup ] [ pca ] [ bigwig ]"
	echo "fastqc: performs FastQC and MultiQC on fq.gz files"
	echo "align: aligns fq.gz files to selected genome (no deduplication) using HISAT2"
	echo -e "Available genomes for alignment:\n\t<human>: hg38\n\t<mouse>: mm9"
	echo "dedup: removes duplicates using PICARD"
	echo "pca: performs PCA analysis on alignment files"
	echo "bigwig: creates BigWig files"
	exit 2
}

while getopts 'g:?h' c
do
	case $c in
		h|?) usage;;
	esac
done

#determines CPU thread count
max_threads=$(nproc --all)

if [[ $* == *"fastqc"* ]];
then
    echo "Performing FastQC/MultiQC"
	mkdir fastqc    	
	fastqc --threads $max_threads -o fastqc/ raw-data/*fastq.gz 2> chip-seq.log
	multiqc -o "fastqc/" "fastqc/" . 2>> chip-seq.log
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
fi

if [[ $* == *"align"* ]];
then
    	if [[ $* == *"human"* ]];
	then 
		index_path="/home/niek/Documents/references/bowtie2-index/GRCh38.p13.genome-index/GRCh38.p13.genome-index" #hg38
		blacklist_path="/home/niek/Documents/references/blacklists/Human/hg38-blacklist.v2.bed" #hg38
		echo "Human genome (hg38) and blacklist selected"
	elif [[ $* == *"mouse"* ]];
	then	
		index_path="/home/niek/Documents/references/bowtie2-index/mm9/mm9.genome-index" #mm9
		blacklist_path="/home/niek/Documents/references/blacklists/Mouse/mm9-blacklist.bed" #mm9
		echo "Mouse genome (mm9) and blacklist selected"
	fi
	echo "Performing alignment"
	mkdir -p {bam,trim_galore}
	for read1_fastq_gz in raw-data/*_1.fastq.gz
	do
		read2_fastq_gz="${read1_fastq_gz%_1.fastq.gz}_2.fastq.gz"
		echo "Trimming adapter sequences $read1_fastq_gz"
		trim_galore -j 4 -o ./trim_galore --paired $read1_fastq_gz $read2_fastq_gz 2>> chip-seq.log
		read1_val_1_fq_gz="${read1_fastq_gz%_1.fastq.gz}_1_val_1.fq.gz"
		read1_val_1_fq_gz="trim_galore/${read1_val_1_fq_gz##*/}"
		read2_val_2_fq_gz="${read2_fastq_gz%_2.fastq.gz}_2_val_2.fq.gz"
		read2_val_2_fq_gz="trim_galore/${read2_val_2_fq_gz##*/}"
		bowtie2_output=${read1_fastq_gz%_1.fastq.gz}
		bowtie2_output=${bowtie2_output##*/}
		echo "Aligning reads to reference $read1_fastq_gz"
		bowtie2 -p 40 -x $index_path -1 $read1_val_1_fq_gz -2 $read2_val_2_fq_gz 2>> chip-seq.log | samtools view -q 15 -bS -@ $max_threads - > "bam/${bowtie2_output}.bam" 2>> chip-seq.log
		echo "Removing blacklisted regions from alignment $read1_fastq_gz"
		bedtools intersect -v -a "bam/${bowtie2_output}.bam" -b $blacklist_path > "bam/${bowtie2_output}-bl.bam" -nonamecheck 2>> chip-seq.log
	done	
	#count mapped reads:	
	echo "Mapped read counts before deduplication:" >> mapped_read_count_no_dedup.txt	
	for file in bam/*bl.bam
	do
		count=$(samtools view -@ $max_threads -c -F 4 $file)  
		echo "Mapped read count $file: $count" >> mapped_read_count_no_dedup.txt
	done
	echo "Alignment completed. Check if deduplication and/or downsampling is required."	
	if [[ $# == 2 ]];
		then
			exit 0
	fi
fi

if [[ $@ == *"dedup"* ]];
then
    	echo "Performing deduplication"
	for file in bam/*-bl.bam
	do
		echo "Removing duplicates $file"
		sort_output=${file%-bl.bam}
		sort_output="$sort_output-sort-bl.bam"
		java -jar $PICARD SortSam INPUT=$file OUTPUT=$sort_output SORT_ORDER=coordinate 2>> chip-seq.log
		dedup_output=${file%-bl.bam}
		dedup_output="$dedup_output-dedupl-sort-bl.bam"
		java -jar $PICARD MarkDuplicates INPUT=$sort_output OUTPUT=$dedup_output REMOVE_DUPLICATES=TRUE METRICS_FILE=metric.txt 2>> chip-seq.log
	done	
	#count mapped reads after deduplication:	
	echo "Mapped read counts after deduplication:" >> mapped_read_count_dedup.txt	
	for file in bam/*dedupl-sort-bl.bam
	do
		count=samtools view -c -F 4 $file  
		echo "Mapped read count $file: $count" >> mapped_read_count_dedup.txt
	done	
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
fi

if [[ $@ == *"pca"* ]];
then
    	echo "Performing PCA analysis"
	mkdir PCA
	for file in bam/*dedupl-sort-bl.bam
	do 
		sort_file="${file%.-dedupl-sort-bl.bam}_sorted.bam"
		sort_file=${sort_file##*/}
		samtools sort -@ $max_threads $file -o "PCA/$sort_file"
		samtools index -@ $max_threads -b "PCA/$sort_file"
	done
	sorted_sample_list=$(ls PCA/*sorted.bam | tr "\n" " ")
	multiBamSummary bins --numberOfProcessors max -b $sorted_sample_list -o multibamsummary.npz #deeptools
	plotPCA -in multibamsummary_all.npz -o PCA_readCounts_all.png -T "PCA of read counts" #deeptools	
	rm PCA/*.bam
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
fi

if [[ $@ == *"bigwig"* ]];
then
    	echo "Creating BigWig files"
	mkdir bigwig	
	for file in bam/*dedupl-sort-bl.bam
	do 
		samtools index -@ $max_threads -b $file 2>> chip-seq.log
		bigwig_output="${file%-dedupl-sort-bl.bam}-norm.bw"
		bigwig_output=${bigwig_output##*/}
		bamCoverage -p $max_threads --normalizeUsing RPKM -b $file -o bigwig/$bigwig_output 2>> chip-seq.log
	done
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
fi
