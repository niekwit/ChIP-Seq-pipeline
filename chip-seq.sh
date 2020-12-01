#!/bin/bash

###Author: Niek Wit (University of Cambridge)###

file_path="/home/niek/Documents/analyses/ChIP-Seq/test-data"
cd $file_path
index_path="/home/niek/Documents/references/bowtie2-index/mm9/mm9.genome-index"
blacklist_path="/home/niek/Documents/references/blacklists/Mouse/mm9-blacklist.bed"
PICARD="/home/niek/Documents/scripts/Picard/picard.jar"

fastqc --threads 40 -o fastqc/ raw-data/*fastq.gz 2> chip-seq.log
multiqc -o "fastqc/" "fastqc/" . 2>> chip-seq.log

for i in raw-data/*fastq.gz
do
	file_name=${i##*/} #substring removal of file location
	output_name=${file_name%.fastq.gz} #substring removal of .fastq
	echo "Trimming and mapping ${file_name}"
	trim_galore -j 4 $i -o trim_galore 2>>chip-seq.log #trims off adapter sequences
	trim_galore_output=${i%.fastq.gz}
	trim_galore_output=${trim_galore_output##*/}
	trim_galore_output="trim_galore/${trim_galore_output}_trimmed.fq.gz"
	bowtie2 -p 40 -x $index_path -U $trim_galore_output 2>>chip-seq.log | samtools view -q 15 -bS - > "bam/${output_name}.bam" #creates BAM file
	echo "Removing blacklisted reads from ${output_name}.bam"
	bedtools intersect -v -a "bam/${output_name}.bam" -b $blacklist_path > "bam/${output_name}-bl.bam" -nonamecheck #removes reads from BAM file that are blacklisted
	echo "Performing sorting and deduplication ${output_name}-bl.bam"
	java -jar $PICARD SortSam INPUT="bam/${output_name}-bl.bam" OUTPUT="bam/${output_name}.sorted-bl.bam" SORT_ORDER=coordinate 2>>chip-seq.log #sorts BAM file
	java -jar $PICARD MarkDuplicates INPUT="bam/${output_name}.sorted-bl.bam" OUTPUT="bam/${output_name}.dedup-bl.bam" METRICS_FILE=metric.txt 2>>chip-seq.log #removes duplicates from BAM file
	samtools index -@ 30 -b "bam/${output_name}.dedup-bl.bam" #creates index for deduplicated BAM file
	
done

