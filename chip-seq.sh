#!/bin/bash

file_path="/home/niek/Documents/analyses/ChIP-Seq/test-data"
cd $file_path
index_path="/home/niek/Documents/references/bowtie2-index/mm9/mm9.genome-index"

fastqc --threads 40 -o fastqc/ raw-data/*fastq.gz 2> chip-seq.log
multiqc -o "fastqc/" "fastqc/" . 2>> chip-seq.log

for i in raw-data/*fastq.gz
do
	file_name=${i##*/} #substring removal of file location
	output_name=${file_name%.fastq.gz} #substring removal of .fastq
	echo "Trimming and mapping ${file_name}"
	trim_galore -j 4 $i -o trim_galore 2>>chip-seq.log
	trim_galore_output=${i%.fastq.gz}
	trim_galore_output=${trim_galore_output##*/}
	trim_galore_output="trim_galore/${trim_galore_output}_trimmed.fq.gz"
	bowtie2 -p 40 -x $index_path -U $trim_galore_output 2>>chip-seq.log | samtools view -q 15 -bS - > "bam/${output_name}.bam"
	echo "Performing sorting and deduplication ${output_name}.bam"
	java -jar /home/niek/Documents/scripts/Picard/picard.jar SortSam INPUT="bam/${output_name}.bam" OUTPUT="bam/${output_name}.sorted.bam" SORT_ORDER=coordinate 2>>chip-seq.log
	java -jar /home/niek/Documents/scripts/Picard/picard.jar MarkDuplicates INPUT="bam/${output_name}.sorted.bam" OUTPUT="bam/${output_name}.dedup.bam" METRICS_FILE=metric.txt 2>>chip-seq.log
done

