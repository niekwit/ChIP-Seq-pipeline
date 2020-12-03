#!/bin/bash

###Author: Niek Wit (University of Cambridge)###

file_path="/home/niek/Documents/analyses/ChIP-Seq/test-data"
cd $file_path
index_path="/home/niek/Documents/references/bowtie2-index/GRCh38.p13.genome-index/GRCh38.p13.genome-index" #hg38
#blacklist_path="/home/niek/Documents/references/blacklists/Mouse/mm9-blacklist.bed" #mm9 
blacklist_path="/home/niek/Documents/references/blacklists/Human/hg38-blacklist.v2.bed" #hg38
PICARD="/home/niek/Documents/scripts/Picard/picard.jar"

echo "Performing QC"
fastqc --threads 40 -o fastqc/ raw-data/*fastq.gz 2> chip-seq.log
multiqc -o "fastqc/" "fastqc/" . 2>> chip-seq.log

for read1_fastq_gz in raw-data/*_1.fastq.gz
do
	
	read2_fastq_gz="${read1_fastq_gz%_1.fastq.gz}_2.fastq.gz"
	trim_galore -j 4 -o ./trim_galore --paired $read1_fastq_gz $read2_fastq_gz
	read1_val_1_fq_gz="${read1_fastq_gz%_1.fastq.gz}_1_val_1.fq.gz"
	read1_val_1_fq_gz="trim_galore/${read1_val_1_fq_gz##*/}"
	read2_val_2_fq_gz="${read2_fastq_gz%_2.fastq.gz}_2_val_2.fq.gz"
	read2_val_2_fq_gz="trim_galore/${read2_val_2_fq_gz##*/}"
	bowtie2_output=${read1_fastq_gz%_1.fastq.gz}
	bowtie2_output=${bowtie2_output##*/}
	bowtie2 -p 40 -x $index_path -1 $read1_val_1_fq_gz -2 $read2_val_2_fq_gz 2>> chip-seq.log | samtools view -q 15 -bS - > "bam/${bowtie2_output}.bam" 2>> chip-seq.log
	bedtools intersect -v -a "bam/${bowtie2_output}.bam" -b $blacklist_path > "bam/${bowtie2_output}-bl.bam" -nonamecheck 2>> chip-seq.log #removes reads from BAM file that are blacklisted
	java -jar $PICARD SortSam INPUT="bam/${bowtie2_output}-bl.bam" OUTPUT="bam/${bowtie2_output}-sort-bl.bam" SORT_ORDER=coordinate 2>> chip-seq.log
	java -jar $PICARD MarkDuplicates INPUT="bam/${bowtie2_output}-sort-bl.bam" OUTPUT="bam/${bowtie2_output}-dedupl-sort-bl.bam" REMOVE_DUPLICATES=TRUE METRICS_FILE=metric.txt 2>> chip-seq.log
	samtools_input="bam/${bowtie2_output}-dedupl-sort-bl.bam"
	samtools index -@ 40 -b $samtools_input 2>> chip-seq.log
	
	bam_input="bam/${bowtie2_output}-dedupl-sort-bl.bam"
	bigwig_output="${samtools_input%-dedupl-sort-bl.bam}-norm.bw"
	bigwig_output=${bigwig_output##*/}
	bamCoverage -p 40 --normalizeUsing RPKM -b $bam_input -o bigwig/$bigwig_output 2>> chip-seq.log
	
done



