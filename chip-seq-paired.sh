#!/bin/bash

###Author: Niek Wit (University of Cambridge)###

################################################################################################################
#Instructions:													   #
#- Create a main folder for the ChIP-Seq experiment with all .fastq.gz files in a sub-folder called raw-data   #
#- From the command line run:$ ./chip-seq-paired.sh <path/to/main/folder> <species>				   #
#	-Species options: mouse										   #
#			  human										   #
################################################################################################################

if (($# != 2));
	then 
		>&2 echo "Error: please provide two arguments as follows: $ ./chip-seq-paired.sh <path/to/main/folder> <species>"
		exit
fi

file_path=$1
cd $file_path
mkdir -p {bam,fastqc,trim_galore,bigwig}

if [ $2 = "human" ];
	then 
		index_path="/home/niek/Documents/references/bowtie2-index/GRCh38.p13.genome-index/GRCh38.p13.genome-index" #hg38
		blacklist_path="/home/niek/Documents/references/blacklists/Human/hg38-blacklist.v2.bed" #hg38
		echo "Human genome (hg38) and blacklist selected"
elif [ $2 = "mouse" ];
	then	
		index_path="/home/niek/Documents/references/bowtie2-index/mm9/mm9.genome-index" #mm9
		blacklist_path="/home/niek/Documents/references/blacklists/Mouse/mm9-blacklist.bed" #mm9
		echo "Mouse genome (mm9) and blacklist selected"
fi

PICARD="/home/niek/Documents/scripts/Picard/picard.jar"

echo "Performing QC"
fastqc --threads 40 -o fastqc/ raw-data/*fastq.gz 2> chip-seq.log
multiqc -o "fastqc/" "fastqc/" . 2>> chip-seq.log

for read1_fastq_gz in raw-data/*_1.fastq.gz
do
	
	read2_fastq_gz="${read1_fastq_gz%_1.fastq.gz}_2.fastq.gz"
	echo "Trimming $read1_fastq_gz"
	trim_galore -j 4 -o ./trim_galore --paired $read1_fastq_gz $read2_fastq_gz 2>> chip-seq.log
	read1_val_1_fq_gz="${read1_fastq_gz%_1.fastq.gz}_1_val_1.fq.gz"
	read1_val_1_fq_gz="trim_galore/${read1_val_1_fq_gz##*/}"
	read2_val_2_fq_gz="${read2_fastq_gz%_2.fastq.gz}_2_val_2.fq.gz"
	read2_val_2_fq_gz="trim_galore/${read2_val_2_fq_gz##*/}"
	bowtie2_output=${read1_fastq_gz%_1.fastq.gz}
	bowtie2_output=${bowtie2_output##*/}
	echo "Aligning reads to reference $read1_fastq_gz"
	bowtie2 -p 40 -x $index_path -1 $read1_val_1_fq_gz -2 $read2_val_2_fq_gz 2>> chip-seq.log | samtools view -q 15 -bS - > "bam/${bowtie2_output}.bam" 2>> chip-seq.log
	echo "Removing blacklisted regions from alignment $read1_fastq_gz"
	bedtools intersect -v -a "bam/${bowtie2_output}.bam" -b $blacklist_path > "bam/${bowtie2_output}-bl.bam" -nonamecheck 2>> chip-seq.log #removes reads from BAM file that are blacklisted
	echo "Removing duplicates $read1_fastq_gz"
	java -jar $PICARD SortSam INPUT="bam/${bowtie2_output}-bl.bam" OUTPUT="bam/${bowtie2_output}-sort-bl.bam" SORT_ORDER=coordinate 2>> chip-seq.log
	java -jar $PICARD MarkDuplicates INPUT="bam/${bowtie2_output}-sort-bl.bam" OUTPUT="bam/${bowtie2_output}-dedupl-sort-bl.bam" REMOVE_DUPLICATES=TRUE METRICS_FILE=metric.txt 2>> chip-seq.log
	samtools_input="bam/${bowtie2_output}-dedupl-sort-bl.bam"
	echo "Creating bigWig file $read1_fastq_gz"
	samtools index -@ 40 -b $samtools_input 2>> chip-seq.log
	bam_input="bam/${bowtie2_output}-dedupl-sort-bl.bam"
	bigwig_output="${samtools_input%-dedupl-sort-bl.bam}-norm.bw"
	bigwig_output=${bigwig_output##*/}
	bamCoverage -p 40 --normalizeUsing RPKM -b $bam_input -o bigwig/$bigwig_output 2>> chip-seq.log
	echo "Finished analysing $read1_fastq_gz"
	
done
