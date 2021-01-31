#!/bin/bash

align(){
if [[ $* == *"hisat"* ]];
	then	
	if [[ $* == *"human"* ]];
	then 
		index_path="/home/niek/Documents/references/hisat2-index/GRCh37-hg19ucsc/hg19-index" #hg19
		blacklist_path="/home/niek/Documents/references/blacklists/Human/hg19/wgEncodeDukeMapabilityRegionsExcludable.bed" #hg19
		touch homer_hg19 #file for loading correct genome with downstream applications
		echo "Human genome (hg19) and blacklist selected"
	elif [[ $* == *"mouse"* ]];
	then	
		index_path="/home/niek/Documents/references/bowtie2-index/mm9/mm9.genome-index" #mm9
		blacklist_path="/home/niek/Documents/references/blacklists/Mouse/mm9-blacklist.bed" #mm9
		touch homer_mm9 #file for loading correct genome with downstream applications
		echo "Mouse genome (mm9) and blacklist selected"
	fi
	
	mkdir -p {bam,trim_galore}
	for read1_fastq_gz in raw-data/*_1.fastq.gz
	do
		read2_fastq_gz="${read1_fastq_gz%_1.fastq.gz}_2.fastq.gz"
		echo "Trimming adapter sequences $read1_fastq_gz"
		trim_galore -j 4 -o ./trim_galore --paired $read1_fastq_gz $read2_fastq_gz 2>> align.log
		read1_val_1_fq_gz="${read1_fastq_gz%_1.fastq.gz}_1_val_1.fq.gz"
		read1_val_1_fq_gz="trim_galore/${read1_val_1_fq_gz##*/}"
		read2_val_2_fq_gz="${read2_fastq_gz%_2.fastq.gz}_2_val_2.fq.gz"
		read2_val_2_fq_gz="trim_galore/${read2_val_2_fq_gz##*/}"
		hisat2_output=${read1_fastq_gz%_1.fastq.gz}
		hisat2_output=${hisat2_output##*/}
		echo "Aligning reads to reference $read1_fastq_gz"
		hisat2 -p $max_threads -x $index_path -1 $read1_val_1_fq_gz -2 $read2_val_2_fq_gz 2>> align.log | samtools view -q 15 -F 260 -bS -@ $max_threads - > "bam/${hisat2_output}.bam" 2>> align.log
		echo "Removing blacklisted regions from alignment $read1_fastq_gz"
		bedtools intersect -v -a "bam/${hisat2_output}.bam" -b $blacklist_path > "bam/${hisat2_output}-bl.bam" -nonamecheck 2>> align.log #combine into new pipe: hisat2 | samtools view | bedtools intersect | samtools sort
	done	
	#count mapped reads:	
	echo "Uniquely mapped read counts before deduplication:" >> mapped_read_count_no_dedup.txt	
	for file in bam/*bl.bam
	do
		count=$(samtools view -@ $max_threads -c -F 260 $file) #bit 3 and 9 SAM FLAGs
		echo "$file: $count" >> mapped_read_count_no_dedup.txt
	done
elif [[ $* == *"bwa"* ]]; ###not working yet###
	then
		then	
	if [[ $* == *"human"* ]];
	then 
		index_path="###" #hg19
		blacklist_path="/home/niek/Documents/references/blacklists/Human/hg19/wgEncodeDukeMapabilityRegionsExcludable.bed" #hg19
		touch homer_hg19 #file for loading correct genome with downstream applications
		echo "Human genome (hg19) and blacklist selected"
	elif [[ $* == *"mouse"* ]];
	then	
		index_path="###" #mm9
		blacklist_path="/home/niek/Documents/references/blacklists/Mouse/mm9-blacklist.bed" #mm9
		touch homer_mm9 #file for loading correct genome with downstream applications
		echo "Mouse genome (mm9) and blacklist selected"
	fi
	
	mkdir -p {bam,trim_galore}
	for read1_fastq_gz in raw-data/*_1.fastq.gz
	do
		read2_fastq_gz="${read1_fastq_gz%_1.fastq.gz}_2.fastq.gz"
		echo "Trimming adapter sequences $read1_fastq_gz"
		trim_galore -j 4 -o ./trim_galore --paired $read1_fastq_gz $read2_fastq_gz 2>> align.log
		read1_val_1_fq_gz="${read1_fastq_gz%_1.fastq.gz}_1_val_1.fq.gz"
		read1_val_1_fq_gz="trim_galore/${read1_val_1_fq_gz##*/}"
		read2_val_2_fq_gz="${read2_fastq_gz%_2.fastq.gz}_2_val_2.fq.gz"
		read2_val_2_fq_gz="trim_galore/${read2_val_2_fq_gz##*/}"
		hisat2_output=${read1_fastq_gz%_1.fastq.gz}
		hisat2_output=${hisat2_output##*/}
		echo "Aligning reads to reference $read1_fastq_gz"
		#bwa mem -t $max_threads 
		#bwa 		
		#hisat2 -p $max_threads -x $index_path -1 $read1_val_1_fq_gz -2 $read2_val_2_fq_gz 2>> align.log | samtools view -q 15 -F 260 -bS -@ $max_threads - > "bam/${hisat2_output}.bam" 2>> align.log
		echo "Removing blacklisted regions from alignment $read1_fastq_gz"
		bedtools intersect -v -a "bam/${hisat2_output}.bam" -b $blacklist_path > "bam/${hisat2_output}-bl.bam" -nonamecheck 2>> align.log
	done	
	#count mapped reads:	
	echo "Uniquely mapped read counts before deduplication:" >> mapped_read_count_no_dedup.txt	
	for file in bam/*bl.bam
	do
		count=$(samtools view -@ $max_threads -c -F 260 $file) #bit 3 and 9 SAM FLAGs
		echo "$file: $count" >> mapped_read_count_no_dedup.txt
	done
fi

echo "Alignment completed. Check if deduplication and/or downsampling is required."	
#add code thats can suggest if samples need downsampling?
if [[ $# == 3 ]];
	then
		exit 0
fi
}
