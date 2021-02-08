#!/bin/bash


mkdir -p {bam,trim_galore}

#selects relevant reference files
if [[ "$*" == *"hisat"* ]] && [[ "$*" == *"human"* ]];
then
	index_path="$index_path_hs1" #hg19
	blacklist_path="$blacklist_path_hg19" #hg19
	touch homer_hg19 #file for loading correct genome with downstream applications
	echo "Human genome (hg19) and blacklist selected"
elif [[ "$*" == *"hisat"* ]] && [[ "$*" == *"mouse"* ]];
then	
	index_path="$index_path_mm1" #mm9
	blacklist_path="$blacklist_path_mm9" #mm9
	touch homer_mm9 #file for loading correct genome with downstream applications
	echo "Mouse genome (mm9) and blacklist selected"
elif [[ "$*" == *"bwa"* ]] && [[ "$*" == *"human"* ]]; ###work in progress###
then
	index_path="$index_path_hs2" #hg19
	blacklist_path="$blacklist_path_hg19" #hg19
	touch homer_hg19 #file for loading correct genome with downstream applications
	echo "Human genome (hg19) and blacklist selected"
elif [[ "$*" == *"bwa"* ]] && [[ "$*" == *"mouse"* ]];
then
	index_path="$index_path_mm2" #mm9
	blacklist_path="$blacklist_path_mm9" #mm9
	touch homer_mm9 #file for loading correct genome with downstream applications
	echo "Mouse genome (mm9) and blacklist selected"
fi

#performs trimming, alignment, removal of blacklisted regions, and sorts	
if [[ $* == *"hisat-pe"* ]];
then
	echo "Paired-end alignment selected with HISAT2"
	echo "Trimming and aligning fastq files"
	for read1_fastq_gz in raw-data/*_1.fastq.gz
	do
		read2_fastq_gz="${read1_fastq_gz%_1.fastq.gz}_2.fastq.gz"
		trim_galore -j 4 -o ./trim_galore --paired $read1_fastq_gz $read2_fastq_gz 2>> align.log
		read1_val_1_fq_gz="${read1_fastq_gz%_1.fastq.gz}_1_val_1.fq.gz"
		read1_val_1_fq_gz="trim_galore/${read1_val_1_fq_gz##*/}"
		read2_val_2_fq_gz="${read2_fastq_gz%_2.fastq.gz}_2_val_2.fq.gz"
		read2_val_2_fq_gz="trim_galore/${read2_val_2_fq_gz##*/}"
		hisat2_output="${read1_fastq_gz%_1.fastq.gz}-sort-bl.bam"
		hisat2_output="bam/${read1_fastq_gz##*/}"
		hisat2 -p "$max_threads" -x $index_path -1 "$read1_val_1_fq_gz" -2 "$read2_val_2_fq_gz" 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$max_threads" - | bedtools intersect -v -a "stdin" -b "$blacklist_path" -nonamecheck | samtools sort -@ "$max_threads" - > "$hisat2_output" #bam file will not contain unmapped and multimapped reads 
	done
elif [[ $* == *"hisat-se"* ]];
then
	echo "Single-end alignment selected with HISAT2"
	echo "Trimming and aligning fastq files"
	for file in raw-data/*.fastq.gz
	do
		trim_galore -j 4 -o ./trim_galore "$file" 2>> align.log #output extension will be "_trimmed.fq.gz"
		hisat2_input="${file%.fastq.gz}_trimmed.fq.gz"
		hisat2_input="trim_galore/${hisat2_input##*/}"
		hisat2_output="bam/${file%_1.fastq.gz}-sort-bl.bam"
		hisat2 -p "$max_threads" -x "$index_path" "$hisat2_input" 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$max_threads" - | bedtools intersect -v -a "stdin" -b "$blacklist_path" -nonamecheck | samtools sort -@ "$max_threads" - > "$hisat2_output" #bam file will not contain unmapped and multimapped reads
	done
elif [[ $* == *"bwa-pe"* ]];
then
	echo "Paired-end alignment selected with bwa"
	echo "Trimming and aligning fastq files"
	for read1_fastq_gz in raw-data/*_1.fastq.gz
	do
		read2_fastq_gz="${read1_fastq_gz%_1.fastq.gz}_2.fastq.gz"
		trim_galore -j 4 -o ./trim_galore --paired $read1_fastq_gz $read2_fastq_gz 2>> align.log
		read1_val_1_fq_gz="${read1_fastq_gz%_1.fastq.gz}_1_val_1.fq.gz"
		read1_val_1_fq_gz="trim_galore/${read1_val_1_fq_gz##*/}"
		read2_val_2_fq_gz="${read2_fastq_gz%_2.fastq.gz}_2_val_2.fq.gz"
		read2_val_2_fq_gz="trim_galore/${read2_val_2_fq_gz##*/}"
		bwa_output="${read1_fastq_gz%_1.fastq.gz}-sort-bl.bam"
		bwa_output="bam/${bwa_output##*/}"
		bwa mem "$index_path" -t "$max_threads" "$read1_val_1_fq_gz" "$read2_val_2_fq_gz" | samtools view -q 15 -F 260 -bS -@ "$max_threads" - | bedtools intersect -v -a "stdin" -b "$blacklist_path" -nonamecheck | samtools sort -@ "$max_threads" - > "$bwa_output" #bam file will not contain unmapped and multimapped reads
	done
elif [[ $* == *"bwa-se"* ]];
then
	echo "Single-end alignment selected with bwa"
	echo "Trimming and aligning fastq files"
	for file in raw-data/*.fastq.gz
	do
		trim_galore -j 4 -o ./trim_galore $file 2>> align.log #output extension will be "_trimmed.fq.gz"
		bwa_input="${file%.fastq.gz}_trimmed.fq.gz"
		bwa_input="trim_galore/${hisat2_input##*/}"
		bwa_output="bam/${file%_1.fastq.gz}-sort-bl.bam"
		bwa mem $index_path -t $max_threads $bwa_input 2>> align.log | samtools view -q 15 -F 260 -bS -@ $max_threads - | bedtools intersect -v -a "stdin" -b $blacklist_path -nonamecheck | samtools sort -@ $max_threads - > $bwa_output #bam file will not contain unmapped and multimapped reads
	done
fi	

#count mapped reads:	
echo "Uniquely mapped read counts before deduplication:" >> mapped_read_count_no_dedup.txt	
for file in bam/*-sort-bl.bam
do
	count=$(samtools view -@ $max_threads -c -F 260 $file) #bit 3 and 9 SAM FLAGs
	echo "$file: $count" >> mapped_read_count_no_dedup.txt
done

echo "Alignment completed. Check if deduplication and/or downsampling is required."	
#add code thats can suggest if samples need downsampling?

if [[ $# == 3 ]];
	then
		exit 0
fi
