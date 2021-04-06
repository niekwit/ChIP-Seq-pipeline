#!/bin/bash
align_setting=""

threads=$1
align_prog=$2
genome=$3

#mkdir -p {bam,trim_galore}

#selects relevant reference files
if [[ $align_prog == *"hisat2"* ]];then
	align_setting="hisat2"
elif [[ $align_prog == *"bwa"* ]];then
	align_setting="bwa"
fi

index_path=$(cat "$SCRIPT_DIR/settings.yaml" | shyaml get-value $align_setting.$genome)
blacklist_path=$(cat "$SCRIPT_DIR/settings.yaml" | shyaml get-value blacklist.$genome)

#performs trimming, alignment, removal of blacklisted regions, and sorts
#to do: perform trimming seperately; fix se-alignments
if [[ $align_prog == *"-se"* ]];then
	if [[ ! -d ./trim_galore ]];then
		echo "Trimming fastq files"
		for file in raw-data/*fastq.gz 
		do
			trim_galore -j 4 -o ./trim_galore $file  #output extension will be "_trimmed.fq.gz"
		done
	else
		echo "Skipping trimming (already performed)"
	fi
elif [[ $align_prog == *"-pe"* ]];then #not proper yet
	if [[ ! -d ./trim_galore ]];then
		for file in raw-data/*fastq.gz 
		do
			echo "Trimming fastq files"
			trim_galore -j 4 -o ./trim_galore $file  #output extension will be "_trimmed.fq.gz"
		done
	else
		echo "Skipping trimming (already performed)"
	fi
fi

if [[ $align_prog == *"hisat2-pe"* ]]; then
	echo "Paired-end alignment selected with HISAT2"
	echo "Aligning fastq files"
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
		hisat2 -p "$threads" -x $index_path -1 "$read1_val_1_fq_gz" -2 "$read2_val_2_fq_gz" 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$threads" - | bedtools intersect -v -a "stdin" -b "$blacklist_path" -nonamecheck | samtools sort -@ "$threads" - > "$hisat2_output" #bam file will not contain unmapped and multimapped reads 
	done
elif [[ $align_prog == *"hisat2-se"* ]]; then
	echo "Single-end alignment selected with HISAT2"
	echo "Aligning fastq files"
	mkdir -p bam/
	for file in trim_galore/*trimmed.fq.gz
	do
		echo $file
		hisat2_output="${file%_trimmed.fq.gz}-sort-bl.bam"
		hisat2_output="bam/${hisat2_output##*/}"
		zcat $file | hisat2 -p "$threads" -x "$index_path" - 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$threads" - | bedtools intersect -v -a "stdin" -b "$blacklist_path" -nonamecheck | samtools sort -@ "$threads" - > "$hisat2_output" #bam file will not contain unmapped and multimapped reads
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

: <<'END'
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
elif [[ $* == *"hisat-se"* ]]; #not working
then
	echo "Single-end alignment selected with HISAT2"
	echo "Aligning fastq files"
	for file in trim_galore/*trimmed.fq.gz
	do
		hisat2_output="${file%_trimmed.fq.gz}-sort-bl.bam"
		hisat2_output="bam/${hisat2_output##*/}"
		hisat2 -p "$max_threads" -x "$index_path" $file 2>> align.log | samtools view -q 15 -F 260 -bS -@ "$max_threads" - | bedtools intersect -v -a "stdin" -b "$blacklist_path" -nonamecheck | samtools sort -@ "$max_threads" - > "$hisat2_output" #bam file will not contain unmapped and multimapped reads
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
elif [[ $* == *"bwa-se"* ]]; #not working
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
END
