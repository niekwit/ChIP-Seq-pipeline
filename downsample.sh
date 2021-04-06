#!/bin/bash

echo "Performing downsampling"
mkdir -p "downsample"
dedup_bam=$(ls bam/*dedupl-sort-bl.bam 2> /dev/null | head -1 ) #returns empty variable without *dedupl-sort-bl.bam files and no error message
if [[ -f $dedup_bam ]]; 
then
	read_count_file="mapped_read_count_dedup.txt"
	indexed_bam=$(ls bam/*.bam.bai 2> /dev/null | head -1 ) #returns empty variable if bam files haven't been indexed and no error message
	if [[ -z "$indexed_bam" ]]; 
	then
		echo "No indexed bam files found, indexing..."
		for file in bam/*dedupl-sort-bl.bam
		do
			samtools index -@ "$max_threads" -b "$file"
		done
	elif [[ -f $indexed_bam ]];
	then
		echo "Bam files already indexed"
	fi
elif [[ -z $dedup_bam ]]; 
then
	read_count_file="mapped_read_count_no_dedup.txt"
	indexed_bam=$(ls bam/*.bam.bai 2> /dev/null | head -1 ) #returns empty variable if bam files haven't been indexed and no error message
	if [[ -z "$indexed_bam" ]]; 
	then
		echo "No indexed bam files found, indexing..."
		for file in bam/*-sort-bl.bam
		do
			samtools index -@ "$max_threads" -b "$file"
		done	
	elif [[ -f $indexed_bam ]];
	then
		echo "Bam files already indexed"
	fi
fi

working_dir=$(pwd)
Rscript "${SCRIPT_DIR}sample-size.R" "$working_dir" "$read_count_file" #creates file that contains scaling factors for each input file and file:chip pair
sed '1d' downsample/scaling_factors.txt > downsample/scaling_factors-temp.txt #removes header
input="downsample/scaling_factors-temp.txt"
count=0
totalfiles=$(wc -l < "downsample/scaling_factors-temp.txt")
while IFS= read -r line
do
	count=$(($count + 1))
	echo -ne "\rDownscaling file $count/$totalfiles"
	scaling_factor=$(echo $line | cut -d "," -f 1)
	scaling_input=$(echo $line | cut -d "," -f 2)
	scaling_output_name=$(echo $line | cut -d "," -f 3)
	samtools view -@ "$max_threads" -s "$scaling_factor" -b "$scaling_input" -o "$scaling_output_name"
done < "$input"
echo -e "\n"

rm "downsample/scaling_factors-temp.txt"

if [[ $# == 1 ]];
	then
		exit 0
fi
