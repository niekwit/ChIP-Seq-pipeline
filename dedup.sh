#!/bin/bash

PICARD=$1

echo "Performing deduplication"
for file in bam/*-sort-bl.bam
do
	echo "Removing duplicates $file"
	dedup_output=${file%-sort-bl.bam}
	dedup_output="$dedup_output-dedupl-sort-bl.bam"
	java -jar $PICARD MarkDuplicates INPUT=$file OUTPUT=$dedup_output REMOVE_DUPLICATES=TRUE METRICS_FILE=metric.txt 2>> deduplication.log
done
echo "Uniquely mapped read counts after deduplication:" >> mapped_read_count_dedup.txt	
for file in bam/*dedupl-sort-bl.bam
do
	count=$(samtools view -c -F 260 $file) #bit 3 and 9 SAM FLAGs
	echo "$file: $count" >> mapped_read_count_dedup.txt
done
