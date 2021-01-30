#!/bin/bash

dedup() {
echo "Performing deduplication"
	for file in bam/*-bl.bam
	do
		echo "Removing duplicates $file"
		sort_output=${file%-bl.bam}
		sort_output="$sort_output-sort-bl.bam"
		java -jar $PICARD SortSam INPUT=$file OUTPUT=$sort_output SORT_ORDER=coordinate 2>> deduplication.log
		dedup_output=${file%-bl.bam}
		dedup_output="$dedup_output-dedupl-sort-bl.bam"
		java -jar $PICARD MarkDuplicates INPUT=$sort_output OUTPUT=$dedup_output REMOVE_DUPLICATES=TRUE METRICS_FILE=metric.txt 2>> deduplication.log
	done	
	#count mapped reads after deduplication:	
	echo "Uniquely mapped read counts after deduplication:" >> mapped_read_count_dedup.txt	
	for file in bam/*dedupl-sort-bl.bam
	do
		count=$(samtools view -c -F 260 $file) #bit 3 and 9 SAM FLAGs
		echo "$file: $count" >> mapped_read_count_dedup.txt
	done	
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
}
