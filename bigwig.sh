#!/bin/bash

echo "Creating BigWig files"
mkdir bigwig	
sorted_bam=$(ls bam/*dedupl-sort-bl.bam 2> /dev/null | head -1 ) #returns empty variable without *dedupl-sort-bl.bam files
if [ -f $sorted_bam ]; 
then
	for file in bam/*dedupl-sort-bl.bam
	do 
		samtools index -@ $max_threads -b $file 2>> bigwig.log
		bigwig_output="${file%-dedupl-sort-bl.bam}-norm.bw"
		bigwig_output=${bigwig_output##*/}
		bamCoverage -p $max_threads --binSize 10 --normalizeUsing RPKM --extendReads 200 --effectiveGenomeSize 2827437033 -b $file -o bigwig/$bigwig_output 2>> bigwig.log
	done
elif [ -z $sorted_bam ]; 
then	
	for file in bam/*sort-bl.bam
	do 
		samtools index -@ $max_threads -b $file 2>> bigwig.log
		bigwig_output="${file%-dedupl-sort-bl.bam}-norm.bw"
		bigwig_output=${bigwig_output##*/}
		bamCoverage -p $max_threads --binSize 10 --normalizeUsing RPKM --extendReads 200 --effectiveGenomeSize 2827437033 -b $file -o bigwig/$bigwig_output 2>> bigwig.log
	done

fi
	
if [[ $# == 1 ]];
	then
		exit 0
