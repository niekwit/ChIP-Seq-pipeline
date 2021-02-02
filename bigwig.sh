#!/bin/bash

echo "Creating BigWig files"
mkdir bigwig	
for file in bam/*dedupl-sort-bl.bam
do 
	samtools index -@ $max_threads -b $file 2>> bigwig.log
	bigwig_output="${file%-dedupl-sort-bl.bam}-norm.bw"
	bigwig_output=${bigwig_output##*/}
	bamCoverage -p $max_threads --binSize 10 --normalizeUsing RPKM --extendReads 200 --effectiveGenomeSize 2827437033 -b $file -o bigwig/$bigwig_output 2>> bigwig.log
done
if [[ $# == 1 ]];
	then
		exit 0
fi
