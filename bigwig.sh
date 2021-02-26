#!/bin/bash

#load bamCoverage settings:
binsize=$(cat settings.yaml | shyaml get-value BigWig.binSize)
normalizeusing=$(cat settings.yaml | shyaml get-value BigWig.normalizeUsing)
extendreads=$(cat settings.yaml | shyaml get-value BigWig.extendReads)
effectivegenomesize=$(cat settings.yaml | shyaml get-value BigWig.effectiveGenomeSize)

#checking bamCoverage settings:
if [[ ! $binsize =~ ^-?[0-9]+$ ]];
	then
		echo "ERROR: binSize should be an integer."
		exit 1
fi

if [[ $normalizeusing != "RPKM" ]] && [[ $normalizeusing != "CPM" ]] && [[ $normalizeusing != "BPM" ]] && [[ $normalizeusing != "RPGC" ]] && [[ $normalizeusing != "None" ]];
	then
		echo "ERROR: invalid normalisation method chosen."
		echo "Available methods: RPKM, CPM, BPM, RPGC or None"
		exit 1
fi

if [[ ! $extendreads =~ ^-?[0-9]+$ ]];
	then
		echo "ERROR: extendReads should be an integer."
		exit 1
fi

if [[ ! $effectivegenomesize =~ ^-?[0-9]+$ ]];
	then
		echo "ERROR: effectiveGenomeSize should be an integer."
		echo "See https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html"
		exit 1
fi

#creates bigWig files
bigwig_dir=bigwig_bin$binsize
mkdir -p "$bigwig_dir"

function big_wig {
	bigwig_output="${file%-dedupl-sort-bl.bam}-norm.bw"
	bigwig_output=${bigwig_output##*/}
	bamCoverage -p $max_threads --binSize "$binsize" --normalizeUsing "$normalizeusing" --extendReads "$extendreads" --effectiveGenomeSize "$effectivegenomesize" -b $file -o $bigwig_dir/$bigwig_output 2>> bigwig.log
}  

echo "Creating BigWig files"
sorted_bam=$(ls -l bam/*dedupl-sort-bl.bam 2> /dev/null | wc -l) #returns zero without *dedupl-sort-bl.bam files
if [[ $sorted_bam != 0 ]]; 
then
	index_bam=$(ls -l bam/*dedupl-sort-bl.bam.bai 2> /dev/null | wc -l)
	if [[ $index_bam == 0 ]];
	then
		for file in bam/*dedupl-sort-bl.bam
		do 
			samtools index -@ $max_threads -b $file 2>> bigwig.log
			big_wig
		done
	elif [[ $index_bam != 0 ]];
	then
		for file in bam/*dedupl-sort-bl.bam
		do 
			bigwig_output="${file%-dedupl-sort-bl.bam}-norm.bw"
			big_wig
		done
	fi
elif [[ $sorted_bam == 0 ]]; 
then	
	index_bam=$(ls -l bam/*sort-bl.bam.bai 2> /dev/null | wc -l)
	if [[ $index_bam == 0 ]];
	then
		for file in bam/*sort-bl.bam
		do 
			samtools index -@ $max_threads -b $file 2>> bigwig.log
			big_wig
		done
	elif [[ $index_bam != 0 ]];
	then
		for file in bam/*sort-bl.bam
		do 
			big_wig
		done
	fi		
fi
	
if [[ $# == 1 ]];
	then
		exit 0
fi
