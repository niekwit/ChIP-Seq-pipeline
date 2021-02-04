#!/bin/bash

echo "Generating metagene plots and heatmaps with ngs.plot.r"
mkdir ngsplot 
if [[ -f homer_hg19 ]]; 
	then
		ngsplot_genome=hg19
elif [[ -f homer_mm_9 ]]; 
	then
		ngsplot_genome=mm9
fi

for file in bam/*-dedupl-sort-bl.bam #generate plots only for non-input files
do
	case $file in 
		*input*) continue;;
		*) 
			ngs_output=${file%-dedupl-sort-bl.bam}
			ngs_output=${ngs_output##*/}
			mkdir "ngsplot/$ngs_output"
			ngs.plot.r -G $ngsplot_genome -R tss -C $file -O "ngsplot/$ngs_output/${ngs_output}_tss" -T $ngs_output -L 5000 2>> ngsplot.log
	esac
done

if [[ -f downsample/scaling_factors.txt ]]; 
then
	sed '1d' downsample/scaling_factors.txt > downsample/scaling_factors-temp.txt #removes header
	input="downsample/scaling_factors-temp.txt"
	count=0
	while IFS= read -r line
	do
		input=$(echo $line | cut -d "," -f 2)
		chip=$(echo $line | cut -d "," -f 4)
		scaling_output_name=$(echo $line | cut -d "," -f 3)
		scaling_output_name=${scaling_output_name##*/}
		scaling_output_name=${scaling_output_name%-ds.bam}
		ngs.plot.r -G $ngsplot_genome -R tss -C "$input:$chip" -O "ngsplot/$scaling_output_name/${scaling_output_name}_tss" -T $scaling_output_name -L 5000 2>> ngsplot.log
	done < "$input"


if [[ -f ngsplot.conf ]] && [[ ! -f downsample/scaling_factors.txt ]]; 
then
	sed '1d' ngsplot.conf > ngsplot-temp.conf #removes header
	input="ngsplot-temp.conf"
	while IFS= read -r line
	do
		ngsplot_chip=$(echo $line | cut -d "," -f 1)
		ngsplot_input=$(echo $line | cut -d "," -f 2)
		ngsplot_output_name=$(echo $line | cut -d "," -f 3)
		mkdir "ngsplot/$ngsplot_output_name"
		ngs.plot.r -G $ngsplot_genome -R tss -C "$ngsplot_chip:$ngsplot_input" -O "ngsplot/$ngsplot_output_name/${ngsplot_output_name}_tss" -T $ngsplot_output_name -L 5000 2>> ngsplot.log
	done < "$input"
fi

rm ngsplot-temp.conf *.cnt
if [[ $# == 1 ]];
	then
		exit 0
fi

