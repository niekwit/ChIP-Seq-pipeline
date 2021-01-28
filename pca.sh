#!/bin/bash

pca(){
echo "Performing PCA analysis"
mkdir pca
sorted_bam=$(ls bam/*dedupl-sort-bl.bam | head -1)
if [ -f $sorted_bam ]; 
then
	for file in bam/*dedupl-sort-bl.bam
	do
		samtools index -@ $max_threads -b $file
	done
	sorted_sample_list=$(ls bam/*dedupl-sort-bl.bam | tr "\n" " ")
	multiBamSummary bins --numberOfProcessors max -b $sorted_sample_list -o pca/multibamsummary.npz 2>> chip-seq.log #deeptools
	plotPCA -in pca/multibamsummary.npz -o pca/PCA_readCounts_all.png -T "PCA of BAM files" 2>> chip-seq.log #deeptools
else
	for file in bam/*-bl.bam
	do
		sort_file="${file%-bl.bam}_sort-bl.bam"
		samtools sort -@ $max_threads $file -o $sort_file
		samtools index -@ $max_threads -b $file
	done
	sorted_sample_list=$(ls bam/*_sort-bl.bam | tr "\n" " ")
	multiBamSummary bins --numberOfProcessors max -b $sorted_sample_list -o pca/multibamsummary.npz 2>> chip-seq.log #deeptools
	plotPCA -in pca/multibamsummary.npz -o pca/PCA_readCounts_all.png -T "Principle Component Analysis" 2>> chip-seq.log #deeptools
fi
if [[ $# == 1 ]];
	then
		exit 0
fi
}
