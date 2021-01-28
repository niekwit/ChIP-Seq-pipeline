#!/bin/bash

dedup(){
echo "Performing PCA analysis"
mkdir pca
sorted_bam="bam/*dedupl-sort-bl.bam"
if [ -f $sorted_bam ]; 
then
	for file in bam/*dedupl-sort-bl.bam
	do
		samtools index -@ $max_threads -b "bam/$file"
		sorted_sample_list=$(ls bam/*dedupl-sort-bl.bam | tr "\n" " ")
		multiBamSummary bins --numberOfProcessors max -b $sorted_sample_list -o pca/multibamsummary.npz 2>> chip-seq.log #deeptools
		plotPCA -in pca/multibamsummary.npz -o pca/PCA_readCounts_all.png -T "PCA of BAM files" 2>> chip-seq.log #deeptools
	done
else
	for file in bam/*-bl.bam
	do
		sort_file="${file%-bl.bam}_sort-bl.bam"
		samtools sort -@ $max_threads $file -o "bam/$sort_file"
		samtools index -@ $max_threads -b "bam/$file"
		sorted_sample_list=$(ls bam/*_sort-bl.bam | tr "\n" " ")
		multiBamSummary bins --numberOfProcessors max -b $sorted_sample_list -o pca/multibamsummary.npz 2>> chip-seq.log #deeptools
		plotPCA -in pca/multibamsummary.npz -o pca/PCA_readCounts_all.png -T "PCA of BAM files" 2>> chip-seq.log #deeptools
	done
if [[ $# == 1 ]];
	then
		exit 0
fi
}
