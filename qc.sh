#!/bin/bash

echo "QC analysis of ChIP-Seq data"
mkdir chip-qc
sorted_bam=$(ls bam/*dedupl-sort-bl.bam | head -1)
if [ -f $sorted_bam ]; 
then
	for file in bam/*dedupl-sort-bl.bam
	do
		samtools index -@ $max_threads -b $file
	done
	sorted_sample_list=$(ls bam/*dedupl-sort-bl.bam | tr "\n" " ")
	multiBamSummary bins --numberOfProcessors max -b $sorted_sample_list -o chip-qc/multibamsummary.npz 2>> chip-qc.log #deeptools
	plotPCA -in chip-qc/multibamsummary.npz -o chip-qc/PCA_readCounts_all.png -T "PCA of BAM files" 2>> chip-qc.log #deeptools
else
	for file in bam/*-bl.bam
	do
		sort_file="${file%-bl.bam}_sort-bl.bam"
		samtools sort -@ $max_threads $file -o $sort_file
		samtools index -@ $max_threads -b $file
	done
	sorted_sample_list=$(ls bam/*_sort-bl.bam | tr "\n" " ")
	multiBamSummary bins --numberOfProcessors max -b $sorted_sample_list -o chip-qc/multibamsummary.npz 2>> chip-qc.log #deeptools
	plotPCA -in chip-qc/multibamsummary.npz -o chip-qc/PCA_readCounts_all.png -T "Principle Component Analysis" 2>> chip-qc.log #deeptools
fi
plotCorrelation -in chip-qc/multibamsummary.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Average Scores Per Read" --whatToPlot scatterplot -o chip-qc/scatterplot_PearsonCorr_bam-Scores.png --outFileCorMatrix PearsonCorr_bam-Scores.tab #deeptools
plotCorrelation -in chip-qc/multibamsummary.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --colorMap viridis --whatToPlot heatmap -o chip-qc/heatmap_SpearmanCorr_readCounts.png --outFileCorMatrix SpearmanCorr_readCounts.tab #deeptools
if [[ $# == 1 ]];
	then
		exit 0
fi
