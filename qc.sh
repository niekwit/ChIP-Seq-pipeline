#!/bin/bash

threads=$1

echo "QC analysis of ChIP-Seq data"
mkdir chip-qc
sorted_bam=$(ls bam/*dedupl-sort-bl.bam 2> /dev/null | head -1 ) #returns empty variable without *dedupl-sort-bl.bam files and no error message
if [ -f $sorted_bam ]; 
then
	for file in bam/*dedupl-sort-bl.bam
	do
		samtools index -@ $threads -b $file
	done
	sorted_sample_list=$(ls bam/*dedupl-sort-bl.bam | tr "\n" " ")
	multiBamSummary bins --numberOfProcessors $threads -b $sorted_sample_list -o chip-qc/multibamsummary.npz 2>> chip-qc.log #deeptools
	plotPCA -in chip-qc/multibamsummary.npz -o chip-qc/PCA_readCounts_all.png -T "Principle Component Analysis" 2>> chip-qc.log #deeptools
	#add phantompeakqualtools (run_spp.R)
elif [ -z $sorted_bam ]; 
then
	for file in bam/*sort-bl.bam
	do
		samtools index -@ $threads -b $file
	done
	sorted_sample_list=$(ls bam/*_sort-bl.bam | tr "\n" " ")
	multiBamSummary bins --numberOfProcessors $threads -b $sorted_sample_list -o chip-qc/multibamsummary.npz 2>> chip-qc.log #deeptools
	plotPCA -in chip-qc/multibamsummary.npz -o chip-qc/PCA_readCounts_all.png -T "Principle Component Analysis" 2>> chip-qc.log #deeptools
	#add phantompeakqualtools (run_spp.R)
fi
plotCorrelation -in chip-qc/multibamsummary.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Average Scores Per Read" --whatToPlot scatterplot -o chip-qc/scatterplot_PearsonCorr_bam-Scores.png --outFileCorMatrix PearsonCorr_bam-Scores.tab #deeptools
plotCorrelation -in chip-qc/multibamsummary.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --colorMap viridis --whatToPlot heatmap -o chip-qc/heatmap_SpearmanCorr_readCounts.png --outFileCorMatrix SpearmanCorr_readCounts.tab #deeptools
