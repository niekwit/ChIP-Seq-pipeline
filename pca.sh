#!/bin/bash

#samtools sort -@ 30 ../bam/H3K4me3_HeLa.bam -o ../bam/H3K4me3_HeLa_sorted.bam
#samtools index -@ 30 -b /home/niek/Documents/analyses/ChIP-Seq/ChIP-Seq_vs_RNA-Seq_Brian-SET1B/bam/H3K4me3_HeLa_sorted.bam

#this loop executes the above commands for all .bam files available
for file in ../bam/*
do 
	file2=${file%.bam} #substring removal
	sort_extension="_sorted.bam"
	sort_output_file=$file2$sort_extension
	samtools sort -@ 30 $file -o $sort_output_file
	samtools index -@ 30 -b $sort_output_file
done	

	
	



multiBamSummary bins --numberOfProcessors max -b ../bam/H3K4me3_HeLa_sorted.bam ../bam/wt_negative-control_HeLa_sorted.bam ../bam/S02_ControlHypoxiaChIP-H3K4Me3_Brian_sorted.bam ../bam/S13_ControlNormoxiaInput-H3K4Me3_Brian_sorted.bam ../bam/S14_ControlHypoxiaInput-H3K4Me3_Brian_sorted.bam ../bam/S01_ControlNormoxiaChIP-H3K4Me3_Brian_sorted.bam -o multibamsummary.npz

plotPCA -in multibamsummary.npz -o PCA_readCounts.png -T "PCA of read counts"


