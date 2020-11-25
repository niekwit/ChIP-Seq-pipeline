#!/bin/bash

#samtools sort -@ 30 ../bam/H3K4me3_HeLa.bam -o ../bam/H3K4me3_HeLa_sorted.bam
#samtools index -@ 30 -b /home/niek/Documents/analyses/ChIP-Seq/ChIP-Seq_vs_RNA-Seq_Brian-SET1B/bam/H3K4me3_HeLa_sorted.bam

#this loop executes the above commands for all .bam files available
for file in ../bam/*
do 
	file2=${file%.bam} #substring removal of .bam 
	sort_extension="_sorted.bam"
	sort_output_file=$file2$sort_extension
	samtools sort -@ 30 $file -o $sort_output_file
	samtools index -@ 30 -b $sort_output_file
done	

sorted_sample_list=$(ls ../bam/*sorted.bam | tr "\n" " ")

multiBamSummary bins --numberOfProcessors max -b $sorted_sample_list -o multibamsummary.npz

plotPCA -in multibamsummary_all.npz -o PCA_readCounts_all.png -T "PCA of read counts"


