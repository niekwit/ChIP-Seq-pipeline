#!/bin/bash

samtools sort -@ 30 ../bam/H3K4me3_HeLa.bam -o ../bam/H3K4me3_HeLa_sorted.bam
samtools index -@ 30 -b /home/niek/Documents/analyses/ChIP-Seq/ChIP-Seq_vs_RNA-Seq_Brian-SET1B/bam/H3K4me3_HeLa_sorted.bam

for file in ../bam/*
do 
	file2=${file%.bam}
	$sort_extension="_sorted.bam"
	$sort_output_file=$file2$sort_extension
	samtools sort -@ 30 ../bam/"$file" -o ../bam/"$sort_output_file"
	
	



#multiBamSummary bins -b ../bam/H3K4me3_HeLa.bam ../bam/wt_negative-control_HeLa.bam ../bam/S02_ControlHypoxiaChIP-H3K4Me3_Brian.bam ../bam/S13_ControlNormoxiaInput-H3K4Me3_Brian.bam ../bam/S14_ControlHypoxiaInput-H3K4Me3_Brian.bam -o multibamsummary.npz

