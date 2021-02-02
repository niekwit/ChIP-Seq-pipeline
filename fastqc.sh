#!/bin/bash

echo "Performing FastQC/MultiQC"
	mkdir fastqc    	
	fastqc --threads $max_threads -o fastqc/ raw-data/*fastq.gz 2>> fastqc.log
	multiqc -o "fastqc/" "fastqc/" . 2>> fastqc.log
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
