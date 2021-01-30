#!/bin/bash

#This pipeline performs ChIP-Seq analsyis for paired-end data (Niek Wit, University of Cambridge, 2021)
SCRIPT_DIR="/home/niek/Documents/scripts/ChIP-Seq-pipeline/"
PICARD="/home/niek/Documents/scripts/Picard/picard.jar"
max_threads=40 #$(nproc --all) #determines CPU thread count

usage() {                                    
	echo "Usage: $0 [ rename ] [ fastqc ] [ align <species> ] [ dedup ] [ pca ] [ bigwig ] [ peaks ] [ ngsplot ]"
	echo "rename: renames fastq files from configuration file (rename.config)"
	echo "fastqc: performs FastQC and MultiQC on fq.gz files"
	echo "align: aligns fq.gz files to selected genome (no deduplication) using HISAT2"
	echo -e "Available genomes for alignment:\n\t<human>: hg38\n\t<mouse>: mm9"
	echo "dedup: removes duplicates using PICARD"
	echo "qc: performs quality control on alignment files"
	echo "bigwig: creates bigWig files"
	echo "peaks: calls peaks with MACS2 using selected BAM files (enter samples in macs2-input.csv)"
	echo "ngsplot: generates metagene plots and heatmaps with ngsplot"
	exit 2
}

while getopts '?h' c
do
	case $c in
		h|?) usage;;
	esac
done

#renames files
if [[ $* == *"rename"* ]];
	then
		sed '1d' rename.config > raw-data/rename.config
		input="raw-data/rename.config"
		while IFS= read -r line
		do
			original_file=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into original file name
			new_file=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into new file name
			mv "raw-data/${original_file}" "raw-data/${new_file}"
		done < "$input"
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
fi

if [[ $* == *"fastqc"* ]];
then
	echo "Performing FastQC/MultiQC"
	mkdir fastqc    	
	fastqc --threads $max_threads -o fastqc/ raw-data/*fastq.gz 2>> fastqc.log
	multiqc -o "fastqc/" "fastqc/" . 2>> fastqc.log
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
fi

if [[ $* == *"align"* ]];
then
    	if [[ $* == *"human"* ]];
	then 
		index_path="/home/niek/Documents/references/hisat2-index/GRCh37-hg19ucsc/hg19-index" #hg19
		blacklist_path="/home/niek/Documents/references/blacklists/Human/hg19/wgEncodeDukeMapabilityRegionsExcludable.bed" #hg19
		touch homer_hg19 #temp file for loading correct genome into HOMER for peak annotation
		echo "Human genome (hg19) and blacklist selected"
	elif [[ $* == *"mouse"* ]];
	then	
		index_path="/home/niek/Documents/references/bowtie2-index/mm9/mm9.genome-index" #mm9
		blacklist_path="/home/niek/Documents/references/blacklists/Mouse/mm9-blacklist.bed" #mm9
		touch homer_mm9 #temp file for loading correct genome into HOMER for peak annotation
		echo "Mouse genome (mm9) and blacklist selected"
	fi
	echo "Performing alignment"
	mkdir -p {bam,trim_galore}
	for read1_fastq_gz in raw-data/*_1.fastq.gz
	do
		read2_fastq_gz="${read1_fastq_gz%_1.fastq.gz}_2.fastq.gz"
		echo "Trimming adapter sequences $read1_fastq_gz"
		trim_galore -j 4 -o ./trim_galore --paired $read1_fastq_gz $read2_fastq_gz 2>> align.log
		read1_val_1_fq_gz="${read1_fastq_gz%_1.fastq.gz}_1_val_1.fq.gz"
		read1_val_1_fq_gz="trim_galore/${read1_val_1_fq_gz##*/}"
		read2_val_2_fq_gz="${read2_fastq_gz%_2.fastq.gz}_2_val_2.fq.gz"
		read2_val_2_fq_gz="trim_galore/${read2_val_2_fq_gz##*/}"
		hisat2_output=${read1_fastq_gz%_1.fastq.gz}
		hisat2_output=${hisat2_output##*/}
		echo "Aligning reads to reference $read1_fastq_gz"
		hisat2 -p $max_threads -x $index_path -1 $read1_val_1_fq_gz -2 $read2_val_2_fq_gz 2>> align.log | samtools view -q 15 -F 260 -bS -@ $max_threads - > "bam/${hisat2_output}.bam" 2>> align.log
		echo "Removing blacklisted regions from alignment $read1_fastq_gz"
		bedtools intersect -v -a "bam/${hisat2_output}.bam" -b $blacklist_path > "bam/${hisat2_output}-bl.bam" -nonamecheck 2>> align.log
	done	
	#count mapped reads:	
	echo "Uniquely mapped read counts before deduplication:" >> mapped_read_count_no_dedup.txt	
	for file in bam/*bl.bam
	do
		count=$(samtools view -@ $max_threads -c -F 260 $file) #bit 3 and 9 SAM FLAGs
		echo "$file: $count" >> mapped_read_count_no_dedup.txt
	done
	echo "Alignment completed. Check if deduplication and/or downsampling is required."	
	#add code thats can suggest if samples need downsampling
	if [[ $# == 2 ]];
		then
			exit 0
	fi
fi

if [[ $@ == *"dedup"* ]] || [[ $@ == *"deduplication"* ]];
then
    	echo "Performing deduplication"
	for file in bam/*-bl.bam
	do
		echo "Removing duplicates $file"
		sort_output=${file%-bl.bam}
		sort_output="$sort_output-sort-bl.bam"
		java -jar $PICARD SortSam INPUT=$file OUTPUT=$sort_output SORT_ORDER=coordinate 2>> deduplication.log
		dedup_output=${file%-bl.bam}
		dedup_output="$dedup_output-dedupl-sort-bl.bam"
		java -jar $PICARD MarkDuplicates INPUT=$sort_output OUTPUT=$dedup_output REMOVE_DUPLICATES=TRUE METRICS_FILE=metric.txt 2>> deduplication.log
	done	
	#count mapped reads after deduplication:	
	echo "Uniquely mapped read counts after deduplication:" >> mapped_read_count_dedup.txt	
	for file in bam/*dedupl-sort-bl.bam
	do
		count=$(samtools view -c -F 260 $file) #bit 3 and 9 SAM FLAGs
		echo "$file: $count" >> mapped_read_count_dedup.txt
	done	
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
fi

if [[ $@ == *"bigwig"* ]];
then
    	echo "Creating BigWig files"
	mkdir bigwig	
	for file in bam/*dedupl-sort-bl.bam
	do 
		samtools index -@ $max_threads -b $file 2>> bigwig.log
		bigwig_output="${file%-dedupl-sort-bl.bam}-norm.bw"
		bigwig_output=${bigwig_output##*/}
		bamCoverage -p $max_threads --binSize 10 --normalizeUsing RPKM --extendReads 200 --effectiveGenomeSize 2827437033 -b $file -o bigwig/$bigwig_output 2>> bigwig.log
	done
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
fi

if [[ $@ == *"peaks"* ]];
then
    	echo "Calling/annotating peaks with MACS2/HOMER"
	sed '1d' macs2-input.csv > macs2-input-temp.csv #removes header from settings part
	settings=$(head -5 macs2-input-temp.csv) #the first 5 lines contain the MACS2 settings
	samples_length=`expr $(wc -l < macs2-input-temp.csv) - 6` 
	tail -n $samples_length macs2-input-temp.csv > macs2-samples-temp.csv
	mkdir peaks 
	if [[ -f homer_hg19 ]]; 
		then
			homer_genome=hg19
	elif [[ -f homer_mm_9 ]]; 
		then
			homer_genome=mm9
	fi
	input="macs2-samples-temp.csv"
	while IFS= read -r line
	do
		macs2_sample=$(echo $line | cut -d " " -f 1)
		macs2_input=$(echo $line | cut -d " " -f 2)
		macs2_output_name=$(echo $line | cut -d " " -f 3)
		macs2 callpeak -t "bam/$macs2_sample" -c "bam/$macs2_input" -n $macs2_output_name --outdir "peaks/$macs2_output_name" $settings 2>> peak.log
		annotatePeaks.pl "peaks/${macs2_output_name}_summits.bed" $homer_genome > "peaks/${macs2_output_name}_annotated_peaks.txt" 2>> peak.log #HOMER
	done < "$input"
	rm macs2-input-temp.csv macs2-samples-temp.csv
	if [[ $# == 1 ]];
		then
    			exit 0
	fi
fi

if [[ $@ == *"ngsplot"* ]];
then
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
	
	if [[ -f ngsplot.conf ]]; 
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
fi

if [[ $@ == *"qc"* ]] || [[ $@ == *"QC"* ]];
then
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
fi
