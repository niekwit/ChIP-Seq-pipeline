#!/bin/bash

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
