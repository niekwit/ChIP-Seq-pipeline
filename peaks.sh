#!/bin/bash

echo "Calling/annotating peaks with MACS2/HOMER"
sed '1d' macs2-input.csv > macs2-input-temp.csv #removes header from settings part
settings=$(head -5 macs2-input-temp.csv) #the first 5 lines contain the MACS2 settings
samples_length=`expr $(wc -l < macs2-input-temp.csv) - 6` 
tail -n $samples_length macs2-input-temp.csv > macs2-samples-temp.csv
mkdir peaks 

#Load MACS2 settings
format=$(cat settings.yaml | shyaml get-value MACS2.format)
macs2_genome=$(cat settings.yaml | shyaml get-value MACS2.genome)
genome_size=$(cat settings.yaml | shyaml get-value MACS2.genome-size)
qvalue=$(cat settings.yaml | shyaml get-value MACS2.qvalue)
extsize=$(cat settings.yaml | shyaml get-value MACS2.extsize)

#checking MACS2 settings
if [[ *$format* != "ELAND" ]] && [[ $format != "BED" ]] && [[ $format != "ELANDMULTI" ]] && [[ $format != "ELANDEXPORT" ]] && [[ $format != "SAM" ]] && [[ $format != "BAM" ]] && [[ $format != "BOWTIE" ]] && [[ $format != "BAMPE" ]] && [[ $format != "BEDPE" ]] && [[ $format != "AUTO" ]];
	then
		echo "ERROR: invalid file format chosen."
		echo "Compatible formats: ELAND, BED, ELANDMULTI, ELANDEXPORT, SAM, BAM, BOWTIE, BAMPE or BEDPE"
		echo "Automatic detection can be set with AUTO (does not work with BEDPE or BEDPE formats)"
		exit 1
fi



exit 1


#Load Homer settings
homer_genome=$(cat settings.yaml | shyaml get-value Homer.genome)

input="macs2-samples-temp.csv"
while IFS= read -r line
do
	macs2_sample=$(echo $line | cut -d " " -f 1)
	macs2_input=$(echo $line | cut -d " " -f 2)
	macs2_output_name=$(echo $line | cut -d " " -f 3)
	macs2 callpeak -t "bam/$macs2_sample" -c "bam/$macs2_input" -n $macs2_output_name --outdir "peaks/$macs2_output_name" -g "$macs2_genome" -f "$format" --qvalue "$qvalue" --extsize "$extsize" -B 2>> peak.log
	annotatePeaks.pl "peaks/${macs2_output_name}_summits.bed" $homer_genome > "peaks/${macs2_output_name}_annotated_peaks.txt" 2>> peak.log #HOMER
done < "$input"
rm macs2-input-temp.csv macs2-samples-temp.csv

if [[ $# == 1 ]];
	then
		exit 0
fi
