#!/bin/bash

rename(){
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
}
