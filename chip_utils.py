#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os
import subprocess
import yaml
import sys
import pkg_resources
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy.plot import gseaplot
from  builtins import any as b_any
import pysam
import re
 


def file_exists(file): #check if file exists/is not size zero
    if os.path.exists(file):
        if os.path.getsize(file) > 0:
            print("Skipping "+file+" (already exists/analysed)")
            return(True)
    else:
        return(False)
    
def rename(work_dir):
    file=open(os.path.join(work_dir,"rename.config"), "r")
    lines=file.readlines()
    count=0
    for line in lines: #removes newline characters
        lines[count]=line.replace("\n","")
        count+=1

    for line in lines:#rename files
        old_name,new_name=line.split(";")
        os.rename(os.path.join(work_dir,
                    "raw-data",
                    old_name),os.path.join(work_dir,
                                "raw-data",
                                new_name))    
    
def fastqc(work_dir,threads,file_extension):
    
    if not os.path.isdir(os.path.join(work_dir,"fastqc")) or len(os.listdir(os.path.join(work_dir,"fastqc"))) == 0:
        os.makedirs(os.path.join(work_dir,"fastqc"),
                    exist_ok=True)
        fastqc_command="fastqc --threads " + str(threads) + " --quiet -o fastqc/ raw-data/*" + file_extension
        multiqc_command=["multiqc","-o","fastqc/","fastqc/"]
        #log commands
        with open(os.path.join(work_dir,"commands.log"),"w") as file:
            file.write("FastQC: ")
            print(fastqc_command, file=file)
            file.write("MultiQC: ")
            print(*multiqc_command, sep=" ", file=file)
        print("Running FastQC on raw data")
        subprocess.run(fastqc_command, shell=True)
        print("Running MultiQC")
        subprocess.run(multiqc_command)
    else:
        print("Skipping FastQC/MultiQC (already performed)")
    
def getExtension(work_dir):
    file_list = glob.glob(os.path.join(work_dir,"raw-data","*.gz"))
    test_file = file_list[0]
    extension_index = test_file.index(".",0)
    file_extension = test_file[extension_index:]
    return(file_extension)

def write2log(work_dir,command,name):
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep =" ",file=file)

def getEND(work_dir):
    '''
    Determine whether samples are single-end of paired-end
    '''    
    
    file_list = glob.glob(os.path.join(work_dir,
                                       "raw-data",
                                       "*.gz")) 
    
    PE_tag = "R2_001.fastq.gz"
    PE = b_any(PE_tag in x for x in file_list)
    
    ##based on Illumina naming convention: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
    
    if PE == True:
        return("PE")
    else:
        return("SE")

def trim(threads, work_dir): #make compatible with both SE and PE data (should detect automatically)
    
#cap threads at 4 for trim_galore
    if int(threads) > 4:
        threads = "4"

      
    def trimPE(work_dir, threads):
        print("Trimming paired-end fastq files")
        fastq_list = glob.glob(os.path.join(work_dir,"raw-data","*R1_001.fastq.gz"))
        for read1 in fastq_list:
            out_dir = os.path.dirname(read1)
            out_dir = out_dir.replace("raw-data","trim")
            out_file1 = read1.split(".",1)[0] + "_val_1.fq.gz"
            out_file1 = os.path.basename(out_file1)
            out_file1 = os.path.join(out_dir, out_file1)
            if not file_exists(out_file1):
                read2 = read1.replace("R1","R2")
                trim_galore_command = ["trim_galore","-j", threads, "-o", 
                                       "./trim", "--paired", read1, read2]
                #log commands
                with open(work_dir+"/commands.log", "a") as file:
                    file.write("Trim Galore: ")
                    print(*trim_galore_command, sep = " ",file=file)
                subprocess.run(trim_galore_command)
                
    def trimSE(work_dir, threads):
        print("Trimming single-end fastq files")
        fastq_list = glob.glob(os.path.join(work_dir,
                                            "raw-data",
                                            "*.fastq.gz"))
        
        for file in fastq_list:
            out_file = os.path.join(work_dir,
                                    "trim_galore",
                                    file.replace(".fastq.gz", "_trimmed.fq.gz"))
            if not file_exists(out_file):
                trim_galore_command = ["trim_galore", "-j", threads, 
                                       "-o", "./trim_galore",                                                                                                                                                                                                                                                                                                
                                       file]
                #log command
                with open(os.path.join(work_dir,"commands.log"), "a") as file:
                    file.write("Trim Galore: ")
                    print(*trim_galore_command, sep = " ", file=file)
                try:
                    subprocess.run(trim_galore_command)
                except:
                    sys.exit("ERROR: trimming error. Check commands log.")
    
    #Run appropriate trim function
    if getEND(work_dir) == "PE":
        trimPE(work_dir, threads)
    elif getEND(work_dir) == "SE":
        trimSE(work_dir, threads)

def align(script_dir,work_dir, threads, align, genome):
    #read reference yaml
    with open(os.path.join(script_dir, "references.yaml")) as f:
        references = yaml.safe_load(f)
    
    index = references[align][genome]
    blacklist = references["blacklist"][genome]
    
    def hisat2(work_dir, threads, index, blacklist):
        
        if getEND(work_dir) == "SE":
            os.makedirs(os.path.join(work_dir,"bam"),
                    exist_ok=True)
            trim_list = glob.glob(os.path.join(work_dir,"trim_galore","*_trimmed.fq.gz"))
            for file in trim_list:
                out_file = os.path.basename(file).replace("_trimmed.fq.gz", "-sort-bl.bam")
                out_file = os.path.join(work_dir, "bam", out_file)
                if not file_exists(out_file):
                    hisat2_command = ["zcat", file, "|", "hisat2 -p", threads, "-x",
                                      index, "- 2>> align.log | samtools view -q 15 -F 260 -bS -@",
                                      threads, "- | bedtools intersect -v -a 'stdin' -b", 
                                      blacklist, "-nonamecheck | samtools sort -@",
                                      threads, "- >", out_file
                                      ]
                    write2log(work_dir,hisat2_command,"HISAT2 SE: ")
                    try:
                        subprocess.run(hisat2_command)
                    except:
                        sys.exit("ERROR: HISAT2 alignment error. Check commands log.")
        elif getEND(work_dir) == "PE":
            os.makedirs(os.path.join(work_dir,"bam"),
                    exist_ok=True)
            read1_trim_list = glob.glob(os.path.join(work_dir,
                                                     "trim_galore",
                                                     "*R1_001_val_1.fq.gz"))
            read2_trim_list = [i.replace("R1_001_val_1",
                                         "R2_001_val_2") for i in read1_trim_list]
            
            for read1,read2 in zip(read1_trim_list, read2_trim_list):
                out_file = os.path.basename(read1).replace("R1_001_val_1.fq.gz", "-sort-bl.bam")
                if not file_exists(out_file):
                    hisat2_command = ["hisat2 -p", threads, "-x", index, 
                                      "-1", read1, "-2", read2, "2>> align.log | samtools view -q 15 -F 260 -bS -@",
                                      threads, "- | bedtools intersect -v -a 'stdin' -b", 
                                      blacklist, "-nonamecheck | samtools sort -@",
                                      threads, "- >", out_file
                                      ]
                    write2log(work_dir,hisat2_command,"HISAT2 PE: ")
                    try:
                        subprocess.run(hisat2_command)
                    except:
                        sys.exit("ERROR: HISAT2 alignment error. Check commands log.")
   
    def bwa(work_dir, threads, index, blacklist):
        pass
    
    
    
    if align == "hisat2":
        hisat2(work_dir, threads, index, blacklist)
    elif align == "bwa":
        bwa(work_dir, threads, index, blacklist)
 
def plotAlignmentRates(work_dir):
    pass


