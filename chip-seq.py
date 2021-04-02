#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import subprocess
import multiprocessing
import yaml

ap = argparse.ArgumentParser()

ap.add_argument("-r", "--rename", required=False, action='store_true',
   help="Rename fq files")
ap.add_argument("-g", "--genome", required=True,
   help="Genome build")
ap.add_argument("-t", "--threads", required=False, default=1,
   help="<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
ap.add_argument("-f", "--fastqc", required=False, action='store_true',
   help="Perform FASTQC")
ap.add_argument("-a", "--align", required=False,
   help="Align data to reference")
ap.add_argument("-d", "--deduplication", required=False, action='store_true',
   help="Perform deduplication of BAM files")
ap.add_argument("-s", "--downsample", required=False, action='store_true',
   help="Perform downsampling of BAM files")
ap.add_argument("-b", "--bigwig", required=False, action='store_true',
   help="Create BigWig files")
args = vars(ap.parse_args())

script_dir=os.path.abspath(os.path.dirname(__file__))
current_dir=os.getcwd()
picard=[line[0:] for line in subprocess.check_output("find $HOME -name picard.jar", shell=True).splitlines()]
picard=picard[0].decode("utf-8")


''''
with open(script_dir+"/settings.yaml") as file:
        settings=yaml.full_load(file)
genome= #load from yaml and not command line
'''

max_threads=multiprocessing.cpu_count()
threads=args["threads"]
if threads == "max":
    threads=max_threads

rename=args["rename"]
if rename == True:
    os.system(script_dir + "/rename.sh")

fastqc=args["fastqc"]
if fastqc == True:
    subprocess.call(["bash",script_dir + "/fastqc.sh",{threads}])

align=args["align"]
align_options=["hisat2-se","hisat2-pe","bwa-se","bwa-pe"]
if align in align_options:
    subprocess.call(["bash",script_dir + "/align.sh",{threads},{align},{genome}])
    
dedup=args["deduplication"]
if dedup == True:
    subprocess.call(["bash",script_dir + "/dedup.sh",{picard}])
