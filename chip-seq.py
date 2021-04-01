#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

ap = argparse.ArgumentParser()

ap.add_argument("-r", "--rename", required=False, action='store_true',
   help="Rename fq files")
ap.add_argument("-g", "--genome", required=True,
   help="Genome build")
ap.add_argument("-t", "--threads", required=False,
   help="<INT> number of CPU threads to use")
ap.add_argument("-f", "--fastqc", required=False, action='store_true',
   help="Perform FASTQC")
ap.add_argument("-a", "--align", required=False,
   help="Align data to reference")
ap.add_argument("-d", "--deduplication", required=False, action='store_true',
   help="Perform deduplication of BAM files")
ap.add_argument("-s", "--downsample", required=False, action='store_true',
   help="Perform downsampling of BAM files")
ap.add_argument("-b", "--bigwig", required=False, action='store_true',
   help="Perform deduplication of BAM files")

args = vars(ap.parse_args())

script_dir=os.path.abspath(os.path.dirname(__file__))
current_dir=os.getcwd()

fastqc=args["rename"]
if fastqc == True:
    os.system(script_dir + "/rename.sh")

fastqc=args["fastqc"]
if fastqc == True:
    os.system(script_dir + "/fastqc.sh")

