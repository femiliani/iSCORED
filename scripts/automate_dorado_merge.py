#!/usr/bin/python3

import glob 
import sys 
import subprocess



outputdir = sys.argv[1]
reference = sys.argv[2]
threads = sys.argv[3]


samfiles = ""
if reference == 'hg19':
    found = glob.glob(f"{outputdir}/**/*_hg19_filt.sam", recursive=True)
elif reference == 'hg38':
    found = glob.glob(f"{outputdir}/**/*_hg38.bam", recursive=True)

for i in found:
    samfiles = samfiles + " " + i
        
if reference == 'hg19':
    bashCommand = f"samtools merge -@ {threads} all_parts_hg19_filt.sam {samfiles}"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
elif reference == 'hg38':
    bashCommand = f"samtools merge -@ {threads} all_parts_hg38.bam {samfiles}"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()