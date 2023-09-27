#!/usr/bin/python3

import glob 
import sys 
import shutil


workdir = sys.argv[1]
part = sys.argv[2]
fast5dir = sys.argv[3]

#get the total list of basecalled fastqs
fast5s = glob.glob(f"{fast5dir}/*.fast5")

#import the list of used fastqs 
usedfast5 = []
with open(f"{workdir}/used_fast5.txt") as f:
    for line in f.readlines():
        usedfast5.append(line.rstrip())
        
new_fastq = [x for x in fast5s if x not in usedfast5]

#cp the fast5 to the new part dir 
for f in new_fastq:
    name = f.split('/')[-1]
    destination = f'{workdir}/{part}/fast5/'
    finalfile = f'{destination}/{name}'
    shutil.copy(f, finalfile)

            
with open(f"{workdir}/used_fast5.txt", 'w') as f:
    for file in fast5s:
        f.write(f"{file}\n")