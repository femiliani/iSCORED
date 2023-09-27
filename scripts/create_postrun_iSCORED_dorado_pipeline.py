#!/usr/bin/python3

import sys 
import os
import itertools
import glob
import gzip
import re
import statistics
import pysam
import pandas as pd

### models = dna_r9.4.1_e8_sup@v3.3
### models = dna_r10.4.1_e8.2_400bps_sup@v3.5.2
configs = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        if not line.startswith('#'):
            info, path = line.rstrip().split('=')
            configs[info.strip()] = path.strip()
            


with open(f'{configs["processdir"]}/run_script.sh', 'w') as f:
    f.write('#!/bin/bash \n\n')
    f.write(f'cd {configs["processdir"]} \n')
    f.write(f'mkdir {configs["outputdirname"]} \n')
    f.write(f'cd {configs["outputdirname"]} \n')
    outdir = f'{configs["processdir"]}/{configs["outputdirname"]}/'
    #f.write(f'touch used_fast5.txt \n')
    f.write(f'touch duplex_reads.txt \n')
    

    f.write(f'mkdir -p ./pod5 \n')

    # convert fast5
    f.write(f'pod5 convert fast5 -r -t {configs["threads"]} -o ./pod5/all_pod5.pod5  {configs["rundir"]}  \n')
    # basecall 
    f.write(f'{configs["dorado"]} basecaller {configs["modeldir"]}/{configs["basecallmodel"]} ./pod5 --modified-bases 5mCG -b 350 --min-qscore 10 --emit-moves > {configs["patientsample"]}_dorado.sam \n')
    # nanofilt 
    f.write(f'{configs["samtools"]} fastq -@ {configs["threads"]} -TMM,ML {configs["patientsample"]}_dorado.sam | {configs["nanofilt"]} --maxlength 15000 > {configs["patientsample"]}_filt.fq \n')
    # align to hg38 for methylation
    f.write(f'{configs["bwa"]} mem -x ont2d -t {configs["threads"]} -k 12 -W 12  -A 4 -B 10 -O 6 -E 3 -T 120 -C -Y {configs["hg38ref"]} {configs["patientsample"]}_filt.fq | {configs["samtools"]} sort -@ {configs["threads"]} -o {configs["patientsample"]}_hg38.bam \n')
    # porechop
    f.write(f'{configs["porechop"]} -i {configs["patientsample"]}_filt.fq -o {configs["patientsample"]}_pore.fq -t {configs["threads"]} --check_reads 2000 \n')
    # align to hg19 for cnv 
    f.write(f'{configs["bwa"]} mem -x ont2d -t {configs["threads"]} -k 12 -W 12  -A 4 -B 10 -O 6 -E 3 -T 120 {configs["hg19ref"]} {configs["patientsample"]}_pore.fq > {configs["patientsample"]}_hg19.sam \n')
    f.write(f'{configs["modifiedScripts"]}/filterAlnScoreAndQual.py -i {configs["patientsample"]}_hg19.sam  -o {configs["patientsample"]}_hg19_filt.sam -s 120 -q 1  \n' )
    # find duplex pairs 
    f.write(f'. {configs["duplex_tools"]} \n')
    f.write(f'duplex_tools pair --output_dir {configs["patientsample"]}_pairs_from_bam {configs["patientsample"]}_dorado.sam \n' )
    f.write(f'deactivate \n')
    f.write(f"cut -f2 -d \' \'  {configs['patientsample']}_pairs_from_bam/pair_ids_filtered.txt >> {outdir}/duplex_reads.txt \n")    
    
    #now process the methylation
    # first combine the bams 
    f.write(f'{configs["mbtools"]} reference-frequency {configs["patientsample"]}_hg38.bam > {configs["patientsample"]}_all_mbtools.bed \n')
    f.write(f'python {configs["modifiedScripts"]}/convert_mbtools_to_modbam2bed.py {configs["patientsample"]}_all_mbtools.bed \n')
    f.write(f'Rscript {configs["modifiedScripts"]}/methylation_classification_nanodx.R --sample={configs["patientsample"]} --out_dir {outdir}/{configs["patientsample"]}_classification/ --methylation_file={configs["patientsample"]}_all_mbtools_converted.bed --probes_file=/home/cclalien/Documents/ont-dependencies/Rapid_CNS2/top_probes_hm450.Rdata --array_file=/home/cclalien/Documents/ont-dependencies/Rapid_CNS2/HM450.hg38.manifest.gencode.v22.Rdata --training_data=/home/cclalien/Documents/ont-dependencies/Rapid_CNS2/capper_top_100k_betas_binarised.Rdata --threads={configs["threads"]} \n')



    #picard to remove duplicates 
    f.write(f'java -jar {configs["picard"]} FilterSamReads -I {configs["patientsample"]}_hg19_filt.sam -O {configs["patientsample"]}_hg19_no_duplex.sam -READ_LIST_FILE {outdir}/duplex_reads.txt -FILTER excludeReadList \n')
    f.write(f' \n')

    f.write(f'python {configs["modifiedScripts"]}/getBinCounts.py -i {configs["patientsample"]}_hg19_no_duplex.sam -c {configs["modifiedScripts"]}/hg19.chrom.sizes -b {configs["modifiedScripts"]}/bins_50k_hg19.bed -o {outdir}/{configs["patientsample"]}_hg19_50k_bin_counts.bed -s {outdir}/{configs["patientsample"]}_hg19_50k_bin_stats.txt \n' )
    f.write(f'python {configs["modifiedScripts"]}/getBinCounts.py -i {configs["patientsample"]}_hg19_no_duplex.sam -c {configs["modifiedScripts"]}/hg19.chrom.sizes -b {configs["modifiedScripts"]}/bins_5k_hg19.bed -o {outdir}/{configs["patientsample"]}_hg19_5k_bin_counts.bed -s {outdir}/{configs["patientsample"]}_hg19_5k_bin_stats.txt \n' )
    f.write(f'Rscript {configs["modifiedScripts"]}/cnvAnalysis.R {outdir}/{configs["patientsample"]}_hg19_5k_bin_counts.bed {outdir}/{configs["patientsample"]}_hg19_5k_bins_cnv {configs["modifiedScripts"]}/bins_5k_hg19_gc.txt {configs["modifiedScripts"]}/bins_5k_hg19_exclude_noY.txt \n' )   
    f.write(f'mkdir {outdir}/results_all_samples/ \n')
    f.write(f'python {configs["modifiedScripts"]}/aziz_pipeline_FEedit.py {outdir}/{configs["patientsample"]}_hg19_50k_bin_counts.bed {configs["patientsample"]} {configs["modifiedScripts"]} {outdir} \n')
    f.write(f'\n\n\n\n \n')   

