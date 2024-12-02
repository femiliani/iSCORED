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
            



######################
# function to write the methylation classification code 
######################

def choose_methylation(methOptions):
    chosen = methOptions.split(',')
    for moption in chosen:
        if moption == "RapidCNS":
    #nanodx 
            f.write(f'        {configs["mbtools"]} reference-frequency all_parts_hg38.bam > {configs["patientsample"]}_all_mbtools.bed \n')
            f.write(f'        python {configs["modifiedScripts"]}/convert_mbtools_to_modbam2bed.py {configs["patientsample"]}_all_mbtools.bed \n')
            f.write(f'        Rscript {configs["modifiedScripts"]}/methylation_classification_nanodx.R --sample={configs["patientsample"]} --out_dir {outdir}/{configs["patientsample"]}_classification/ --methylation_file={configs["patientsample"]}_all_mbtools_converted.bed --probes_file={configs["modifiedScripts"]}/top_probes_hm450.Rdata --array_file={configs["modifiedScripts"]}/HM450.hg38.manifest.gencode.v22.Rdata --training_data={configs["modifiedScripts"]}/capper_top_100k_betas_binarised.Rdata --threads={configs["threads"]} \n')
    
        elif moption == "Sturgeon38":
        #sturgeon hg38
            f.write(f'        samtools sort -O BAM -@ {configs["threads"]} all_parts_hg38.bam -o all_parts_hg38.sorted.bam  \n')
            f.write(f'        samtools index -@ {configs["threads"]} all_parts_hg38.sorted.bam  \n')
            f.write(f'        /home/linlab_ws/analysis/dependencies/modkit/installation/bin/modkit extract --allow-non-primary -t 32 all_parts_hg38.sorted.bam  {configs["patientsample"]}_hg38_modkit_extract.txt \n')
            f.write(f'        . ~/analysis/Francesco/sturgeon_setup/sturgeon/venv/bin/activate \n')
            f.write(f'        sturgeon inputtobed -i {configs["patientsample"]}_hg38_modkit_extract.txt -o {configs["patientsample"]}_hg38_sturg_convert -s modkit --probes-file {configs["modifiedScripts"]}/probes_hg38.bed \n')
        
            f.write(f'        sturgeon predict -i {configs["patientsample"]}_hg38_sturg_convert/merged_probes_methyl_calls.bed -o {configs["patientsample"]}_hg38_sturg_out --model-files {configs["modifiedScripts"]}/sturgeon91.zip --plot-results \n')
            f.write(f'        deactivate \n')
        
        elif moption == "SturgeonT2T":
        #sturgeon T2T
            f.write(f'        samtools sort -O BAM -@ {configs["threads"]} all_parts_t2t.bam -o all_parts_t2t.sorted.bam  \n')
            f.write(f'        samtools index -@ {configs["threads"]} all_parts_t2t.sorted.bam  \n')
            f.write(f'        /home/linlab_ws/analysis/dependencies/modkit/installation/bin/modkit extract --allow-non-primary -t 32 all_parts_t2t.sorted.bam {configs["patientsample"]}_t2t_modkit_extract.txt \n')
            f.write(f'        . ~/analysis/Francesco/sturgeon_setup/sturgeon/venv/bin/activate \n')
            f.write(f'        sturgeon inputtobed -i {configs["patientsample"]}_t2t_modkit_extract.txt -o {configs["patientsample"]}_t2t_sturg_convert -s modkit --probes-file {configs["modifiedScripts"]}/probes_chm13v2.bed \n')
        
            f.write(f'        sturgeon predict -i {configs["patientsample"]}_t2t_sturg_convert/merged_probes_methyl_calls.bed -o {configs["patientsample"]}_t2t_sturg_out --model-files {configs["modifiedScripts"]}/sturgeon91.zip --plot-results \n')
            f.write(f'        deactivate \n')




#####################
# bulk of pipeline
#####################

with open(f'{configs["processdir"]}/run_script.sh', 'w') as f:
    f.write('#!/bin/bash \n\n')
    f.write(f'cd {configs["processdir"]} \n')
    f.write(f'mkdir {configs["outputdirname"]} \n')
    f.write(f'cd {configs["outputdirname"]} \n')
    outdir = f'{configs["processdir"]}/{configs["outputdirname"]}/'
    f.write(f'mkdir readids  \n')
    f.write(f'mkdir basecalled  \n')
    f.write(f'mkdir temporary  \n')
    f.write(f'mkdir aligned  \n')
    f.write(f'mkdir cnvproc  \n')
    f.write(f'mkdir methylproc  \n')
    f.write(f'touch duplex_reads.txt \n')
    
    f.write('export start_time=$(date +%s) \n')
    currpart = 'part1'
    f.write(f'export timepoint="{currpart}" \n\n')
    
    f.write('while true; do  \n')
    f.write('    sleep 30s \n')
    f.write('    export current_time=$(date +%s)  \n')
    f.write('    export elapsed_time=$((current_time - start_time))  \n')
    
    runtime = float(configs["runtime"])*60*60 # turns from h to seconds
    runsubsets = int(configs["subsets"])
    # create the dict of parts to times
    earlykey = {k:v for k,v in zip([f'part{str(i)}' for i in range(1,runsubsets+1)], [timep*int(runtime/runsubsets) for timep in range(1,runsubsets+1)] )} 

    def early_processing(p1):
        # creates the info on which timepoint we are on and then moves to the next one 
        ppast, pcurrent, pnext = [f'part{int(p1[4:])+i}' for i in [-1, 0, 1]]

        if pcurrent == "part1":
             f.write(f'    if [ $elapsed_time -ge {earlykey[pcurrent]} ] && [ $timepoint == "{pcurrent}" ]; then  \n')
        else:
             f.write(f'    elif [ $elapsed_time -ge {earlykey[pcurrent]} ] && [ $timepoint == "{pcurrent}" ]; then  \n')
        f.write(f'        export timepoint="{pnext}" \n\n')
        
        f.write(f'        pod5 view -I --no-header -r -t 12 -o ./readids/{pcurrent}_readnames.txt {configs["rundir"]}/pod5/  \n')
        if p1 == "part1":
            f.write(f'        cp ./readids/{pcurrent}_readnames.txt ./readids/{pcurrent}_only_readnames.txt  \n')
        else: 
            f.write(f'        comm -13 <(sort ./readids/{ppast}_readnames.txt) <(sort ./readids/{pcurrent}_readnames.txt) > ./readids/{pcurrent}_only_readnames.txt  \n')
        #
        #############################
        # basecalling 
        #############################
        
        ### basecall either the full set 
        if configs["timepoints"] == "all" and configs["prebasecalled"] == "none":
            f.write(f'        {configs["dorado"]} basecaller {configs["modeldir"]}/{configs["basecallmodel"]} {configs["rundir"]}/pod5/ --modified-bases 5mCG_5hmCG -l ./readids/{pcurrent}_only_readnames.txt --min-qscore 10 --emit-moves | samtools sort -n -@ 48 -o ./basecalled/{p1}_dorado.sam \n')

        ### or basecall a subset 
        elif "," in configs["timepoints"] and configs["prebasecalled"] == "none":
            start,stop = configs["timepoints"].split(',')
            f.write(f'        {configs["dorado"]} basecaller {configs["modeldir"]}/{configs["basecallmodel"]} {configs["rundir"]}/pod5/ --modified-bases 5mCG_5hmCG -l ./readids/{pcurrent}_only_readnames.txt --min-qscore 10 --emit-moves | samtools sort -n -@ 48 -o ./basecalled/{p1}_dorado_alldata.sam \n')

        #### subset out the reads by time of sequencing 
            f.write(f'        python {configs["modifiedScripts"]}/dorado_subset_time.py ./basecalled/{p1}_dorado_alldata.sam {start} {stop} ./basecalled/{p1}_timed_readids.txt \n')
            f.write(f'        {configs["samtools"]} view -N ./basecalled/{p1}_timed_readids.txt -o ./basecalled/{p1}_dorado.sam ./basecalled/{p1}_dorado_alldata.sam \n')

        #########################
        # for the pre-basecalled data
        #########################
        elif configs["timepoints"] == "all" and configs["prebasecalled"] != "none":
            f.write(f'        cp {configs["prebasecalled"]} ./basecalled/{p1}_dorado.sam \n')

        ### or basecall a subset 
        elif "," in configs["timepoints"] and configs["prebasecalled"] != "none":
            start,stop = configs["timepoints"].split(',')

        #### subset out the reads by time of sequencing 
            f.write(f'        python {configs["modifiedScripts"]}/dorado_subset_time.py {configs["prebasecalled"]} {start} {stop} ./basecalled/{p1}_timed_readids.txt \n')
            f.write(f'        {configs["samtools"]} view -N ./basecalled/{p1}_timed_readids.txt -o ./basecalled/{p1}_dorado.sam {configs["prebasecalled"]} \n')
            
        #########################

        
        # nanofilt 
        f.write(f'        {configs["samtools"]} fastq -@ {configs["threads"]} -TMN ./basecalled/{p1}_dorado.sam | {configs["nanofilt"]} --maxlength 15000 > ./temporary/{p1}_filt.fq \n')
        
        ###################
        # Methylation alignment choice
        ###################

        if "38" in configs['methylation']:
            # align to hg38 for methylation
            f.write(f'        {configs["bwa"]} mem -x ont2d -t {configs["threads"]} -k 12 -W 12  -A 4 -B 10 -O 6 -E 3 -T 120 -C -Y {configs["hg38ref"]} ./temporary/{p1}_filt.fq | {configs["samtools"]} sort -n -@ {configs["threads"]} -o ./temporary/{p1}_hg38.bam \n')
            # fix the hg38 bam methylation 
            #/home/linlab_ws/analysis/dependencies/mbtools/target/release/mbtools reference-frequency part1_hg38_repair_5monly.bam
            f.write(f'        {configs["samtools"]} sort -n -@ {configs["threads"]} ./basecalled/{p1}_dorado.sam -o ./temporary/{p1}_dorado.nsorted.sam \n')
            f.write(f'        modkit repair -d ./temporary/{p1}_dorado.nsorted.sam -a ./temporary/{p1}_hg38.bam -o ./temporary/{p1}_hg38_repaired.bam \n')
            f.write(f'        modkit adjust-mods --convert h m ./temporary/{p1}_hg38_repaired.bam ./temporary/{p1}_hg38_5mCG.bam \n')
            f.write(f'         \n')

        if "T2T" in configs['methylation']:
            # align to hg38 for methylation
            f.write(f'        {configs["bwa"]} mem -x ont2d -t {configs["threads"]} -k 12 -W 12  -A 4 -B 10 -O 6 -E 3 -T 120 -C -Y {configs["T2Tref"]} ./temporary/{p1}_filt.fq | {configs["samtools"]} sort -n -@ {configs["threads"]} -o ./temporary/{p1}_t2t.bam \n')
            # fix the hg38 bam methylation 
            #/home/linlab_ws/analysis/dependencies/mbtools/target/release/mbtools reference-frequency part1_hg38_repair_5monly.bam
            f.write(f'        {configs["samtools"]} sort -n -@ {configs["threads"]} ./basecalled/{p1}_dorado.sam -o ./temporary/{p1}_dorado.nsorted.sam \n')
            f.write(f'        modkit repair -d ./temporary/{p1}_dorado.nsorted.sam -a ./temporary/{p1}_t2t.bam -o ./temporary/{p1}_t2t_repaired.bam \n')
            f.write(f'        modkit adjust-mods --convert h m ./temporary/{p1}_t2t_repaired.bam ./temporary/{p1}_t2t_5mCG.bam \n')
            f.write(f'         \n')
        
        ##################
        # CNV analysis 
        ##################
        
        # porechop
        f.write(f'        {configs["porechop"]}-runner.py -i ./temporary/{p1}_filt.fq -o ./temporary/{p1}_pore.fq -t {configs["threads"]} --check_reads 2000 \n')
        # align to hg19 for cnv 
        f.write(f'        {configs["bwa"]} mem -x ont2d -t {configs["threads"]} -k 12 -W 12  -A 4 -B 10 -O 6 -E 3 -T 120 {configs["hg19ref"]} ./temporary/{p1}_pore.fq >./temporary/{p1}_hg19.sam \n')
        f.write(f'        {configs["modifiedScripts"]}/filterAlnScoreAndQual.py -i ./temporary/{p1}_hg19.sam  -o ./temporary/{p1}_hg19_filt.sam -s 120 -q 1  \n' )
        # find duplex pairs 
        f.write(f'        . {configs["duplex_tools"]} \n')
        f.write(f'        duplex_tools pair --output_dir {p1}_pairs_from_bam ./basecalled/{p1}_dorado.sam \n' )
        f.write(f'        deactivate \n')
        f.write(f"        cut -f2 -d \' \'  {p1}_pairs_from_bam/pair_ids_filtered.txt >> {outdir}/duplex_reads.txt \n")
        
        f.write(f'        cd {outdir} \n\n\n')
    
    for k in earlykey.keys():
        early_processing(k)

    
    #now process the methylation
    # first combine the bams 
    if "38" in configs['methylation']:
        f.write(f'        python {configs["modifiedScripts"]}/automate_dorado_merge_240131.py {outdir} _hg38_5mCG.bam {configs["threads"]} \n')

    if "T2T" in configs['methylation']:
        f.write(f'        python {configs["modifiedScripts"]}/automate_dorado_merge_240131.py {outdir} _t2t_5mCG.bam {configs["threads"]} \n')

    ##################
    # type of methylation classification chosen
    ################## 
    
    choose_methylation(configs['methylation'])
    
    #################
    # complete the CNV analysis 
    ###############
    
    #now gather the sams
    f.write(f'        python {configs["modifiedScripts"]}/automate_dorado_merge_240131.py {outdir} _hg19_filt.sam {configs["threads"]} \n')
    #picard to remove duplicates 
    f.write(f'        java -jar {configs["picard"]} FilterSamReads -I all_parts_hg19_filt.sam -O {configs["patientsample"]}_hg19_no_duplex.sam -READ_LIST_FILE {outdir}/duplex_reads.txt -FILTER excludeReadList \n')
    f.write(f' \n')

    f.write(f'        mkdir {outdir}/results_all_samples/ \n')
    f.write(f'        python {configs["modifiedScripts"]}/getBinCounts.py -i {configs["patientsample"]}_hg19_no_duplex.sam -c {configs["modifiedScripts"]}/hg19.chrom.sizes -b {configs["modifiedScripts"]}/bins_50k_hg19.bed -o {outdir}/results_all_samples/{configs["patientsample"]}_hg19_50k_bin_counts.bed -s {outdir}results_all_samples/{configs["patientsample"]}_hg19_50k_bin_stats.txt \n' )
    f.write(f'        python {configs["modifiedScripts"]}/getBinCounts.py -i {configs["patientsample"]}_hg19_no_duplex.sam -c {configs["modifiedScripts"]}/hg19.chrom.sizes -b {configs["modifiedScripts"]}/bins_5k_hg19.bed -o {outdir}/results_all_samples/{configs["patientsample"]}_hg19_5k_bin_counts.bed -s {outdir}/results_all_samples/{configs["patientsample"]}_hg19_5k_bin_stats.txt \n' )
    f.write(f'        Rscript {configs["modifiedScripts"]}/cnvAnalysis.R {outdir}/results_all_samples/{configs["patientsample"]}_hg19_5k_bin_counts.bed {outdir}/results_all_samples/{configs["patientsample"]}_hg19_5k_bins_cnv {configs["modifiedScripts"]}/bins_5k_hg19_gc.txt {configs["modifiedScripts"]}/bins_5k_hg19_exclude_noY.txt \n' )   
    f.write(f'        python {configs["modifiedScripts"]}/aziz_pipeline_FEedit.py {outdir}/results_all_samples/{configs["patientsample"]}_hg19_50k_bin_counts.bed {configs["patientsample"]} {configs["modifiedScripts"]} {outdir} \n')
    f.write(f'\n\n \n')  

    if configs["cleanup"]  == "TRUE":
        f.write(f'        rm -rf {outdir}/temporary \n\n\n') 
    
    f.write('        break\n')

    
    
    #done 
    f.write('    fi \n')
    f.write('done \n')
