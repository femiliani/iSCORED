# Running iSCORED

## Table of contents:
-  The basic premise of iSCORED pipeline
-  Running the post-run analysis
-  Running the live analysis
-  Description of config file parameters 
-  Brief explanation of where the data ends up


### The basic premise

To run iSCORED there are three options, either running live during the sequencing, or post-run on your machine. Alternatively, you can reach out to our group and have the data processed in a clinical informatics cloud platform that will simplify your computational needs. 

For your first run we recommend starting with a small dataset so that you can ensure all the tools are installed correctly.  
  
The first step, that works for both pipelines is to create an analysis directory, this should be sample specific.  

For example `/home/user/analysis/230130_new_sample_analysis` 

The next step is to copy over the configuration file from our repository. This has to be modified for each run. This can be found:  
`/<yourpath>/iSCORED/iSCORED.config`

It explains which fields have to be changed. Please avoid adding unnecessary spaces, special characters, and new lines  
Also its best practice for paths to be absolute and end in a trailing '/' eg: `/home/your/account/` instead of `/home/your/account`   
  
  
*** WHEN RUNNING LIVE *** 
The `rundir` parameter looks for the directory where MinKNOW is saving the data, this should be the lowest directory, i.e. the one that contains the pod5 output directory, the run stats, etc. The pipeline automatically looks in that directory for the pod5 directory.  

You are now ready to run the pipelines, below is an explanation on how to run the two variants plus some additional information

*** Before starting you need to obtain the `sturgeon91.zip` and `probes_chm13v2.bed` from the sturgeon github at: https://github.com/marcpaga/sturgeon and place them in the scripts direcotry. 

## Running the post-run analysis  
The post run analysis assumes that you have all the data already generated. It will process all of the data you give it. The process for this is:
1. Make an analysis directory as above, this will be both the `processdir` and possibly the `rundir` too depending on where your pod5s are and change directory into it. eg:
```
cd /home/user/analysis/
mkdir 230130_new_sample_analysis
cd 230130_new_sample_analysis
```

2. Copy over the configuration file from `/<yourpath>/iSCORED/iSCORED.config` by running:
``` 
#you can rename it to be specific for this run, for example since we are running this on 23/01/30 i will name it 230130_iSCORED.config
cp /<yourpath>/iSCORED/iSCORED.config ./230130_iSCORED.config 
```
The rundir is not necessary in this analysis, the processdir needs to point to the directory you just created, and the prebasecalled needs to point to a BAM file that was processed with the sup model with 5mCG_5hmCG modifications. You can set subsets to 1 and runtime to 0.01 (this affects how long the script waits to starts processing data).  This mode lets you control how much data you want to process, `all` will process the entire bam file,  `10,50` will process reads produced at from minutes 10 to 50. 

3. Create the analysis script by running:
```
#make sure the name of the config file at the end of this command matches what you named it 
python /<yourpath>/iSCORED/create_iSCORED_dorado_pipeline.py 230130_iSCORED.config
```

4. Run the analysis with: 
```
bash run_script.sh
```

5. Come back and enjoy the results 


## Running the live analysis  
This assumes you are about to start generating data try to set up the things you can during the incubations just before running the flow cell   

1. Make an analysis directory as above, this will be the `processdir` and change directory into it. eg:
```
cd /home/user/analysis/
mkdir 230130_new_sample_analysis
cd 230130_new_sample_analysis
```
2. Copy over the configuration file from `/<yourpath>/iSCORED/iSCORED.config` and edit the appropriate fields by running:
``` 
#you can rename it to be specific for this run, for example since we are running this on 23/01/30 i will name it 230130_iSCORED.config
### *** NB: the rundir parameter does not exist yet, since minknow hasnt started running, you will have to go back and add it later 
cp /<yourpath>/iSCORED/iSCORED.config ./230130_iSCORED.config 
```
3. Start running the flow cell, and give minknow a few minutes to generate the `rundir` 
4. Go back and add the `rundir` path to the config file 
5. Create the analysis script by running:
```
#make sure the name of the config file at the end of this command matches what you named it 
python /<yourpath>/iSCORED/create_iSCORED_dorado_pipeline.py 230130_iSCORED.config
```

6. Run the analysis with: 
```
bash run_script.sh
```

7. Come back in an hour and enjoy the results 

## Description of config file parameters 

**rundir**: in a live processing run, this is where dorado is saving the data eg: '/nanopore/data/yourrun/yoursample/20230920_1812_P2S-00522-A_PAM29382_6aef5d61/' the script automatically looks for the pod5 directory in here. In a post-run analysis, this field is not necessary, you can use the 'none' placeholder \
\
**processdir**: this is where you want the processing to happen, and where the output files will be saved, you need to create this directory. \
\
**patientsample**: this is what you want the output to be named (best not to use patient identifiers, a codename is best. \
\
**outputdirname**: what you want the output dir to be called. \
\
**threads**: the more threads, the faster this will run. We recommend a minimum of 12 threads. Numerical, eg 12\
\
**prebasecalled**: if you are running a post-run analysis, this is the path to the BAM file from basecalling, it has to be basecalled with the sup model and 5mCG_5hmCG model, otherwise we cannot extract the methylation information. In a live run analysis, this filed is not necessary, you can use the 'none' placeholder\
\
**basecallmodel**: in a live analysis this will determine which model dorado uses, we recommend having this model downloaded in the 'modeldir' directory (parameter below) as downloading it every time will slow down the pipeline. \
\
**runtime**: in a live-run this will determine how much data will be collected for analysis, it is in hours so 20 minutes (what we use) is 0.3. If you are doing a post-run analysis, there is no need to wait for data, you can set this to 0.01. \
\
**subsets**: In a live-run you can subset data to help reduce processing bottlenecking. We generally use 4 subsets, in a 20 minute run this means it will process data in ~5 minute chunks as they are produced. In a post-run analysis, this can be set to 1. \
\
**cleanup: TRUE or FALSE, if TRUE will remove some of the intermediate files that just take up space. If you are having errors in the pipeline you can set this to FALSE and troubleshoot. \
\
**timepoints**: In a live-run analysis set this to 'all', in a post-run you can set to 'all' to process everything, or, if you are only interested in a subset of the data you can set this to '10,50' this will process data generated between minutes 10 and 50. \
\
**methylation**: which methyaltion model you want to use, options are RapidCNS, Sturgeon38, SturgeonT2T as a comma-delimited string, eg. RapidCNS,Sturgeon38 it should be apparent that the more options you choose the longer the processing will take, so we do not recommend using more than 1 in a live-run analysis, but in post-run analysis this only adds 5-10 minutes. \



## Brief explanation of where the data ends up
The data will all appear in the output directory (automatically called output_automated, but you can change its name in the config file) the important files are:
-  The basecalled data called `samplename_dorado.sam` this file contains the basecall information (sequence + quality), the methylation info, and the moves required for duplex_basecalling. It also contains the time the read was sequenced (important if you are trying to later separate it by time)
-  The 5k bin CNV named `samplename_hg19_5k_bins_cnv.pdf`
-  The 5k and 50k bedfiles named `samplename_hg19_5k_bin_counts.bed` and `samplename_hg19_50k_bin_counts.bed`
-  Methylation classification in a directory called `samplename_classification` the important file for us is the `samplename_calibrated_classification.tsv`
-  The output of Aziz's pipeline in a directory called `results_all_samples`
