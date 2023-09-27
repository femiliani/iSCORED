# Running iSCORED

## Table of contents:
-  The basic premise of iSCORED pipeline
-  Running the post-run analysis
-  Running the live analysis 
-  Brief explanation of where the data ends up


### The basic premise

To run iSCORED there are two options, either running live during the sequencing, or post-run.  
For your first run we recommend starting with a small dataset so that you can ensure all the tools are installed correctly.  
  
The first step, that works for both pipelines is to create an analysis directory, this should be sample specific.  

For example `/home/user/analysis/230130_new_sample_analysis` 

The next step is to copy over the configuration file from our repository. This has to be modified for each run. This can be found:  
`/<yourpath>/iSCORED/iSCORED.config`

It explains which fields have to be changed. Please avoid adding unnecessary spaces, special characters, and new lines  
Also its best practice for paths to be absolute and end in a trailing '/' eg: `/home/your/account/` instead of `/home/your/account`   
  
  
*** WHEN RUNNING LIVE *** 
The `rundir` parameter looks for the directory where MinKNOW is saving the data, this should be the lowest directory, i.e. the one that contains the fast5 output directory, the run stats, etc. The pipeline automatically looks in that directory for the fast5 directory.  

You are now ready to run the pipelines, below is an explanation on how to run the two variants plus some additional information

## Running the post-run analysis  
The post run analysis assumes that you have all the data already generated. It will process all of the data you give it. The process for this is:
1. Make an analysis directory as above, this will be both the `processdir` and possibly the `rundir` too depending on where your fast5s are and change directory into it. eg:
```
cd /home/user/analysis/
mkdir 230130_new_sample_analysis
cd 230130_new_sample_analysis
```
Note: When you are running the post-run analysis pipeline, if the data is backed up on a different drive, e.g. an HDD you have to temporarily copy the fast5s to a directory on the computers SSD because the conversion to pod5 stalls if the data is on an HDD. The simplest solution is to set the `rundir` to be the same as the `processdir` location and do add a fast5 directory in the `processdir` that contains the fast5s 

2. Copy over the configuration file from `/<yourpath>/iSCORED/iSCORED.config` and edit the appropriate fields by running:
``` 
#you can rename it to be specific for this run, for example since we are running this on 23/01/30 i will name it 230130_iSCORED.config
cp /<yourpath>/iSCORED/iSCORED.config ./230130_iSCORED.config 
```
3. Create the analysis script by running:
```
#make sure the name of the config file at the end of this command matches what you named it 
python /<yourpath>/iSCORED/create_postrun_iSCORED_dorado_pipeline.py 230130_iSCORED.config
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
4. Go back and add the `rundir` path 
5. Create the analysis script by running:
```
#make sure the name of the config file at the end of this command matches what you named it 
python /<yourpath>/iSCORED/create_postrun_iSCORED_dorado_pipeline.py 230130_iSCORED.config
```

6. Run the analysis with: 
```
bash run_script.sh
```

7. Come back in an hour and enjoy the results 


## Brief explanation of where the data ends up
The data will all appear in the output directory (automatically called output_automated, but you can change its name in the config file) the important files are:
-  The basecalled data called `samplename_dorado.sam` this file contains the basecall information (sequence + quality), the methylation info, and the moves required for duplex_basecalling. It also contains the time the read was sequenced (important if you are trying to later separate it by time)
-  The 5k bin CNV named `samplename_hg19_5k_bins_cnv.pdf`
-  The 5k and 50k bedfiles named `samplename_hg19_5k_bin_counts.bed` and `samplename_hg19_50k_bin_counts.bed`
-  Methylation classification in a directory called `samplename_classification` the important file for us is the `samplename_calibrated_classification.tsv`
-  The output of Aziz's pipeline in a directory called `results_all_samples`
