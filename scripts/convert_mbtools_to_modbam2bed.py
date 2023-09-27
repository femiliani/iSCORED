#!/usr/bin/python

##This script takes the output of https://github.com/jts/mbtools and converts it to 
## https://github.com/epi2me-labs/modbam2bed


import pandas as pd
import sys 

df = pd.read_csv(sys.argv[1], sep='\t')

chrom = list(df['chromosome'])
start = list(df['position'])
stop = [int(x)+1 for x in df['position']]
mc5 = ['5mC']*len(chrom)
score = [1000*int(x) for x in df['modified_frequency']]
strand =  list(df['strand'])
#start
#stop
empty = ["0,0,0"]*len(chrom)
readcount = list(df['total_reads'])
freq = ['{0:.2f}'.format(int(x)*100) for x in df['modified_frequency']]

dfnew = pd.DataFrame(list(zip(chrom, start, stop, mc5, score, strand, start, stop, empty, readcount, freq)))

#name = sys.argv[1].split('/')[-1].split('.')[0]
name = sys.argv[1].split('.')[0]

dfnew.to_csv(f'{name}_converted.bed', index=False, header=None, sep='\t')