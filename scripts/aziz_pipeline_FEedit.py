#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:57:15 2022

@author: Aziz
"""
### sys argvs:
# bedfile 
# sample name 
# modifiedscript dir 
# outdir 


print (">>> uploading modules... <<<")
import numpy as np
import uuid
import pandas as pd
import warnings
from scipy.stats import stats
import matplotlib.pyplot as plt
import sys
import os
from itertools import groupby
from operator import itemgetter

warnings.filterwarnings("ignore")
pd.set_option('display.max_colwidth', 0)

def squarify(fig):
    w, h = fig.get_size_inches()
    if w > h:
        t = fig.subplotpars.top
        b = fig.subplotpars.bottom
        axs = h*(t-b)
        l = (1.-axs/w)/2
        fig.subplots_adjust(left=l, right=1-l)
    else:
        t = fig.subplotpars.right
        b = fig.subplotpars.left
        axs = w*(t-b)
        l = (1.-axs/h)/2
        fig.subplots_adjust(bottom=l, top=1-l)
        
def read_bed(myBed):
    content = []
    with open(myBed)as f:
        for line in f:
            content.append(line.strip().split())
    return content

raw_data=pd.DataFrame(pd.read_csv(f'{sys.argv[3]}/Template.csv'))
raw_dataG=pd.DataFrame(read_bed(f'{sys.argv[3]}/template.bed')) # from Francesco
raw_data.columns=['chromosome', 'start', 'stop', 'reads', 'proportion','Gene']
raw_dataG.columns=['chromosome', 'start', 'stop','Gene']
raw_data['Gene']=raw_dataG.Gene
raw_data['DataSource']='LGWT'
raw_data=raw_data[raw_data.chromosome!="chrY"]
raw_data=raw_data[raw_data.chromosome!="chrX"]

raw_data = raw_data.fillna('')
geneOI=pd.DataFrame(pd.read_csv(f'{sys.argv[3]}/Gene amplification_deletion.csv'))

raw_data['DoubleGene']=""
words = [rf'\b{string}\b' for string in geneOI.Gene]
nn=raw_data[raw_data['Gene'].str.contains('|'.join(words))]
nn=list(nn.index.values)
raw_data['DoubleGene'][nn]="Important"

c5=np.percentile(np.float16(raw_data.proportion),0.2)
raw_data_new=raw_data[np.float16(raw_data.proportion)>c5]
raw_data_new['LGWT-P']=raw_data_new['proportion']
raw_data_new['LGWT-R']=raw_data_new['reads']

print("This/these important bin(s) will be missed: \n",
raw_data.loc[list(set(nn) - 
set(list(raw_data_new.index.values[raw_data_new.DoubleGene=="Important"]))),
                   ['chromosome','Gene']])

df=raw_data_new

WTs={}
ff=str(uuid.uuid4())
ff=ff[0:7]

ListWTs=[''.join(list(sys.argv[1]))]

for bb_, bb in enumerate(ListWTs):    
    bb_read=read_bed(bb)
    bb_df=pd.DataFrame(bb_read)
    bb_df.columns=['chromosome', 'start', 'stop', 'reads', 'proportion']
    bb_df=bb_df[bb_df.chromosome!="chrY"]
    bb_df=bb_df[bb_df.chromosome!="chrX"]
    bb_df['DataSource']=bb
    bb_df['LGWT-P']=raw_data['proportion']
    bb_df['LGWT-R']=raw_data['reads']
    bb_df['Gene']=raw_data['Gene']
    bb_df['DoubleGene']=raw_data['DoubleGene']
    bb_df=bb_df[np.float16(raw_data.proportion)>c5]
    WTs[bb]=bb_df
    frames = [df, bb_df]
    df = pd.concat(frames)

df['proportion'] = df['proportion'].astype(float)
df['reads'] = df['reads'].astype(float)

df.index = np.arange(len(df))

chros=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 
       'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 
       'chr14', 'chr15','chr16', 'chr17', 'chr18', 'chr19','chr20', 
       'chr21','chr22']

df['ratio']=np.nan
df['Zval']=np.nan
df['TableOI']=np.nan

ff=str(uuid.uuid4())

listOIZ=[]
base_=0
FP=[]
list_samples=list(np.unique(df.DataSource))[:-1]
for DS in ListWTs:
    print(DS)
    listOIZ=[]
    Zs=[]
    readsOI=[]
    readsOIL=[]
    TableOIs=[]
    base_=0
    for i, chro in enumerate(chros):
        raw_data=df[df.DataSource==DS]
        base_=base_+len(raw_data[raw_data.chromosome==chro])
        readsOI.extend(np.float64
                       (raw_data.reads[raw_data.chromosome==chro]))
        readsOIL.extend(np.float64
                        (raw_data[raw_data.chromosome==chro]['LGWT-R']))
        val16=np.float16(df[df.DataSource==DS]['LGWT-P'][raw_data.chromosome==chro])
        ratio_=np.array(df[df.DataSource==DS]['proportion'][raw_data.chromosome==chro])/val16
        listOIZ.extend(ratio_)
        TabOI=ratio_.copy()
        Zs.extend(stats.zscore(ratio_))
        c=stats.zscore(list(np.array(raw_data['LGWT-P']
                                 [raw_data.chromosome==chro],dtype='float')))
        idx_5=list(np.where(ratio_>5)[0])
        cons_idx=list(np.where(np.diff(idx_5)==1)[0])
        new_idx5=[]
        for k, g in groupby(enumerate(idx_5), lambda ix : ix[0] - ix[1]):
            l_cons=list(map(itemgetter(1), g))
            if len(l_cons)>1:
               new_idx5.extend(l_cons) 


        strt=sys.argv[2]+"_"+chro
        del idx_5
        if len(new_idx5)>0:
            idx_5=new_idx5
            TabOI[idx_5]=666666
            
            ggOrg=raw_data['Gene'][raw_data.chromosome==chro][idx_5[0]-10:idx_5[-1]+10]
            gg=raw_data['DoubleGene'][raw_data.chromosome==chro][idx_5[0]-10:idx_5[-1]+10]
            Pos_Idx=raw_data['start'][raw_data.chromosome==chro][idx_5[0]-10:idx_5[-1]+10]
            
            Pos_Idx=Pos_Idx.values.astype(str).astype(int)
            
            fig, ax1 = plt.subplots(figsize = (7,20));ax2 = ax1.twinx()
            
            ax1.stackplot(Pos_Idx,ratio_[idx_5[0]-10:idx_5[-1]+10],color='grey',alpha=0.3)
            im=ax1.scatter(Pos_Idx,ratio_[idx_5[0]-10:idx_5[-1]+10],
                           c=list(c[idx_5[0]-10:idx_5[-1]+10]),
                           cmap='jet',alpha=1,vmin=-2, vmax=2,s=100);
            ax2.scatter(Pos_Idx,list(stats.zscore(ratio_)
                                     [idx_5[0]-10:idx_5[-1]+10]),
                        c=list(c[idx_5[0]-10:idx_5[-1]+10]),
                        cmap='jet',alpha=0);
            if len(Pos_Idx[np.where(gg=="Important")[0]])>0:
                axT = ax1.twiny()
                axT.scatter(Pos_Idx,ratio_[idx_5[0]-10:idx_5[-1]+10],
                            c=list(c[idx_5[0]-10:idx_5[-1]+10]),
                            cmap='jet',alpha=0,vmin=-2, vmax=2,s=100);
                mask=ratio_[idx_5[0]-10:idx_5[-1]+10][list(np.where(gg=="Important")[0])]>5
                axT.set_xticks(list(np.array(list(Pos_Idx
                          [list(np.where(gg=="Important")[0])]))[mask]))
                    
                axT.set_xticklabels(list(np.array(list([ggOrg.values
                [np.where(gg=="Important")[0]]][0]))[mask]),color='Darkred',
                                        rotation=90,fontsize=4)
            ax2.set_ylabel("Z value");
            ax1.set_ylabel("Ratio")
            ax1.set_xlabel("Position")
            
            squarify(fig)
            fig.savefig(f"{sys.argv[4]}/results_all_samples/"+ff[0:7]+"_{}.pdf".format(strt))
        
        TableOIs.extend(TabOI)
        
    raw_data.ratio=listOIZ
    raw_data.Zval=Zs
    raw_data.TableOI=TableOIs

    raw_data.index = np.arange(len(raw_data))
    df1 = raw_data[['chromosome', 'start','stop','ratio','Zval','Gene','TableOI']]    
    df1=df1[raw_data.DataSource!='LGWT']
    df1=df1[raw_data.TableOI==666666]
    df2=df1[df1.ratio>=5]
    df2['Gene'] = df2['Gene'].str.replace('[','')
    df2['Gene'] = df2['Gene'].str.replace(']','')
    df2.index = np.arange(len(df2))
    geneOI=pd.DataFrame(pd.read_csv(f'{sys.argv[3]}/Gene amplification_deletion.csv'))

    df2['DoubleGene']=df2.Gene; 
    
    words = [rf'\b{string}\b' for string in geneOI.Gene]
    nn=df2[df2['Gene'].str.contains('|'.join(words))]
    nn=list(nn.index.values)
    df2['DoubleGene'][nn]="Important"
    df2=df2.round(2)
    for chro in list(df2.chromosome.unique()):
        if len(np.intersect1d(df2[df2['chromosome']==chro].start.values,
                          df2[df2['chromosome']==chro].stop.values))<1:
            df2=df2[df2['chromosome']!=chro]
            
    def select_col(x):
        c1 = 'color: #B22222; font-weight: bold;font-family:Courier New, Background-color:#DCDCDC'
        c2 = 'color: #000000;font-family:Courier New, Background-color:#F0F8FF'
        #compare columns
        mask = x['DoubleGene'] == "Important"
        #DataFrame with same index and columns names as original filled empty strings
        df1_ =  pd.DataFrame(c2, index=x.index, columns=x.columns)
        #modify values of df1 column by boolean mask
        df1_.loc[mask, 'Gene'] = c1
        return df1_

    stylish_df=df2.style.apply(select_col, axis=None)
    
    strt=sys.argv[2]
    if not os.path.exists(f"{sys.argv[4]}/results_all_samples/"):
        os.mkdir(f"{sys.argv[4]}/results_all_samples/")
        
    stylish_df.to_excel(f"{sys.argv[4]}/results_all_samples/"+ff[0:7]+'_{}.xlsx'.format(strt), 
            engine='openpyxl',
            columns=['chromosome', 'start', 'stop', 'ratio', 'Zval', 'Gene'])

    FP.append((len(df2)))
    del df1, df2



