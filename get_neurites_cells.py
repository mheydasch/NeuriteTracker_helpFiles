#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 10:39:31 2018

@author: max
"""

import pandas as pd
#import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import os
import glob
import re
import seaborn as sns
sns.set(color_codes=True)
from statsmodels.stats.multicomp import pairwise_tukeyhsd   
plt.ioff()
#import math
#import paramiko
#import sys
#%% inputs
identifier='_folder'
#%% loading and formatting data
l=[]
wellnames=[]
#sns.set(font_scale=3) #enable for figureprinting
#plt.figure(figsize=(20,15)) #enable for figureprint
summary=pd.DataFrame()
cell_number=pd.DataFrame()
df_complex=pd.DataFrame()
#for filepath in glob.glob('/home/mheydasch/myimaging/N1E115/SiRNA/SiRNA_20/3_hours_timelapse_40x/adjusted/test/*_folder'): 
for filepath in glob.glob('/Users/max/Desktop/Office/Phd/Data/test/*_folder'):
   try:
    identifier=vars()['filepath'].split('/')[-1]
    identifier=identifier.strip('.ini_folder')
    file = (filepath + '/csv/2018_09_15_2018_09_15_B3_0011_0018_Farred_nuclei/TotalNeuritesLength.csv')
    complexity=(filepath + '/csv/2018_09_15_2018_09_15_B3_0011_0018_Farred_nuclei/NbBranchesPerNeuriteMean.csv')
    tempcomplex=pd.read_csv(complexity)
    tempcomplex=tempcomplex.dropna(axis=0, how='all', subset=[x for x in list(tempcomplex) if x.startswith('cell ')])
    temp=pd.read_csv(file)
    temp=temp.dropna(axis=0, how='all', subset=[x for x in list(temp) if x.startswith('cell ')])
    timer=0
    for x in range(len(temp)):
            temp.at[x, 'time']=timer
            tempcomplex.at[x, 'time']=timer
            timer+=15 #Adjust for respective time interval
    cell_cols= [col for col in temp.columns if 'cell' in col]  
    temp=temp.drop(['frame'], axis=1) 
    tempcomplex=tempcomplex.drop(['frame'], axis=1) 
    tempcomplex['identifier']=identifier 
    temp['identifier']=identifier  
    tempcomplex= pd.melt(tempcomplex, id_vars=['identifier', 'time'])           
    temp= pd.melt(temp, id_vars=['identifier', 'time'])
    for i in tempcomplex.index.values:
        tempcomplex['value'][i]=float(tempcomplex['value'][i])
    tot_length=temp['value'].sum()  
    tot_complexity=tempcomplex['value'].sum()
    #df=pd.DataFrame({identifier:[tot_length]})
    summary.loc[1, identifier]=tot_length
    cell_number.loc[2, identifier]=len(cell_cols)
    df_complex.loc[1, identifier]=tot_complexity
    #summary=pd.concat(df)
   except FileNotFoundError:
       pass
filepath=filepath.strip(identifier+'.ini_folder') 
longsummary=pd.melt(summary)
cell_number=pd.melt(cell_number)
df_complex=pd.melt(df_complex)
df_complex['value']=df_complex['value']/longsummary['value']
plt.figure(figsize=(80,48))
ax= sns.barplot(x=longsummary['variable'], y=longsummary['value'])
fig1=ax.get_figure()
fig1.savefig('{}_length_for_settings.png'.format(filepath))
plt.close()
plt.figure(figsize=(80,48))
ax_cell=sns.barplot(x=cell_number['variable'], y=cell_number['value'])  
fig2=ax_cell.get_figure()
fig2.savefig('{}_cellnumber_for_settings.png'.format(filepath))
plt.close()
plt.figure(figsize=(80,48))
ax_complex=sns.barplot(x=df_complex['variable'], y=df_complex['value'])  
fig3=ax_complex.get_figure()
fig3.savefig('{}_complexity_for_settings.png'.format(filepath))