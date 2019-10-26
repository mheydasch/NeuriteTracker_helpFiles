#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 16:59:38 2018
This Script computes the neurite length for the different conditions from the 
outut of a two Cell profiler CSV file.
@author: max
"""

# =============================================================================
# import matplotlib
# matplotlib.use('Agg')
# =============================================================================
import pandas as pd
#import numpy as np
# =============================================================================
import scipy.stats as ss
import matplotlib.pyplot as plt
# import os
# import glob
import re
# import seaborn as sns
# sns.set(color_codes=True)
from statsmodels.stats.multicomp import pairwise_tukeyhsd   
import seaborn as sns
# =============================================================================

path='/Users/max/Desktop/Office/Phd/Data/N1E_115/SiRNA/SiRNA_14/SIRNA14_objCells.csv'
data=pd.read_csv(path, skiprows=1)
expname=path.split('/')[-1].replace('.csv', '')
wellid=pd.read_csv('~/Desktop/Office/Phd/Data/N1E_115/SiRNA/SiRNA_14/SIRNA14_Image.csv')
data=data.rename(columns={'Neuron_TotalNeuriteLength_imCellsSkel': 'Total_neurite_length'})
data=data[data.Total_neurite_length != 0]
data.set_index('ImageNumber', inplace=True)
wellid=wellid[['ImageNumber', 'Metadata_Well']]
wellid.set_index('ImageNumber', inplace=True)
merged=pd.merge(data, wellid, how='outer', left_index=True, right_index=True)
merged=merged[['Metadata_Well', 'Total_neurite_length']]
means=merged.groupby('Metadata_Well')['Total_neurite_length'].mean()
median=merged.groupby('Metadata_Well')['Total_neurite_length'].median().values
merged['Metadata_Well']=merged['Metadata_Well'].map({'B2':'Ctrl', 'B3' : 'Rhoa', 'B4' : 'ARHGEF3'}) #, 'B5' : 'Ctrl'})
my_pal={'Rhoa':'r', 'ARHGEF3':'skyblue', 'ARHGEF28':'skyblue', 'Ctrl':'cyan'}
sns.boxplot(x=merged['Metadata_Well'], y=merged['Total_neurite_length'], palette=my_pal)
ax= sns.stripplot(x=merged['Metadata_Well'], y=merged['Total_neurite_length'], color='orange', jitter=0.2, size=2.5)
sns.set_style('ticks') #removes background and adds ticks
sns.despine() #removes top and right
counts=merged['Metadata_Well'].value_counts()
nobs = [str(x) for x in counts.tolist()]
nobs = ["n: " + i for i in nobs]
pos = range(len(nobs)) # add number of observations (currently aded to wrong plots)
for tick,label in zip(pos,ax.get_xticklabels()):
    ax.text(pos[tick], median[tick] + median[tick]*0.2, nobs[tick],
    horizontalalignment='center',  color='black', weight='semibold')
# =============================================================================
# plt.figure()
# ax=plt.subplot(111)
# i=1
# for group in merged.groupby('Metadata_Well'):
#     ax.boxplot(group[1]['Total_neurite_length'], positions = [i])
#     plt.xticks([i], [group[0]])
#     i+=1
# ax.set_xlim(0, 5)    
# =============================================================================
#plt.bar(merged.groupby('Metadata_Well').mean(), width=0.3, color='b') #creating bars for baseline above expression levels
#plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='on') #removing ticks

samples = [condition[1] for condition in merged.groupby('Metadata_Well')['Total_neurite_length']]
f_val, p_val = ss.f_oneway(*samples) 
fig1=plt.gcf()
for name, grouped_df in merged.groupby('Metadata_Well'):
    print('Well {}'.format(name), pairwise_tukeyhsd(merged['Total_neurite_length'], merged['Metadata_Well']))
fig1.savefig('{}{}_neurite_length.jpg'.format(path.strip(('{}{}').format(expname,'.csv')), expname),dpi=150)


