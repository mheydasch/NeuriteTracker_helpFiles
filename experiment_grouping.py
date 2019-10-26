#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 12:25:14 2018

@author: max
"""
import pandas as pd
import re
import os
import seaborn as sns
import matplotlib.pyplot as plt
#creating a dictionary, assigning SiRNAS to GEFS and GAPS

GEF_dict={'ARGHEF28': 'GEF', 'ARGHEF3': 'GEF', 'ARHGEF1':'GEF', 'TRIO':'GEF', 'DOCK9':'GEF',
      'DOCK10':'GEF', 'FGD2':'GEF', 'ITSN1':'GEF', 'VAV2':'GEF', 'ARGHAP1':'GAP',
      'ARGHAP10': 'GAP', 'ARGHAP17':'GAP', 'ARGHAP21':'GAP', 'ARGHAP26':'GAP', 
      'RHOA':'GEF', 'CDC42':'GEF', 'RAC':'GEF'}
#creating a dictionary, assigning SiRNAS to RhoGTPases
GTP_dict={'ARGHEF28': 'RHOA', 'ARGHEF3': 'RHOA', 'ARHGEF1':'RHOA', 'TRIO':'RAC', 'DOCK9':'CDC42',
      'DOCK10':'RAC', 'FGD2':'RHO', 'ITSN1':'CDC42', 'VAV2':'RAC', 'ARGHAP1':'CDC42',
      'ARGHAP10': ('RHOA', 'CDC42'), 'ARGHAP17':'CDC42', 'ARGHAP21':'RHOA', 'ARGHAP26':'RHOA',
      'RHOA':'RHOA', 'CDC42':'CDC42', 'RAC':'RAC'}
#extract all GEFs

#pattern for zscore files
file_pattern=re.compile('SiRNA_[0-9]+_z_scores.csv')
#path for all SiRNA folders
path='/Users/max/Desktop/Office/Phd/Data/N1E_115/SiRNA/'
#pattern for SiRNA folders
folder_pattern=re.compile('SiRNA_[0-9]+')
#finds all folders
dirs=os.listdir(path)
#filters out those folders that correspond to pattern
dirs=list(filter(folder_pattern.search, dirs))
#initiates variables for loop
file=[]
z_scores=[]
#joins the list of folders with the path and looks if zscore file exists
#if it does it will be read and appended to the list of z_scores.
for folder in dirs:
    file_path=os.path.join(path, folder)
    files=os.listdir(file_path)
    zscore_file=list(filter(file_pattern.search, files))
    try:
        file = os.path.join(file_path, zscore_file[0])	
    except IndexError:
        print('no file in directory: ', file_path)
        pass
    if len(file)!=0:
        temp_z=pd.read_csv(file)
        z_scores.append(temp_z)
#z_score list is concatinated to a dataframe
z_scores=pd.concat(z_scores)
z_scores=z_scores.drop(columns='Unnamed: 0')   
#fixing a naming error introduced in a previous calulation
z_scores['well'][z_scores['well']=='ARGHEF1']='ARHGEF1' 
wells=z_scores['well'].unique() 
#%% Defining lineplot function and style
sns.set_style('ticks') #removes background and adds ticks
sns.despine() #removes top and right
#initiates dash style
#dash_list=['-', '.', '']
dash_dict={'GEF': '-' , 'GAP' : '.', 'Baseline': ''}
 

def lineplot (df, measurement, wellnames, title=''):
    '''
    lineplot (df, measurement, wellnames, title='')
    df: dataframe to get the data from
    measurement: variable that should be plotted
    wellnames: Knock downs that should be plotted
    title: title of the plot
    '''
    #plt.subplots(figsize=(20, 15))
    ax=sns.lineplot(x='time', y='value', data=df[(df['variable']==measurement) & (df['well'].isin(wellnames))] , hue='well', style='Classifier', style_order=['GEF', 'GAP', 'Baseline'], )    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[1:], labels=labels[1:])
    new_title='SiRNAs'
   
# =============================================================================
#     legend=ax.get_legend()   
#     new_labels=['']
#     for well in wellnames:
#         count=0
#         #for element in the legend, the counter will increase, and if the element contains
#         #the wellname, a list will be created with the wellname +
#         #the number of observations for that well
#         for name in legend.texts:
#             count+=1
#             #string conversion is required, because the 'Text' object cannot otherwise be 
#             #searched
#             if well in str(name):
#                 #the counter already increased by one before the condition, so 
#                 #it needs to be decreased again to adress the correct element
#                 count-=1
#                 #ax.legend.texts[count].set_text('{} n: {}'.format(well, counts[well]))
#                 #print('legend.texts: ',legend.texts[count])
#                 temp=('{}; n: {}'.format(well, df[df['well']==well]['cells'].values[0]))
#                 new_labels.append(temp)
#                 break
# =============================================================================

    
    ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #setting labels
    ax.set(xlabel='Time (minutes)', ylabel='Z-Score: '+measurement)
    art=[]
    art=art.append(ax.legend)
    plt.title(title)
    plt.show()   
    ax=ax.get_figure()
    plt.close()
    return ax, art
#%%
GEFs=[]
GAPs=[]
RHOA_GEFs=[]
RHOA_GAPs=[]
CDC42_GEFs=[]
CDC42_GAPs=[]
RAC_GEFs=[]
RAC_GEFs=[]
RHOA_m=[]
RAC_m=[]
CDC42_m=[]

#for each well, except Ctrl, if the SiRNA is a GEF it is put
#in the list GEFs, if it is a GAP, it is put in the list GAPs
for well in wells[wells!='Ctrl']:
    if GTP_dict[well]=='CDC42':
        CDC42_m.append(well)
    if GTP_dict[well]=='RAC':
        RAC_m.append(well)
    if GTP_dict[well]=='RHOA':
        RHOA_m.append(well)        
    if GEF_dict[well]=='GEF':
        GEFs.append(well)
    if GEF_dict[well]=='GAP':
        GAPs.append(well)
#df.loc[df['c1'] == 'Value', 'c2']
#creating a style variable for plotting, telling if the SiRNA
#is a GEF, GAP or baseline

z_scores['Classifier']=''
z_scores['Classifier'].loc[z_scores['well'].isin(GEFs)]='GEF'
z_scores['Classifier'].loc[z_scores['well'].isin(GAPs)]='GAP'
z_scores['Classifier'].loc[z_scores['well']=='Ctrl']='Baseline'
#%% creating combined lineplots
lineplot(z_scores, 'n_neurites', wells)
lineplot(z_scores, 'length', wells)
lineplot(z_scores, 'n_branch', wells)
lineplot(z_scores, 'length/complexity', wells)
lineplot(z_scores, 'corrected_length', wells)
#%% creating RhoGTPASE subsetted lineplots
figure1, art1 = lineplot(z_scores, 'length/neurites', RAC_m+['Ctrl'], 'Rac')
figure2, art2 =lineplot(z_scores, 'length/neurites', RHOA_m+['Ctrl'], 'RhoA')
figure3, art3 =lineplot(z_scores, 'length/neurites', CDC42_m+['Ctrl'], 'CDC42')
#%%
#% saving figures
figure1.savefig('{}Rac_affectors.png'.format(path), additional_artists=art1, bbox_inches='tight', dpi=500)
figure2.savefig('{}RhoA_affectors.png'.format(path), additional_artists=art2, bbox_inches='tight', dpi=500)
figure3.savefig('{}CDC42_affectors.png'.format(path), additional_artists=art3, bbox_inches='tight', dpi=500)
#%% n branch figures
figure3, art3 = lineplot(z_scores, 'n_branch', RAC_m+['Ctrl'], 'Rac')
figure4, art4 =lineplot(z_scores, 'n_branch', RHOA_m+['Ctrl'], 'RhoA')
figure5, art5 =lineplot(z_scores, 'n_branch', CDC42_m+['Ctrl'], 'CDC42')
#% saving figures
figure3.savefig('{}Rac_branch.png'.format(path), additional_artists=art3, bbox_inches='tight', dpi=500)
figure4.savefig('{}RhoA_branch.png'.format(path), additional_artists=art4, bbox_inches='tight', dpi=500)
figure5.savefig('{}CDC42_branch.png'.format(path), additional_artists=art5, bbox_inches='tight', dpi=500)

   