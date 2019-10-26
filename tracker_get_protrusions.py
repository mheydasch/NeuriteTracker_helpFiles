#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 12:13:04 2018

@author: max
"""

import pandas as pd
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import os
import glob
import re
import seaborn as sns
sns.set(color_codes=True)   
plt.ioff()
#%% inputs
identifier='Farred_nuclei'
Si_Well={'B2':'ARGHEF1', 'B3':'DOCK10', 'B4': 'FGD2', 'B5':'Ctrl'}
filepath='/Users/max/Desktop/Office/Phd/Data/N1E_115/SiRNA/SiRNA_21/csv/'
timeinterval=15
#set to true if you want all NaNs and 0s dropped
#set true if you want to use the robust version of z score calculation
robust=False
#%% defining loading function
length_l=[]
branch_l=[]
neurites_l=[]
wellnames=[]
#defining the patterns of the well and fov descriptors
wellpattern=re.compile('_[A-Z][0-9]_')
fovpattern=re.compile('[0-9]{4}_' + identifier)
def file_load(path, file_name, parameter, nandrop=False):
    """
    file_load(path, file_name, parameter, nandrop=False):
        path= path where the files are located
        file_name= name of the csv file to be opened
        parameter= name of the measured parameter, this will rename the value column
        nandrop= False by default, will replaces NaNs by 0s and interpolates 0 values
                if True, all nan values will be dropped from the data
        
    """
    for filepath in glob.glob(path + '*{}'.format(identifier)):
        file = os.path.join(filepath, file_name)	
        filename = vars()['filepath'].split('/')[-1] 
        print(filename)
        try: 
            temp=pd.read_csv(file)
            rows, columns = temp.shape
            timer=1 
            temp=temp.dropna(axis=0, how='all', subset=[x for x in list(temp) if x.startswith('cell ')]) 
            for x in range(len(temp)):
                temp.at[x, 'time']=timer
                timer+=1 #Adjust for respective time interval  
            well=re.search(wellpattern, file).group().strip('_')
            temp['well']=well
            fov = re.search(fovpattern, file).group().strip('_' + identifier)
            temp['fov']=fov
            temp=temp.drop(['frame'], axis=1)
            temp['filename']=filename
            length_l.append(temp) #appends current dataframe to list length_l
        except: 
            print('{} could not be found'.format(file)) 
            next 
    allwells = pd.concat(length_l, sort=True)
    allwells= pd.melt(allwells, id_vars=['well', 'fov', 'time', 'filename'])
    allwells['unique']=allwells['well']+allwells['fov']+allwells['variable']
    allwells.rename(index=str, columns={'value' : parameter}, inplace=True)
    if nandrop == False:
        allwells=allwells.fillna(0) #filling all NA values
        for cell in allwells.groupby('unique'):
            #for each parameter, for each row, if the value is 0 it will be replaced by the average of it's two
            #adjacent values    
            for i in range(0, len(cell[1])):
                if i < (len(cell[1])-1):
                    i_plus=i+1
                else:
                    i_plus=i
                if i >0:
                    i_minus=i-1
                else:
                    i_minus=i
                    #if the current value of the measurement is 0, but the two adjacent ones are not
                if cell[1][parameter][cell[1][parameter].index[i]]==0 and cell[1][parameter][cell[1][parameter].index[i_minus]]!=0 and cell[1][parameter][cell[1][parameter].index[i_plus]]!=0:
                    new_value=(cell[1][parameter][cell[1][parameter].index[i_minus]]+cell[1][parameter][cell[1][parameter].index[i_plus]])/2
                    #the current value in the dataframe (not the grouped object) will be replaced
                    #by the mean of the two adjacent values
                    allwells.loc[[cell[1][parameter].index[i]],[parameter]]=new_value
    if nandrop == True:
        allwells=allwells.dropna() #removing all lines where no value could be obtained, this could probably already be implemented at the "allwells" level
        allwells=allwells[allwells[parameter]!=0] #activate if you want to ignore zero values        
    return allwells
#%% loading data
neurite_length=file_load(filepath, 'TotalNeuritesLength.csv', 'length')  
#%% output_identification
#pattern recognition for the SiRNA experiment
outputpattern=re.compile('SiRNA_[0-9]+')
#splitting the filepath into a list based on slashes 
outputidentifier=vars()['filepath'].split('/')
#filtering the list by the outputpattern.
outputidentifier=list(filter(outputpattern.search, outputidentifier))
#%%
file_list=neurite_length['filename'].unique()
#%%
protrusions=neurite_length[neurite_length['variable']=='cell 1']
protrusions=protrusions[['time', 'filename']]
protrusions['protrusion_events']=0
#%% calculates protrusion events
index=0

# groups the dataframe by cells
for cell in neurite_length.groupby('unique'):
    #for each timepoint of the current cell, it will be calculated if the 
    #calculation of neurite length current time- neurite length next timepoint
    # is bigger than 0. As long as the next time point still falls in the range of total
    #timepoints. 
    #if this is the case the protrusion value for the file the cell belongs to
    #will be increased at the corresponding timepoint
    for time in cell[1]['time']:
        t_next=time+1
        print('calculating for timepoint: ', time, 'of cell: ', cell[1]['variable'][0], 'file: ', cell[1]['filename'][0])
        if t_next <= len(cell[1]) and cell[1].loc[:, 'length'][(cell[1].loc[:, 'time']==time)].values[0]-cell[1]['length'][(cell[1].loc[:, 'time']==t_next)].values[0] > 0:
            protrusions.loc[:, 'protrusion_events'][(protrusions.loc[:, 'filename']==cell[1].loc[:, 'filename'][0]) & (protrusions.loc[:, 'time']==time)]+=1            
            index+=1
#%% definint output folder
print(outputidentifier[0])
path=vars()['filepath'].split('/')
path=os.path.join(path[1], path[2], path[3], path[4], path[5], path[6], path[7], path[8], path[9])
path='/'+ path + '/' 
           
#%%
max_protrusions={}
#for each filename, checks if there is more than one slice with the max timepoint
#if this is the, it will check that the max_protrusion slice is neither the first nor last
#slice 
for i in protrusions.groupby('filename'):
    if len(i[1].loc[:, 'time'][i[1].loc[:, 'protrusion_events']==i[1].loc[:, 'protrusion_events'].max()])==1:
        max_prot=i[1].loc[:, 'time'][i[1].loc[:, 'protrusion_events']==i[1].loc[:, 'protrusion_events'].max()][0]
    elif i[1].loc[:, 'time'][i[1].loc[:, 'protrusion_events']==i[1].loc[:, 'protrusion_events'].max()][0]== 1 and i[1].loc[:, 'time'][i[1].loc[:, 'protrusion_events']==i[1].loc[:, 'protrusion_events'].max()][1]!= len(neurite_length.loc[:, 'time'].unique()):
        max_prot=i[1].loc[:, 'time'][i[1].loc[:, 'protrusion_events']==i[1].loc[:, 'protrusion_events'].max()][1]
    temp={i[0]: max_prot}            
    max_protrusions.update(temp)
max_protrusions=pd.DataFrame.from_dict(max_protrusions, orient='index')
max_protrusions.reset_index(level=0, inplace=True)        
max_protrusions.columns=['filename', 'slice']    
max_protrusions.to_csv('{}{}_protrusions.csv'.format(path, outputidentifier[0]))        
    

