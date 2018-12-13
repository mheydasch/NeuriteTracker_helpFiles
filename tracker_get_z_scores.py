#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 13:30:16 2018

Reads in CSV files from the neurite_tracker output and computes individual zscores
for each timepoint, each cell and each loaded measurement

@author: max
"""

import pandas as pd
import numpy as np
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
#%%
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
            timer=0 
            temp=temp.dropna(axis=0, how='all', subset=[x for x in list(temp) if x.startswith('cell ')]) 
            for x in range(len(temp)):
                temp.at[x, 'time']=timer
                timer+=15 #Adjust for respective time interval  
            well=re.search(wellpattern, file).group().strip('_')
            temp['well']=well
            fov = re.search(fovpattern, file).group().strip('_' + identifier)
            temp['fov']=fov
            temp=temp.drop(['frame'], axis=1)
            length_l.append(temp) #appends current dataframe to list length_l
        except: 
            print('{} could not be found'.format(file)) 
            next 
    allwells = pd.concat(length_l, sort=True)
    allwells= pd.melt(allwells, id_vars=['well', 'fov', 'time'])
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
    return(allwells)
#%% loading data
neurite_length=file_load(filepath, 'TotalNeuritesLength.csv', 'length')  
branch_numbers = file_load(filepath, 'TotalNeuritesBranches.csv', 'n_branch')	
neurite_numbers = file_load(filepath, 'NumberOfNeurites.csv', 'n_neurites')
#%% merging data
def merge_df(df1, df2, parameter):
    merged=df1.merge(df2[[parameter,'unique', 'time']], on=['unique', 'time'])
    return merged
data=merge_df(neurite_length, neurite_numbers, 'n_neurites')
data=merge_df(data, branch_numbers, 'n_branch')
#%%creating new measurements based on the three loaded measurements
data['length/complexity']=data['length']/data['n_branch']
data['length/neurites']=data['length']/data['n_neurites']
data['corrected_length']=data['length/complexity']/data['n_neurites']
data['well']=data['well'].map(Si_Well)
data['branch/length']=data['n_branch']/data['length']
#replaces infinite numbers by nans
data=data.replace([np.inf, -np.inf], np.nan)
#replaces nans by 0
data=data.fillna(0)
wellnames=data['well'].unique()
#%% output_identification
#pattern recognition for the SiRNA experiment
outputpattern=re.compile('SiRNA_[0-9]+')
#splitting the filepath into a list based on slashes 
outputidentifier=vars()['filepath'].split('/')
#filtering the list by the outputpattern.
outputidentifier=list(filter(outputpattern.search, outputidentifier))
#%% z_score calculation
"""
calculation of z_scores starts here, everything above is reusable for other analysis

"""
#either computes the zscore using the MAD (robust==True), or the standard deviation(robust==False)
def MAD_robust(x):
    if robust==True:
        med=np.median(x)
        dif=[np.abs(i-med) for i in x]
        return np.median(dif)
    else:
        return np.std(x)
def Mean_robust(x):
    if robust==True:
        return np.median(x)
    else:
        return np.mean(x)
    
        
#creating the number of observations(cells) per condition
cell_observations={}
for group in data.groupby(['time','well']):
    temp={group[1]['well'].unique()[0]:len(group[1]['unique'].unique())}
    cell_observations.update(temp)    


#groups merged by well and time, and calculates the median for each parameter
mean=data.groupby(['well', 'time'], as_index=False).agg({'length': [MAD_robust, Mean_robust], 'n_neurites':[MAD_robust, Mean_robust],'n_branch':[MAD_robust, Mean_robust],'length/complexity':[MAD_robust, Mean_robust], 'length/neurites':[MAD_robust, Mean_robust], 'corrected_length':[MAD_robust, Mean_robust], 'branch/length': [MAD_robust, Mean_robust]})
#groups merged together by time, without the control and calculates the MAD and median for each parameter. All conditions combined
CTRL_MAD=data[data['well']=='Ctrl'].groupby(('time'), as_index=False).agg({'length': [MAD_robust, Mean_robust], 'n_neurites':[MAD_robust, Mean_robust],'n_branch':[MAD_robust, Mean_robust],'length/complexity':[MAD_robust, Mean_robust], 'length/neurites':[MAD_robust, Mean_robust], 'corrected_length':[MAD_robust, Mean_robust], 'branch/length': [MAD_robust, Mean_robust] })

#creates more reasonable names for the columns
mean.columns = ["_".join(x) for x in mean.columns.ravel()]
CTRL_MAD.columns = ["_".join(x) for x in CTRL_MAD.columns.ravel()]
#removes 'Ctrl' from the wellnames
wellnames_2=wellnames[wellnames!='Ctrl']
#initiates a list to store the zscore dataframes
zscore=[]
#creates a list of all parameters
parameters=list(data.columns.values)
#creates a dictionary of unwanted parameters
unwanted={'well', 'time', 'fov', 'variable', 'unique'}
#creates parameters as everythin that is not in the unwanted dictionary
parameters=[e for e in parameters if e not in unwanted]
#initates an empty list for stripped parameters
stripped_parameters=[]
#to be able to use the parameters from the list for both, mean and MAD we need to strip
#'_mean' from them 
# =============================================================================
# for i in parameters:
#     if i.endswith('Mean_robust'):
#         temp=i[:-len('Mean_robust')]
#         stripped_parameters.append(temp)
# =============================================================================
unique_ids=[]
for cell_id in data['unique'].unique():
    cell_ids=cell_id
    unique_ids.append(cell_ids)
number_of_cells=len(unique_ids)
loop_iteration=1
#%%
#groups dataframe by time and loops through all timepoints
for timepoint in data.groupby(['time'], as_index=False):
    loop_iteration=1
    print('computing timepoint ', timepoint[0])
#loops through the individual cells
    for cell_id in unique_ids:
         well=data['well'][data['unique']==cell_id].values[0]
         loop_iteration+=1
         print('computing timepoint ', timepoint[0])
         print('computing cell ', loop_iteration, 'of ', number_of_cells)
#loops through the parameters              
         for parameter in parameters:
             #creates a dataframe with the wellname, the timepoint timepoint[0], the unique_cell_ID
             # the experiment name,and the zscore for the values
             #timepoint[1] where well==1 and parameter equals item
             #MAD of control
             #calculates the zscore by using the value from the test condition (dataframe mean) and subtracting
             #from that the mean of the control. this is then divided by the standard deviation of the control
             temp=pd.DataFrame([[well, timepoint[0], parameter, 
                                 (timepoint[1][parameter][timepoint[1]['unique']==cell_id].values[0]
     -CTRL_MAD[parameter+'_Mean_robust'][CTRL_MAD['time_']==timepoint[0]].values[0])/(
        CTRL_MAD[CTRL_MAD['time_']==timepoint[0]][parameter+'_MAD_robust'].values[0]), outputidentifier[0], outputidentifier[0]+'_'+cell_id ]], 
        columns=['well', 'time', 'variable', 'value', 'experiment', 'cell_id' ])
            
             zscore.append(temp)
print('concatonating zscore df')
zscore1=pd.concat(zscore, axis=0, join='outer', ignore_index=True)   
print('Figureprinting...') 
#%% definint output folder
print(outputidentifier[0])
path=vars()['filepath'].split('/')
print (filepath)
path=os.path.join(path[1], path[2], path[3], path[4], path[5], path[6], path[7], path[8], path[9])
path='/'+ path + '/' 
print(path)
#%% saving file to csv   
zscore1.to_csv('{}{}_z_scores.csv'.format(path, outputidentifier[0]))