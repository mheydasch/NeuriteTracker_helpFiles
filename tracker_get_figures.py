#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 15:24:11 2018

@author: max
Attention this script alters your data:
It will either drop all nan and 0 values from your data (nandrop=true),
or replaces them by 0 and tries to fill single timepoint 0 values by the
average of the two adjacent values


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
from statsmodels.stats.multicomp import pairwise_tukeyhsd   
plt.ioff()
#%% inputs
identifier='Farred_nuclei'
Si_Well={'B2':'ARGHEF1', 'B3':'DOCK10', 'B4': 'FGD2', 'B5':'Ctrl'}
timeinterval=15
#set to true if you want all NaNs and 0s dropped
nandrop=False
#set true if you want to use the robust version of z score calculation
robust=False

#%% loading and formatting data
length_l=[]
branch_l=[]
neurites_l=[]
wellnames=[]
#sns.set(font_scale=3) #enable for figureprinting
#plt.figure(figsize=(20,15)) #enable for figureprinting
#for filepath in glob.glob('/home/mheydasch/myimaging/N1E115/SiRNA/SiRNA_11/4days/adjusted/analyzed/csv/*{}'.format(identifier)): #enable on server
#defining the patterns of the well and fov descriptors
wellpattern=re.compile('_[A-Z][0-9]_')
fovpattern=re.compile('[0-9]{4}_' + identifier)
for filepath in glob.glob('/Users/max/Desktop/Office/Phd/Data/N1E_115/SiRNA/SiRNA_21/csv/*{}'.format(identifier)):

    file = os.path.join(filepath, 'TotalNeuritesLength.csv')	
    filename = vars()['filepath'].split('/')[-1]    
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
    #totalneuritesbranches
    branchfile = os.path.join(filepath, 'TotalNeuritesBranches.csv')	
    try: 
        temp=pd.read_csv(branchfile)
        rows, columns = temp.shape
        timer=0 
        temp=temp.dropna(axis=0, how='all', subset=[x for x in list(temp) if x.startswith('cell ')]) 
        for x in range(len(temp)):
            temp.at[x, 'time']=timer
            timer+=timeinterval #Adjust for respective time interval  
        well=re.search(wellpattern, file).group().strip('_')
        temp['well']=well
        pattern=re.compile('[0-9]{4}_' + identifier)
        fov = re.search(fovpattern, file).group().strip('_' + identifier)
        temp['fov']=fov
        temp=temp.drop(['frame'], axis=1)
        branch_l.append(temp)
    except:   
        print('{} could not be found'.format(branchfile))
        next 
    #NumberOfNeurites
    neuritesfile = os.path.join(filepath, 'NumberOfNeurites.csv')	
    try: 
        temp=pd.read_csv(neuritesfile)
        rows, columns = temp.shape
        timer=0 
        temp=temp.dropna(axis=0, how='all', subset=[x for x in list(temp) if x.startswith('cell ')]) 
        for x in range(len(temp)):
            temp.at[x, 'time']=timer
            timer+=15 #Adjust for respective time interval  
        well=re.search(wellpattern, file).group().strip('_')
        temp['well']=well
        pattern=re.compile('[0-9]{4}_' + identifier)
        fov = re.search(fovpattern, file).group().strip('_' + identifier)
        temp['fov']=fov
        temp=temp.drop(['frame'], axis=1)
        neurites_l.append(temp)
    except: 
        print('{} could not be found'.format(neuritesfile))
        next      
#total length
allwells = pd.concat(length_l)
allwells= pd.melt(allwells, id_vars=['well', 'fov', 'time'])
allwells['unique']=allwells['well']+allwells['fov']+allwells['variable']

#total branches
allbranch = pd.concat(branch_l)
allbranch= pd.melt(allbranch, id_vars=['well', 'fov', 'time'])
allbranch['unique']=allbranch['well']+allbranch['fov']+allbranch['variable']

#number of neurites
allneurites = pd.concat(neurites_l)
allneurites= pd.melt(allneurites, id_vars=['well', 'fov', 'time'])
allneurites['unique']=allneurites['well']+allneurites['fov']+allneurites['variable']

#allwells=allwells.dropna()
#allwells['value']=np.log(allwells.value)
#well_mask=allwells['well']=='A3' #creating boolean mask
#timepoint_mask=(allwells['time']==465)

normality=[]

#%% selecting data and checking for normal distribution per well
wellnames=allwells['well'].unique() #extracting unique well names 
length_l=[]
branch_l=[]
neurites_l=[]

df=pd.DataFrame({'well':[], 'n':[],'mean':[], 'STDV':[], 'normal':[],'Test':[], 'F_value':[], 'P_value':[]},)
rowindex=0
#sns.set(font_scale=3) #enable for figureprinting
#plt.figure(figsize=(1, 1)) #enable for figureprinting

#%% Decide: this one drops all NAs and all 0 values
if nandrop == True:
    for i in wellnames:
        time_wellmask=(allwells['well']==i)  #creating bolean mask
        selection =allwells[time_wellmask] #applying bolean mask to select the appropriate wells 
        cleansed=selection.dropna() #removing all lines where no value could be obtained, this could probably already be implemented at the "allwells" level
        cleansed=cleansed[cleansed.value !=0] #activate if you want to ignore zero values
        #graph=sns.distplot(cleansed['value'], hist=True, label=i)
        #data selection for branches
        time_branchmask=(allbranch['well']==i)  #creating bolean mask
        selection_branch =allbranch[time_branchmask] #applying bolean mask to select the appropriate wells 
        cleansed_branch=selection_branch.dropna() #removing all lines where no value could be obtained, this could probably already be implemented at the "allwells" level
        cleansed_branch=cleansed_branch[cleansed_branch.value !=0]
        #data selection for neurites 
        time_neuritemask=(allneurites['well']==i)  #creating bolean mask
        selection_neurites =allneurites[time_neuritemask] #applying bolean mask to select the appropriate wells 
        cleansed_neurites=selection_neurites.dropna() #removing all lines where no value could be obtained, this could probably already be implemented at the "allwells" level
        cleansed_neurites=cleansed_neurites[cleansed_neurites.value !=0]

#%% Decide: This one will replace all Nas by 0
if nandrop == False:
    for i in wellnames:
        time_wellmask=(allwells['well']==i)  #creating bolean mask
        selection =allwells[time_wellmask] #applying bolean mask to select the appropriate wells 
        cleansed=selection.fillna(0) #filling all NA values
        #graph=sns.distplot(cleansed['value'], hist=True, label=i)
        #data selection for branches
        time_branchmask=(allbranch['well']==i)  #creating bolean mask
        selection_branch =allbranch[time_branchmask] #applying bolean mask to select the appropriate wells 
        cleansed_branch=selection_branch.fillna(0) #removing all lines where no value could be obtained, this could probably already be implemented at the "allwells" level
        #data selection for neurites 
        time_neuritemask=(allneurites['well']==i)  #creating bolean mask
        selection_neurites =allneurites[time_neuritemask] #applying bolean mask to select the appropriate wells 
        cleansed_neurites=selection_neurites.fillna(0) #removing all lines where no value could be obtained, this could probably already be implemented at the "allwells" level
        
    
    #%% test for normality
        k2, p = ss.normaltest(cleansed['value'])
    #%% #writing results of test to dataframe  
        print('p =',p)
        df.loc[rowindex, 'well']=i 
        df.loc[rowindex, 'n']=len(cleansed)
        df.loc[rowindex, 'mean']=cleansed['value'].mean()
        df.loc[rowindex, 'STDV']=cleansed['value'].std()
        print('{} has {} observations'.format(i, len(cleansed)))
        print('The mean distance traveled of {} is {}'.format(i, cleansed['value'].mean()))
        print('Standard deviation of traveled distance of {} is {}'.format(i, cleansed['value'].std()))
     
        if p < 0.05 :
            print('The data of well {} is not normally distributed'.format(i))
            df.loc[len(df)-1, 'normal']='No'
            norm=0
        
        else:
            print('The data of well {} is normally distributed'.format(i)) 
            df.loc[len(df)-1, 'normal']='Normal'
            norm=1 
        length_l.append(cleansed) #appending the data of the current well to a list
        neurites_l.append(cleansed_neurites)
        branch_l.append(cleansed_branch)
        rowindex+=1
        normality.append(norm) #create a list of normality for each condition 1=true 0=false
relevant_data=pd.concat(length_l) #creating a dataframe from the list
relevant_neurites=pd.concat(neurites_l)
relevant_branch=pd.concat(branch_l)
relevant_data.rename(index=str, columns={'value' :'length'}, inplace=True)
relevant_neurites.rename(index=str, columns={'value' :'n_neurites'}, inplace=True)
relevant_branch.rename(index=str, columns={'value' :'n_branch'}, inplace=True)
wellnames= relevant_data['well'].unique() #reassigning wellnames incase one disappeared after cleansing
#%% merging data frames for n_neurites, n_branch and length
merged=relevant_data.merge(relevant_neurites[['unique', 'n_neurites', 'time']], on=['unique', 'time'])
merged=merged.merge(relevant_branch[['unique', 'n_branch', 'time']], on=['unique', 'time'])
#%%correcting for minor segmentation errors
#sets the parameters
parameters=['length', 'n_neurites', 'n_branch']
#loops through each unique cell
for cell in merged.groupby('unique'):
    indices=[]
#for each parameter, for each row, if the value is 0 it will be replaced by the average of it's two
#adjacent values    
    for measurement in parameters:
       #for each element from the 0st to the last in the length of the group
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
            if cell[1][measurement][cell[1][measurement].index[i]]==0 and cell[1][measurement][cell[1][measurement].index[i_minus]]!=0 and cell[1][measurement][cell[1][measurement].index[i_plus]]!=0:
                new_value=(cell[1][measurement][cell[1][measurement].index[i_minus]]+cell[1][measurement][cell[1][measurement].index[i_plus]])/2
                #the current value in the dataframe (not the grouped object) will be replaced
                #by the mean of the two adjacent values
                merged.loc[[cell[1][measurement].index[i]],[measurement]]=new_value
            

#%%creating new measurements based on the three loaded measurements
merged['length/complexity']=merged['length']/merged['n_branch']
merged['length/neurites']=merged['length']/merged['n_neurites']
merged['corrected_length']=merged['length/complexity']/merged['n_neurites']
merged['well']=merged['well'].map(Si_Well)
merged['branch/length']=merged['n_branch']/merged['length']
#replaces infinite numbers by nans
merged=merged.replace([np.inf, -np.inf], np.nan)
#replaces nans by 0
merged=merged.fillna(0)
wellnames=merged['well'].unique()

#%% ANOVA test uncomment if needed
# =============================================================================
# threshold=sum(normality)/len(normality) #relative normality
# samples = [condition[1] for condition in merged.groupby('well')['length']]
# if  threshold >= 1: #check if data is normal to decide for a test
#  f_val, p_val = ss.f_oneway(*samples)    
#  df.loc[:, 'Test']= 'ANOVA'
#  
# else:
#     f_val, p_val = ss.kruskal(*samples) 
#     df.loc[:, 'Test']= 'Kruskal'
# print('F value: {:.3f}, p value: {:.3f}'.format(f_val, p_val))    
# print('F value: {:.3f}, p value: {:.3f}'.format(f_val, p_val))    
# df.loc[1, 'F_value']= f_val
# df.loc[1, 'P_value']= p_val
# =============================================================================
#%% Tukey test
# =============================================================================
# for name, grouped_df in relevant_data.groupby('well'):
#     print('Well {}'.format(name), pairwise_tukeyhsd(merged['length'], merged['well']))
# =============================================================================
#%% creating a normalized dataset to compare data from multiple analysis
averaged=pd.DataFrame({'time': [], 'well':[], 'length':[], 'n_neurites' :[], 'n_branch':[], 'length/complexity':[], 'length/neurites':[], 'corrected_length': []})
grouped_df=merged.groupby(['time', 'well'])

#%% z-score calculation
outputpattern=re.compile('SiRNA_[0-9]+')
outputidentifier=vars()['filepath'].split('/')
outputidentifier=list(filter(outputpattern.search, outputidentifier))
#%%
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
for group in merged.groupby(['time','well']):
    temp={group[1]['well'].unique()[0]:len(group[1]['unique'].unique())}
    cell_observations.update(temp)    


#groups merged by well and time, and calculates the median for each parameter
mean=merged.groupby(['well', 'time'], as_index=False).agg({'length': [MAD_robust, Mean_robust], 'n_neurites':[MAD_robust, Mean_robust],'n_branch':[MAD_robust, Mean_robust],'length/complexity':[MAD_robust, Mean_robust], 'length/neurites':[MAD_robust, Mean_robust], 'corrected_length':[MAD_robust, Mean_robust], 'branch/length': [MAD_robust, Mean_robust]})
#groups merged together by time, without the control and calculates the MAD and median for each parameter. All conditions combined
CTRL_MAD=merged[merged['well']=='Ctrl'].groupby(('time'), as_index=False).agg({'length': [MAD_robust, Mean_robust], 'n_neurites':[MAD_robust, Mean_robust],'n_branch':[MAD_robust, Mean_robust],'length/complexity':[MAD_robust, Mean_robust], 'length/neurites':[MAD_robust, Mean_robust], 'corrected_length':[MAD_robust, Mean_robust], 'branch/length': [MAD_robust, Mean_robust] })

#creates more reasonable names for the columns
mean.columns = ["_".join(x) for x in mean.columns.ravel()]
CTRL_MAD.columns = ["_".join(x) for x in CTRL_MAD.columns.ravel()]
#removes 'Ctrl' from the wellnames
wellnames_2=wellnames[wellnames!='Ctrl']
#initiates a list to store the zscore dataframes
zscore=[]
#creates a list of all parameters
parameters=list(merged.columns.values)
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
for cell_id in merged['unique'].unique():
    cell_ids=cell_id
    unique_ids.append(cell_ids)
number_of_cells=len(unique_ids)
loop_iteration=1
#%%
#groups dataframe by time and loops through all timepoints
for timepoint in merged.groupby(['time'], as_index=False):
    loop_iteration=1
    print('computing timepoint ', timepoint[0])
#loops through the individual cells
    for cell_id in unique_ids:
         well=merged['well'][merged['unique']==cell_id].values[0]
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
print('Figureprintin...') 

#%% boxplot length
#plt.ylim(0, 1000)
# =============================================================================
# median=merged.groupby('well')['length'].median().values
# sns.boxplot(x=merged['well'], y=merged['length'])
# ax= sns.stripplot(x=merged['well'], y=merged['length'], color='orange', jitter=0.2, size=2.5)
# sns.set_style('ticks') #removes background and adds ticks
# sns.despine() #removes top and right
# counts=x=merged['well'].value_counts()
# nobs = [str(x) for x in counts.tolist()]
# nobs = ["n: " + i for i in nobs]
# pos = range(len(nobs)) # add number of observations (currently added to wrong box)
# for tick,label in zip(pos,ax.get_xticklabels()):
#     ax.text(pos[tick], median[tick] + median[tick]*0.2, nobs[tick],
#     horizontalalignment='center',  color='black', weight='semibold')
#   
# fig1=ax.get_figure()
# fig1.set_size_inches(12, 8) 
# plt.show()  
# plt.close()
# =============================================================================
#%%boxplot length/complexity
#plt.ylim(0, 40)
# =============================================================================
# median=merged.groupby('well')['length/complexity'].median().values
# sns.boxplot(x=merged['well'], y=merged['length/complexity'])
# ax= sns.stripplot(x=merged['well'], y=merged['length/complexity'], color='orange', jitter=0.2, size=2.5)
# sns.set_style('ticks') #removes background and adds ticks
# sns.despine() #removes top and right
# counts=x=merged['well'].value_counts()
# nobs = [str(x) for x in counts.tolist()]
# nobs = ["n: " + i for i in nobs]
# pos = range(len(nobs)) # add number of observations (currently added to wrong box)
# for tick,label in zip(pos,ax.get_xticklabels()):
#     ax.text(pos[tick], median[tick] + median[tick]*0.2, nobs[tick],
#     horizontalalignment='center',  color='black', weight='semibold')
#   
# fig2=ax.get_figure()
# fig2.set_size_inches(12, 8) 
# plt.show()  
# plt.close()
# =============================================================================
#%%boxplot corrected_length
# =============================================================================
# plt.ylim(0, 10)
# median=merged.groupby('well')['corrected_length'].median().values
# sns.boxplot(x=merged['well'], y=merged['corrected_length'])
# ax= sns.stripplot(x=merged['well'], y=merged['corrected_length'], color='orange', jitter=0.2, size=2.5)
# sns.set_style('ticks') #removes background and adds ticks
# sns.despine() #removes top and right
# counts=x=merged['well'].value_counts()
# nobs = [str(x) for x in counts.tolist()]
# nobs = ["n: " + i for i in nobs]
# pos = range(len(nobs)) # add number of observations (currently added to wrong box)
# for tick,label in zip(pos,ax.get_xticklabels()):
#     ax.text(pos[tick], median[tick] + median[tick]*0.2, nobs[tick],
#     horizontalalignment='center',  color='black', weight='semibold')
#   
# fig3=ax.get_figure()
# fig3.set_size_inches(12, 8) 
# plt.show()  
# =============================================================================
#%%boxplot length/neurite
# =============================================================================
# plt.ylim(0, 2000)
# median=merged.groupby('well')['length/neurites'].median().values
# sns.boxplot(x=merged['well'], y=merged['length/neurites'])
# ax= sns.stripplot(x=merged['well'], y=merged['length/neurites'], color='orange', jitter=0.2, size=2.5)
# sns.set_style('ticks') #removes background and adds ticks
# sns.despine() #removes top and right
# counts=x=merged['well'].value_counts()
# nobs = [str(x) for x in counts.tolist()]
# nobs = ["n: " + i for i in nobs]
# pos = range(len(nobs)) # add number of observations (currently added to wrong box)
# for tick,label in zip(pos,ax.get_xticklabels()):
#     ax.text(pos[tick], median[tick] + median[tick]*0.2, nobs[tick],
#     horizontalalignment='center',  color='black', weight='semibold')
#   
# fig4=ax.get_figure()
# fig4.set_size_inches(12, 8) 
# plt.show()  
# =============================================================================
#%%boxplot number of neurites
#plt.ylim(0, 30)
# =============================================================================
# median=merged.groupby('well')['n_neurites'].median().values
# sns.boxplot(x=merged['well'], y=merged['n_neurites'])
# ax= sns.stripplot(x=merged['well'], y=merged['n_neurites'], color='orange', jitter=0.2, size=2.5)
# sns.set_style('ticks') #removes background and adds ticks
# sns.despine() #removes top and right
# counts=x=merged['well'].value_counts()
# nobs = [str(x) for x in counts.tolist()]
# nobs = ["n: " + i for i in nobs]
# pos = range(len(nobs)) # add number of observations (currently added to wrong box)
# for tick,label in zip(pos,ax.get_xticklabels()):
#     ax.text(pos[tick], median[tick] + median[tick]*0.2, nobs[tick],
#     horizontalalignment='center',  color='black', weight='semibold')
#   
# fig5=ax.get_figure()
# fig5.set_size_inches(12, 8) 
# plt.show()  
# =============================================================================
#%%# neurite length at t=0
#plt.ylim(0, 30)
# =============================================================================
# median=merged.loc[merged['time']==0].groupby('well')['length'].median().values
# sns.boxplot(x=merged['well'], y=merged['length'].loc[merged['time']==0])
# ax= sns.stripplot(x=merged['well'], y=merged['length'], color='orange', jitter=0.2, size=2.5)
# sns.set_style('ticks') #removes background and adds ticks
# sns.despine() #removes top and right
# counts=x=merged['well'].value_counts()
# nobs = [str(x) for x in counts.tolist()]
# nobs = ["n: " + i for i in nobs]
# pos = range(len(nobs)) # add number of observations (currently added to wrong box)
# for tick,label in zip(pos,ax.get_xticklabels()):
#     ax.text(pos[tick], median[tick] + median[tick]*0.2, nobs[tick],
#     horizontalalignment='center',  color='black', weight='semibold')
# ax.set_title('time=0')  
# fig1=ax.get_figure()
# fig1.set_size_inches(12, 8) 
# plt.show()  
# plt.close()
# =============================================================================
#%%timeseries neurite length
ax=sns.lineplot(x='time', y='length', data=merged, hue='well')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[1:], labels=labels[1:])
legend=ax.get_legend()
#creates an empty dictionary
neurite_observations={}
#groups the dataframe by time and well, to get from each well, for one timepoint
#the amount of neurites, adding well and neurites to the dictionary
for group in merged.groupby(['time','well']):
    temp={group[1]['well'].unique()[0]:group[1]['n_neurites'].sum()}
    neurite_observations.update(temp)  
#legend.texts[0].set_text(' ')
#creates '' at the list new labels as number 1, as this is represents the legend title
new_labels=['']
#loops through all wellnames
for well in wellnames:
    count=0
    #for element in the legend, the counter will increase, and if the element contains
    #the wellname, a list will be created with the wellname +
    #the number of observations for that well
    for name in legend.texts:
        count+=1
        #string conversion is required, because the 'Text' object cannot otherwise be 
        #searched
        if well in str(name):
            #the counter already increased by one before the condition, so 
            #it needs to be decreased again to adress the correct element
            count-=1
            #ax.legend.texts[count].set_text('{} n: {}'.format(well, counts[well]))
            #print('legend.texts: ',legend.texts[count])
            temp=('{}; n: {}'.format(well, neurite_observations[well]))
            new_labels.append(temp)
            break
#updates the labels with the new list of labels
ax.legend(handles=handles, labels=new_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)   
art6=[]
art6.append(ax.legend)
plt.show()
fig6=ax.get_figure()
plt.close()
#%%timeseries neuritelengt/number of neurites 


ax=sns.lineplot(x='time', y='length/neurites', data=merged, hue='well')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[1:], labels=labels[1:])
legend=ax.get_legend()

#legend.texts[0].set_text(' ')
new_labels=['']  
for well in wellnames:
    count=0
    for name in legend.texts:
        count+=1       
        if well in str(name):
            count-=1
            #ax.legend.texts[count].set_text('{} n: {}'.format(well, counts[well]))
        
            temp=('{}; n: {}'.format(well, cell_observations[well]))
            new_labels.append(temp)
            break
#edits the legend, uses the just created new labels list, for labels, places the legend
            #outside of the graph
ax.legend(handles=handles, labels=new_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  
#creates a list for additional artists, so the legend will be saved as well.
art7=[]
art7.append(ax.legend)
#ax.legend()
plt.show()
fig7=ax.get_figure()
plt.close()
#%% timeseries neurite lenght/ number of neurites/ number of branches
ax=sns.lineplot(x='time', y='corrected_length', data=merged, hue='well')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[1:], labels=labels[1:])
legend=ax.get_legend()
#legend.texts[0].set_text(' ')
new_labels=['']
for well in wellnames:
    count=0
    #print('well: ', well)
    for name in legend.texts:
        #print('name: ', name)
        count+=1
        #print('count: ', count)
        if well in str(name):
            count-=1
            #ax.legend.texts[count].set_text('{} n: {}'.format(well, counts[well]))
            #print('legend.texts: ',legend.texts[count])
            #print(count)
            temp=('{}; n: {}'.format(well, cell_observations[well]))
            new_labels.append(temp)
            break
 
ax.legend(handles=handles, labels=new_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)    
art8=[]
art8.append(ax.legend)
plt.show()
fig8=ax.get_figure()

#%% Individual tracks
# =============================================================================
# plt.close()
# ax=sns.lineplot(x='time', y='length', data=merged, hue='well', size='unique')
# ax.get_legend().set_visible(False)
# ax.legend_.remove()
# Fig_it=ax.get_figure()
# Fig_it.set_size_inches(20, 20) 
# plt.show()
# =============================================================================
#%% Saving files
#outputidentifier=vars()['filepath'].split('/')[-3]
print(outputidentifier[0])
path=vars()['filepath'].split('/')
print (filepath)
path=os.path.join(path[1], path[2], path[3], path[4], path[5], path[6], path[7], path[8], path[9])
path='/'+ path + '/' 
print(path)

#%% data to save
#df.to_csv('{}{}_results.csv'.format(filepath, outputidentifier[0])) 
#%% Boxlplots   
# =============================================================================
# fig1.savefig('{}{}_neurite_length'.format(path, outputidentifier[0]))
# fig2.savefig('{}{}_neurite_length_p_complexity'.format(path, outputidentifier[0]))
# fig3.savefig('{}{}_corrected_length'.format(path, outputidentifier[0]))
# fig4.savefig('{}{}_length_p_neurite'.format(path, outputidentifier[0]))
# fig5.savefig('{}{}_neurite_number'.format(path, outputidentifier[0]))
# =============================================================================
#%% Timeseries
fig6.savefig('{}{}_length_over_time.png'.format(path, outputidentifier[0]), additional_artists=art6, bbox_inches='tight')
fig7.savefig('{}{}_length_p_neurite_over_time.png'.format(path, outputidentifier[0]), additional_artists=art7, bbox_inches='tight')
fig8.savefig('{}{}_corrected_length.png'.format(path, outputidentifier[0]), additional_artists=art8, bbox_inches='tight')
#%%    
zscore1.to_csv('{}{}_z_scores.csv'.format(path, outputidentifier[0]))


