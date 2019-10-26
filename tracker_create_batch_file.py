#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 12:49:35 2019

@author: max
"""

#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 12:00:03 2018

@author: max

Gets the folders and creates the matlab call for the GC analyzer
for segmentation.
Takes as arguments the path where the folders for knockdowns are located
a precoded modality parameter, -m, to select the parameter file for segmentation 
automatically for the inputs 40x, 60x and FRET(paths to these three files are
contained in a dictionary and need to be adjusted based on the personal needs).
if the -m input begins with a / it does not use the dictionary to find the parameter
files but expects the full path as an input.
If -m is biosensors there will be no parameter file selected and the input
-s is required as the full path to the folder containing the .tif files
for shade correction. In this case no segmentation will be launched but the 
biosensor analysis. This requires that the segmentation has been launched in the same 
folder before.

Assumes following folder structure for segmentation:
    /KD1/FoV/Channels/
    /KD1/FoV2/Channels/
    /KD2/FoV/Channels/
    /KDn/FoV/Channels

Assumes following folder structure for biosensor analysis:
    
    /KD1/FoV/GrowthConeAnalyzer/movieData.mat  
    /KD2/FoV/GrowthConeAnalyzer/movieData.mat  
    
Everything after 'FoV' is hard coded.
KDs can have the following naming patterns:
  '[A-Za-z0-9]+' and !='Flatfield'  
FoV folders can have the following naming patterns:
  '[0-9]+' 

    
Currently does not support multiple cropped regions per FoV.
"""
import os
import re
import argparse




#%%
def parseArguments():
  # Define the parser and read arguments
  parser = argparse.ArgumentParser(description='Get tags from all files in a directory.')
  parser.add_argument('-d', '--dir', type=str, help='The directory where the knockdown folders are', required=True)

# =============================================================================
#   #note: need to set path to parameter file for these keywords, default file if none is given, path is a path is given, set mult channels
#   #==True if biosensor is given
#   
#   parser.add_argument('-m', '--modality', type=str, 
#                       help='Can be: 40x, 60x, FRET, biosensors, or the direct path. Loads the respective parameter file', required=True)
#   
#   parser.add_argument('-s', '--shade_correction', type=str, 
#                       help='Full path to the shade correction folder', required=False)
# =============================================================================
  args = parser.parse_args()
  return(args)
  

def get_folders(path):
    imagegroups={}
    nucleus_pattern=re.compile('^nucleus')
    body_pattern=re.compile('^body')
    filelist=[f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    for file in filelist:
        if re.search(nucleus_pattern, file) !=None:
            #print(file)
            basestring=file.replace(re.search(nucleus_pattern, file).group(), '')
            if basestring in imagegroups:
                imagegroups[basestring]=[imagegroups[basestring], os.path.join(path,file)]
            else:
                imagegroups.update({basestring:os.path.join(path,file)})
        if re.search(body_pattern, file) !=None:
            basestring=file.replace(re.search(body_pattern, file).group(), '')
            if basestring in imagegroups:
                imagegroups[basestring]=[imagegroups[basestring], os.path.join(path,file)]
            else:
                imagegroups.update({basestring:os.path.join(path,file)})
    for k in imagegroups.keys():
        imagegroups[k].sort()
    return imagegroups

# =============================================================================
#%%            
def create_matlab_call(path):
    '''
    creates the part of the final call that is unique for the file,
    adds the static part of the call, and adds the previously generated
    part of the call from create_parameter_call()
    '''
    filedict=get_folders(path)
    calls=[]
    for i in filedict.keys():
        function_call="bash_fed_fileiteration(\'{}\',\'{}\')\"".format(filedict[i][1], filedict[i][0])
               
        matlab_call= '/opt/local/MATLAB/R2016b/bin/matlab -nodisplay -nosplash -nodesktop -noFigureWindows -r \"addpath(genpath(\'/home/mheydasch/Scripts/neuritetracker-master/trunk1\'));  {}'.format(function_call)
        calls.append(matlab_call)
    return calls

def create_batch_file(dir):
      '''
      creates run files for each file that needs to be fed to the growthcone analyzer
      containing the arguments for slurm
      '''
    
      listfiles = [f for f in os.listdir(dir) if re.search('^job_\d+\.txt$', f)]
      with(open('job_files/batch.sl', 'w+')) as f:
        f.write('#!/bin/bash'+'\n')
        f.write('#SBATCH'+'\n')
        f.write('#SBATCH --array=1-{}'.format(len(listfiles))+'\n')
        f.write('#SBATCH --mem=50000'+'\n')
        f.write('#SBATCH -t 6:00:00'+'\n')
        #f.write('OIFS="$IFS"')
        #f.write('IFS=$\'\n\'')
        #f.write('#SBATCH -o job_files/out/Analyzer-%A_%a.out'+'\n')
        #f.write('#SBATCH -e job_files/err/Analyzer-%A_%a.err'+'\n')
        f.write('\n')
        #f.write('print job_${SLURM_ARRAY_TASK_ID}.txt')
        #f.write('print (cat job_${SLURM_ARRAY_TASK_ID}.txt)' )
        f.write('chmod +x job_files/job_${SLURM_ARRAY_TASK_ID}.txt' +'\n')
        f.write('job_files/job_${SLURM_ARRAY_TASK_ID}.txt' + '\n')
        #f.write('$call')
        #f.write('IFS="$OIFS"')
        
#%%        
if __name__ == '__main__':
    args=parseArguments()
    path=args.dir
    os.makedirs('job_files', exist_ok=True)
    #os.makedirs('job_files/out', exist_ok=True)
    #os.makedirs('job_files/err', exist_ok=True)
   
    for n, call in enumerate(create_matlab_call(path)):
        with(open('job_files/job_{}.txt'.format(n+1), 'w+')) as fjob:
            fjob.write(call)

      # Use all job files to create a batch file
    create_batch_file('job_files/')

    #print(args)









#%% get folder 60x
