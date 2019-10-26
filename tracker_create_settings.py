#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 10:48:06 2018

@author: max
"""
#'/Users/max/Desktop/Office/Phd/Analysis_Software/neuritetracker-master/trunk/settings.ini'

"""
Create a list of settings files for grid search of parameters
Change only one at a time
"""
import re
from shutil import copyfile

def create_settings(template, dict_settings, output_folder):
    """
    template: path to template file
    dict_settings: key=name of parameter as to be matched by regular expr
    value = list of values to take for the given parameter
    """
    idx_file = 1
    for param in dict_settings.keys():
        for value in dict_settings[param]:
            # Copy template file
            name_file = output_folder + 'settings' + str(idx_file) + '.ini'
            copyfile(template, name_file)
            # Store content file in a list
            with open(name_file) as file:
                content = file.readlines()
            # Find line of parameter and change it
            for line_nber, line in enumerate(content):
                if re.match(param, line): 
                    #modify the parameter
                    new_line = re.sub('= \d+(\.\d+)?;', '= '+str(value)+';', line)
                    content[line_nber] = new_line
            with open(name_file, 'w') as file:
                file.writelines(content)
            idx_file += 1
    return None






create_settings('/home/mheydasch/Scripts/neuritetracker-master1/trunk/settings.ini', 
                {'LengthThresh': [20, 8, 6, 4, 3], 'WeightThreshold': [200, 300, 320, 340, 360], 
                 'NeuriteProbabilityThresh': [0.4, 0.02, 0.015, 0.01, 0.005], 
                 'NeuritePruningLengthThsh': [20, 0.8, 0.6, 0.4, 0.1]},
                 '/home/mheydasch/myimaging/N1E115/SiRNA/SiRNA_20/3_hours_timelapse_40x/adjusted/test/')

# =============================================================================
# create_settings('/Users/max/Desktop/Office/Phd/Analysis_Software/neuritetracker-master/trunk/settings.ini',
#                 {'NeuritePruningLengthThsh': [1, 2, 3, 4, 5], 'NeuriteProbabilityThresh': [0.025, 0.05, 0.1, 0.15, 0.2], 'MaxNucleusArea': [1300, 1600, 1800, 1900, 2000],
#                 'MinNucleusArea': [50, 60, 70, 80, 90], 'SmoothingNuc':[1.0, 2.0, 3.0, 4.0, 5.0], 'MSER_MaxVariation':[0.4, 0.5, 0.6, 0.7, 0.8],
#                 'LengthThresh': [10, 15, 20, 25, 30], 'MinDistanceToBoundary': [4, 5, 6, 7, 8], 'MaxEccentricity':[0.8, 0.9, 1, 1.1, 1.2], 'MinCircularity': [0.025, 0.05, 0.1, 0.15, 0.2],
#                 'TemporalWindowSize': [2, 3, 4, 5, 6], 'SpatialWindowSize':[80, 90, 100, 110, 120], 'WeightTime': [60, 80, 100, 120, 140], 'WeightShape': [60, 70, 80, 90, 100], 'WeightThreshold': [120, 160, 200, 240, 280],
#                 'bodyMax':[11234, 16851, 22468, 28085, 44936]},
#                 '/home/mheydasch/myimaging/N1E115/SiRNA/SiRNA_20/3_hours_timelapse_40x/adjusted/test/')
# =============================================================================
# =============================================================================
# create_settings('/home/mheydasch/Scripts/neuritetracker-master1/trunk/settings.ini',
#                 {'FrangiBetaOne': [0.2,0.4,0.7,0.9], 'FrangiBetaTwo': [5, 10, 15, 20, 25], 'minimalSizeOfNeurite': [5, 10, 15, 20, 25], 'GeoDistNeuriteThresh': [0.0001, 0.0002, 0.0003, 0.0004, 0.0005],
#                  'NeuriteProbabilityThresh': [0.2, 0.4, 0.6, 0.8, 0.9], 'NeuritePruningLengthThsh' : [5, 10, 15, 20, 25], 'NeuriteStabLenghtThresh': [20, 40, 60, 80, 100], 
#                  'NeuriteWeightThresh': [200, 400, 600, 900, 1100], 'NeuriteMinTrackLength': [10, 20, 30, 40, 50]}, 
#                 '/home/mheydasch/myimaging/N1E115/SiRNA/SiRNA_20/3_hours_timelapse_40x/adjusted/test/')
# =============================================================================
