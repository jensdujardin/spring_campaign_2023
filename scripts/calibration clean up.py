# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 09:39:33 2022

@authors: maartend & prubbens
"""

#%% STEP 1: 'GETTING STARTED: Packages, working directory, folder structure'

'''
This code should be placed in a folder containing the associated functions-script, a metadata-file (metadata.xlsx), and raw Cytosense data in a folder called 'CSV' 
'''

''' Required packages (install using cmd / powershell)'''  
# conda install seaborn scikit-learn statsmodels numba pytables
# conda install -c conda-forge python-igraph leidenalg
# pip install scanpy
# pip install -U PhenoGraph
# pip install statannotations
# pip install scikit-bio==0.5.6

''' Setting and checking working directory'''
import os
import sys
os.chdir(os.path.abspath(os.path.dirname(sys.argv[0])))    # change the working directory to the directory containing the script file
print(os.getcwd())                                         # check the working directory 

#%% STEP 2: 'PREPROCESSING OF CYTOSENSE DATA'
import numpy as np
import pandas as pd
from os import listdir

'''Import custom functions'''
from functions import remove_cells
 
'''Loop over all files, remove beads, cell fragments and noise'''   
np.random.seed(2703)                        # Set random seed 
PATH = 'Calibration/'                               # location of raw data of FCM
FILENAMES = sorted(listdir(PATH))           # generate list of names of raw datafiles

# Create dataframe to merge calibration bead runs
calibration =  pd.DataFrame()

# Open each file, apply filters, store result
for file in FILENAMES: 
    
    # Read in file; Store original number of cells
    print(file)                                                                 # print name to track progress
    df = pd.read_csv(PATH + file, sep=';', decimal=',', index_col=0, header=0)  # open a datafile
    df['Sample'] = file

    # Define additional parameters to cell fragments and noise in accordance with Thomas et al., PLOS One, 2018 (doi: 10.1371/journal.pone.0196225)
    df['FWS Range'] = df.loc[:,'FWS Maximum'] / df.loc[:,'FWS Minimum']
    df['SWS Range'] = df.loc[:,'SWS Maximum'] / df.loc[:,'SWS Minimum']
    df['FL Red Range'] = df.loc[:,'FL Red Maximum'] / df.loc[:,'FL Red Minimum']
    df['FL Red 2 Range'] = df.loc[:,'FL Red 2 Maximum'] / df.loc[:,'FL Red 2 Minimum']
    df['FL Yellow Range'] = df.loc[:,'FL Yellow Maximum'] / df.loc[:,'FL Yellow Minimum']
    df['FL Orange Range'] = df.loc[:,'FL Orange Maximum'] / df.loc[:,'FL Orange Minimum']
    df['FL Y/R ratio'] = df.loc[:,'FL Yellow Average'] / df.loc[:,'FL Red Average']
    df['FL Y/O ratio'] = df.loc[:,'FL Yellow Average'] / df.loc[:,'FL Orange Average']
    df['FL/FWS ratio'] = df.loc[:,'FL Yellow Average'] / df.loc[:,'FWS Average']
    
    # Clean-up steps based on Thomas et al. 2018 - removal of flat curves which are typical for cell fragments
    df = remove_cells(df, 'FWS Range', 0.2, maximum=False)   
    df = remove_cells(df, 'SWS Range', 0.2, maximum=False)      
    df = remove_cells(df, 'FL Red Range', 0.2, maximum=False)
    df = remove_cells(df, 'FL Red 2 Range', 0.2, maximum=False)
    df = remove_cells(df, 'FL Yellow Range', 0.2, maximum=False)
    df = remove_cells(df, 'FL Orange Range', 0.2, maximum=False)
    
    # Remove custom variables of Thomas et al. 2018
    df = df.drop(['FL Red Range','FL Red 2 Range','FL Yellow Range','FL Orange Range'], axis = 1)
    
    # Isolate beads
    df = remove_cells(df, 'SWS Length', 10, maximum=True)
    df = remove_cells(df, 'FWS Length', 10, maximum=True)
    df = remove_cells(df, 'FWS Average', 10, maximum=False)
    df = remove_cells(df, 'FWS Fill factor', 0.45, maximum=True)
    df = remove_cells(df, 'FL Yellow Average', 100, maximum=False)
    df = remove_cells(df, 'FL Yellow Fill factor', 0.5, maximum=True)
    df = remove_cells(df, 'FL Yellow Fill factor', 0.3, maximum=False)
    df = remove_cells(df, 'FL Yellow Asymmetry', 0.1, maximum=False) 
    df = remove_cells(df, 'FL Y/O ratio', 1.5, maximum=False)
    df = remove_cells(df, 'FL Y/O ratio', 2.5, maximum=True)
    df = remove_cells(df, 'FL Y/R ratio', 5, maximum=False)
    df = remove_cells(df, 'FL Y/R ratio', 20, maximum=True)
    
    # Drop particles that contain NA values
    df = df.replace(0, np.nan)
    df.dropna(axis=0, inplace=True)

    # Merge detected beads into a global dataframe
    calibration = pd.concat([df,calibration])
    
# Store measurements of calibration beads
calibration.to_csv('beads_merged.csv')

