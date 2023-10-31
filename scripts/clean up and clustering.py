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

'''Create folder structure'''
if not os.path.exists('CSV_cleaned/'): os.makedirs('CSV_cleaned/') 
if not os.path.exists('CSV_clustered/'): os.makedirs('CSV_clustered/')
if not os.path.exists('CSV_metaclustered/'): os.makedirs('CSV_metaclustered/')      
if not os.path.exists('Results/'): os.makedirs('Results/') 
if not os.path.exists('Figures/'): os.makedirs('Figures/') 

#%% STEP 2: 'PREPROCESSING OF CYTOSENSE DATA'
import numpy as np
import pandas as pd
import seaborn as sns
from os import listdir
import matplotlib.pyplot as plt

'''Import custom functions'''
from functions import remove_cells, remove_feat, transform_asinh
 
'''Loop over all files, remove beads, cell fragments and noise'''   
np.random.seed(2703)                        # Set random seed 
PATH = 'CSV/'                               # location of raw data of FCM
PATH_PROCESSED_FILES = 'CSV_cleaned/'       # output location for data after clean-up
FILENAMES = sorted(listdir(PATH))           # generate list of names of raw datafiles

# Create dataframe to store fraction of retained particles
retained_cells = pd.DataFrame(columns=['Fraction','Total'])

# Open each file, apply filters, store result
for file in FILENAMES: 
    
    # Read in file; Store original number of cells
    print(file)                                                                 # print name to track progress
    df = pd.read_csv(PATH + file, sep=';', decimal=',', index_col=0, header=0)  # open a datafile
    n_cells = df.shape[0]                                                       # determine the number of measurements
    
    # Drop particles that contain NA values
    df = df.replace(0, np.nan)
    df.dropna(axis=0, inplace=True)
    
    # Define additional parameters to sort out detritus, cell fragments, and electrical noise
    df['FL Red MaxMin'] = df.loc[:,'FL Red Maximum'] / df.loc[:,'FL Red Minimum']
    df['SWS/FLRtot'] = df.loc[:,'SWS Total'] / df.loc[:,'FL Red Total']
    df['FWS/FLRfill'] = df.loc[:,'FWS Total'] / df.loc[:,'FL Red Fill factor']
       
    # Clean-up steps
    df = remove_cells(df, 'FWS Length', 1.0, maximum=False)       # particles must scatter across atleast 1 micrometer
    df = remove_cells(df, 'FL Red Length', 1.0, maximum=False)    # particles must fluoresce across atleast 1 micrometer  
    df = remove_cells(df, 'FL Red MaxMin', 10, maximum=False)     # removes particles with a flat fluorescence curve of Chl-a
    df = remove_cells(df, 'SWS/FLRtot', 500, maximum=True)        # removes cell fragments with relatively high scatter, but hardly any fluorescence
    df = remove_cells(df, 'FWS/FLRfill', 10**6, maximum=True)     # removes large detritus particles (high scatter, low fluorescence fill)    

    # Remove custom clean-up variables
    df = df.drop(['FL Red MaxMin','SWS/FLRtot','FWS/FLRfill'], axis = 1)
    
    # Calculate additional custom parameters
    df['Autotrophy'] = df.loc[:,'FL Red 2 Length'] / df.loc[:,'SWS Length']     # Determines the % of each cell's SWS profile that is fluorescent
    df['Scatter'] = df.loc[:,'FWS Total'] / df.loc[:,'SWS Total']               # Determines the ratio between FWS and SWS total as a function of scattering by the cell wall
    df['O/R ratio'] = df.loc[:,'FL Orange Total'] / df.loc[:,'FL Red 2 Total']  # All possible binary ratios of fluorescence
    df['O/Y ratio'] = df.loc[:,'FL Orange Total'] / df.loc[:,'FL Yellow Total']
    df['O/B ratio'] = df.loc[:,'FL Orange Total'] / df.loc[:,'FL Red Total']
    df['Y/R ratio'] = df.loc[:,'FL Yellow Total'] / df.loc[:,'FL Red 2 Total']
    df['Y/B ratio'] = df.loc[:,'FL Yellow Total'] / df.loc[:,'FL Red Total']
    df['B/R ratio'] = df.loc[:,'FL Red Total'] / df.loc[:,'FL Red 2 Total']

    # Store the number of cells after clean-up and the % of measurements retained
    retained_cells.loc[file[(len(file)-29):len(file)],'Fraction'] = df.shape[0] / n_cells
    retained_cells.loc[file[(len(file)-29):len(file)],'Total'] = df.shape[0]

    # Random subsampling of particles
    df = df.sample(3000, replace=True)
    
    # Removing non-informative columns 
    df = df.drop(['Curvature Length','Curvature Total','Curvature Maximum',
                  'Curvature Average','Curvature Inertia','Curvature Center of gravity',
                  'Curvature Fill factor','Curvature Asymmetry','Curvature Number of cells',
                  'Curvature First','Curvature Last','Curvature Minimum','Curvature SWS covariance'], axis = 1)

    df = df.drop(['FWS First','FWS Last','SWS First','SWS Last', 'FL Yellow First', 'FL Yellow Last',
                  'FL Orange First','FL Orange Last', 'FL Red First', 'FL Red Last',
                  'FL Red 2 First', 'FL Red 2 Last'], axis = 1)
                           
    df = df.drop(['FWS Minimum','SWS Minimum','FL Yellow Minimum','FL Orange Minimum','FL Red Minimum',
                  'FL Red 2 Minimum'], axis = 1)

    df = df.drop(['FWS Asymmetry','SWS Asymmetry','FL Yellow Asymmetry','FL Orange Asymmetry','FL Red Asymmetry',
                  'FL Red 2 Asymmetry'], axis = 1)

    df = df.drop(['FWS SWS covariance','FL Yellow SWS covariance','FL Orange SWS covariance','FL Red SWS covariance', 
                  'FL Red 2 SWS covariance', 'SWS SWS covariance'], axis = 1)

    # Store a new copy of the dataframe after clean-up
    df.to_csv(PATH_PROCESSED_FILES + file[(len(file)-29):len(file)])  

# Store overview of fraction retained particles per sample
retained_cells.to_csv('Results/retained_cells.csv')

#%% STEP 3: 'CLUSTERING INDIVIDUAL SAMPLES'

'''Import all required modules'''
from scanpy.external.tl import phenograph
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.cluster import adjusted_rand_score

'''Locating the filtered data and opening the metadata'''
metadata = pd.read_excel('metadata.xlsx', index_col = 0)
PATH = 'CSV_cleaned/'
FILENAMES = sorted(listdir(PATH))

'''Creating a list of all features, without unwanted features'''
df = pd.read_csv(PATH + FILENAMES[0],index_col=0, header=0)                 # open the first datafile as example
FEATURES = remove_feat(df, ['Sample Length', 'Arrival Time'])               # list all features excluding sample length and arrival time

'''Create a list of desired k-values to iterate'''
#ITERATIONS = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]  # Insert the desired numbers of nearest neighbours (k) to iterate
ITERATIONS = [10,20,30,40,50,60,70,80,90,100] 

'''Creating an empty dataframe to store richness per sample for each k-value'''
rich = pd.DataFrame(columns=[ITERATIONS],index=FILENAMES)

'''Clustering each sample''' 

for x in ITERATIONS:                                                       # looping over each potential k-value
    
    ''' Create empty dataframes to store results '''
    df_cell_counts = pd.DataFrame()                                                     # Empty df to store cell counts per cluster per sample
    df_median_cluster_val = pd.DataFrame()                                              # Empty df to store median value for each parameter per cluster
    sfcm_meta = pd.DataFrame()                                                          # Empty df to store the metadata of each sample
    df_median_fcm_values = pd.DataFrame(index=metadata.index, columns = df.columns)     # Empty df to store the global mean of each parameter across all particles per sample

    '''Transform, normalize, and cluster each file individually using PhenoGraph '''

    if not os.path.exists('CSV_clustered/k=' + str(x) + '/'):                  # Create new folder to store clustered data if needed
        os.makedirs('CSV_clustered/k=' + str(x) + '/')    

    if not os.path.exists('Results/k=' + str(x) + '/'):                        # Create new folder to store Results output if needed
        os.makedirs('Results/k=' + str(x) + '/')        
    
    for file in FILENAMES:                                                     # Loop over all files
    
        # Read in file
        print(x); print(file)                                                                       # print k-value and filename to track progress
        df = pd.read_csv(PATH + file, index_col=0, header=0)                                        # read file
        df_median_fcm_values.loc[file[0:len(file)],:] = df.median(axis=0, numeric_only=True)        # Store median values per detector
        df_sample = transform_asinh(df, FEATURES)                                                   # Log-transform the data 
    
        # Normalize the log-transformed data (mean = 0, std = 1) of all features
        scaler = StandardScaler()
        df_sample = pd.DataFrame(scaler.fit_transform(df_sample), index = df_sample.index, columns = df_sample.columns)
    
        # Add sample name
        df['Sample'] = file

        # Identify clusters using the PhenoGraph algorithm, the output "communities" contains the cluster number of each particle
        communities, graph, Q = phenograph(df_sample.loc[:,FEATURES], k=x, primary_metric='Euclidean')
    
        # Add the assigned clusters to the dataframe
        df.loc[:,'Cluster'] = communities
    
        # Remove unassigned particles, denoted by -1
        df.drop(df.loc[df['Cluster']==-1].index, inplace=True)    
    
        # Determine cell counts per cluster per sample
        df_cell_counts = pd.concat([pd.crosstab(df.loc[:,'Sample'].str[0:len(file)],df.loc[:,'Cluster']), df_cell_counts], axis=0)
    
        # Determine median values per detector per cluster per sample
        df_median_cluster_val_sample  = df.groupby('Cluster').median()
        df_median_cluster_val_sample['Sample'] = file[0:len(file)]
        df_median_cluster_val = pd.concat([df_median_cluster_val_sample,df_median_cluster_val],ignore_index=False, axis = 0)
        
        # Save the number of clusters found in this sample in the global dataframe
        rich.loc[file,x] = len(df_median_cluster_val_sample)
        
        # Save resulting dataframe containing assigned clusters and sample name
        df.to_csv('CSV_clustered/k=' + str(x) + '/' + file[0:(len(file)-4)] + '_PhenoGraph_k=' + str(x) + '.csv')

    # Export the number of cells per cluster, the median values per cluster of all samples, and the median values of all particles per sample
    df_cell_counts.to_csv('Results/k=' + str(x) + '/PhenoGraph_pop_cell_counts.csv')
    df_median_cluster_val.to_csv('Results/k=' + str(x) + '/PhenoGraph_pop_median_values.csv')
    df_median_fcm_values.to_csv('Results/k=' + str(x) + '/Median_SFCM_values_filt.csv')   

# Export the number of clusters found in each sample per iteration of k=
rich.to_csv('Results/NrOfClusters.csv')
print('\007')                                                                   # (on Windows) play an audio alert

'''Investigating the performance using the Adjusted Rand Index (ARI)''' 
rand = pd.DataFrame(columns=[ITERATIONS],index=FILENAMES)                       # Create an empty dataframe with all filenames to store the ARI

for i in range(1,len(ITERATIONS)):                                              # Loop over all files; determine the ARI of all pairs of two neighbouring k-values
   for z in range(0,len(FILENAMES)):
       file = FILENAMES[z]   
       x1 = ITERATIONS[i-1]
       x2 = ITERATIONS[i]
       print('Calculating ARI for', file, 'k=', x1, 'vs', x2)  
       df1 = pd.read_csv('CSV_clustered/k=' + str(x1) + '/' + file[(len(file)-29):len(file)-4] + '_PhenoGraph_k=' + str(x1) + '.csv', header=0) 
       df2 = pd.read_csv('CSV_clustered/k=' + str(x2) + '/' + file[(len(file)-29):len(file)-4] + '_PhenoGraph_k=' + str(x2) + '.csv', header=0)  
       df1 = df1[['Particle ID','Cluster']]
       df2 = df2[['Particle ID','Cluster']]
       df = pd.merge(df1, df2, how='inner', on=['Particle ID'],copy=False) 
       rand.loc[FILENAMES[z],ITERATIONS[i]] = adjusted_rand_score(df['Cluster_x'],df['Cluster_y'])

rand.to_csv('Results/RandScore.csv')                                            # Export the ARI scores    

'''Visualise the results of the single-file clustering''' 
if not os.path.exists('Figures/Performance/'): os.makedirs('Figures/Performance/')                           # Create a new folder if needed

rand = pd.read_csv('Results/RandScore.csv',index_col=0, header=0)               # load the randscores
rich = pd.read_csv('Results/NrOfClusters.csv',index_col=0, header=0)            # load the diversity    
rich.drop(columns=rich.columns[0], axis=1, inplace=True) 
rand.drop(columns=rand.columns[0], axis=1, inplace=True)                          

output = plt.figure()                                                           # plot the number of clusters found in each sample per iteration of k=
rich.boxplot(return_type='axes', rot=45, grid=False) 
plt.xlabel('k =', size=16)
plt.ylabel('# clusters', size=16)
output.savefig('Figures/Performance/Clustering_div.png', bbox_inches='tight', dpi=500)  

output = plt.figure()                                                           # plot the ARI of all pairs of k-values
rand.boxplot(return_type='axes', rot=45, grid=False) 
plt.xlabel('k =', size=16)
plt.ylabel('ARI', size=16)
output.savefig('Figures/Performance/Clustering_ARI.png', bbox_inches='tight', dpi=500)  

# Calculate and plot the derivatives of species richness and the Adjusted rand scores
rand = pd.read_csv('Results/RandScore.csv',index_col=0, header=0)                   # load the randscores
rich = pd.read_csv('Results/NrOfClusters.csv',index_col=0, header=0)                # load the diversity    
dARI = pd.DataFrame(columns=rand.columns[1:len(rand.columns)],index=rand.index)    
dDiv = pd.DataFrame(columns=rich.columns[1:len(rand.columns)],index=rich.index) 

for i in range(1,len(ITERATIONS)):
    for z in range(0,len(dARI.index)):
        dARI.loc[dARI.index[z],rand.columns[i]] = (rand.loc[dARI.index[z],rand.columns[i]] - rand.loc[dARI.index[z],rand.columns[i-1]]) / (int(rand.columns[i])-int(rand.columns[i-1]))
        dDiv.loc[dDiv.index[z],rich.columns[i]] = (rich.loc[dDiv.index[z],rich.columns[i]] - rich.loc[dDiv.index[z],rich.columns[i-1]]) / (int(rich.columns[i])-int(rich.columns[i-1]))

dDivmelt = pd.melt(dDiv.T, var_name= 'sample', value_name='dDIV', ignore_index=True)        
dARImelt = pd.melt(dARI.T, var_name= 'sample', value_name='dARI', ignore_index=True)
dDivmelt['dARI'] = dARImelt['dARI']
dDivmelt['k']=np.resize(ITERATIONS[2:],dARImelt.shape[0])
dDivmelt = dDivmelt.drop('sample',axis=1)
dDivmelt = dDivmelt.apply(pd.to_numeric)

fig, ax = plt.subplots(figsize=(6,6))
g = sns.boxplot(x = 'k', y='dARI', data=dDivmelt)
ax.set_xlabel('k=', size=14)
ax.set_ylabel('d(ARI)/d(k)', size=14)
plt.savefig('Figures/Performance/Clustering_d(ARI).png', bbox_inches='tight', dpi=500)

fig, ax = plt.subplots(figsize=(6,6))
g = sns.boxplot(x = 'k', y='dDIV', data=dDivmelt)
ax.set_xlabel('k=', size=14)
ax.set_ylabel('d(Div)/d(k)', size=14)
plt.savefig('Figures/Performance/Clustering_d(Div).png', bbox_inches='tight', dpi=500)  


