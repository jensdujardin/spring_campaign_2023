# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 09:39:33 2022

@authors: maartend & prubbens
"""

'''
This code should be placed in a folder containing the associated functions-script, 
a metadata-file (metadata.xlsx), and clustered data in a subfolder called 'CSV_clustered' 
'''

''' Required packages (install using cmd / powershell)'''  
# conda install seaborn scikit-learn statsmodels numba pytables
# conda install -c conda-forge python-igraph leidenalg
# pip install scanpy
# pip install -U PhenoGraph

# This script needs to placed in a directory containing clustered data in a subfolder called CSV_clustered and a metadatafile in the wd

#%% 'Identify common clusters (metaclusters) between samples'

x = 40                # Insert the number of nearest neighbours (k1) previously used for single-file clustering
N_cells = 3000       # Insert the number of cells used for subsampling during clean up

''' Setting and checking working directory'''
import os
import sys
os.chdir(os.path.abspath(os.path.dirname(sys.argv[0])))    # change the working directory to the directory containing the script file
print(os.getcwd())                                         # check the working directory 

''' Required modules'''
import datetime
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from os import listdir
from scanpy.external.tl import phenograph
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.manifold import TSNE

'''Import custom functions'''
from functions import remove_feat               # Removes unwanted columns using an input column list
from functions import transform_asinh           # Log-transforms the raw data
from functions import assign_tax                # Classifies metaclusters in large functional groups
from functions import calc_median_mc            # Calculates median values of all the input variables for each metacluster

# locate results of clustering
PATH='CSV_cleaned/'
SAMPLES = sorted(listdir(PATH))
PATH = 'CSV_clustered/k=' + str(x) + '/'
FILENAMES = sorted(listdir(PATH))

''' Read in the metadata to obtain sampling volumes'''
metadata = pd.read_excel('metadata.xlsx', index_col = 0)

'''Metacluster each clustered datafile using Phenograph, repeat for a list of k2-values'''
ITERATIONS = [7,9,11,13,15,17,19,21,23,25]                                 # Provide a list of k-values to attempt in metaclustering    

for y in ITERATIONS:  

    if not os.path.exists('CSV_metaclustered/k1=' + str(x) + '-k2='+str(y)):          # Create new folder to store netaclustered FCM data if needed
        os.makedirs('CSV_metaclustered/k1=' + str(x) + '-k2='+str(y)) 

    if not os.path.exists('Results/k=' + str(x) + '/k2='+str(y)):                     # Create new folder to store netaclustering summaries
        os.makedirs('Results/k=' + str(x) + '/k2='+str(y)) 

    PATHMC = 'CSV_metaclustered/k1=' + str(x) + '-k2='+str(y) + '/'                   # Set output location for the metaclustered data-files
    
    ''' Read in clustering results'''
    df_cell_counts = pd.read_csv('Results/k=' + str(x) + '/PhenoGraph_pop_cell_counts.csv',index_col=0, header=0)               # Contains the number of cells per cluster
    df_median_cluster_val = pd.read_csv('Results/k=' + str(x) + '/PhenoGraph_pop_median_values.csv',index_col=None, header=0)   # Contains the median measured values per parameter of each cluster

    ''' Create list of features ''' 
    FEATURES = remove_feat(df_median_cluster_val, ['Sample Length', 'Arrival Time','Cluster','Sample']) # list all measured features without sample length and arrival time

    ''' Create a df with counts for each cluster present (no NAs) '''
    df_cell_counts.columns = df_cell_counts.columns.astype(int)                                         # copy clustering results in a new df
    df_cell_counts_melt = pd.melt(df_cell_counts.T, var_name='Sample', value_name='Proportion')         # unpivot the new cell counts dataframe
    df_cell_counts_melt.dropna(inplace=True)                                                            # dropping NAs creates a df with only counts for each cluster present
 
    ''' Create a df with relative abundance for each cluster present (no NAs) '''                                              
    df_cell_counts_rel = df_cell_counts.div(df_cell_counts.sum(axis=1), axis=0)                         # divide each value in the cell count df by the sum of that row
    df_cell_counts_rel_melt = pd.melt(df_cell_counts_rel.T, var_name='Sample', value_name='Proportion') # unpivot the new df
    df_cell_counts_rel_melt.dropna(inplace=True)                                                        # dropping NAs creates a df with relative abundance per cluster present

    ''' Create a new dataframe to store results of metaclustering '''
    sfcm_meta = pd.DataFrame()

    ''' Log-transform and normalize the median cluster values '''
    scaler = StandardScaler()
    df_median_cluster_val_asinh = transform_asinh(df_median_cluster_val, FEATURES)
    df_median_cluster_val_norm = pd.DataFrame(scaler.fit_transform(df_median_cluster_val_asinh.loc[:,FEATURES]), index=df_median_cluster_val_asinh.index, columns=FEATURES)

    ''' Group filtered clusters per sample into metaclusters (MCs) '''
    metaclusters, graph, Q = phenograph(df_median_cluster_val_norm.loc[:,FEATURES], k=y, primary_metric='Euclidean', min_cluster_size=1)
    df_median_cluster_val['MC'] = metaclusters

    ''' Store absolute and relative counts of the metaclusters '''
    df_median_cluster_val['Proportion'] = df_cell_counts_rel_melt.loc[:,'Proportion'].values

    ''' Assign metaclusters to functional groups based on traditional thresholds (FWS, SWS and O/R ratio)'''
    df_median_cluster_val = assign_tax(df_median_cluster_val, 'FWS Total', 6000,'FWS Total', 17500,'SWS Maximum', 2000,'O/R ratio', 0.5, 1)

    ''' Store number of populations/clusters (not metaclusters) per sample '''
    n_pop = df_cell_counts.count(axis=1)                    # isolate the sample names in a single column df
    n_pop = pd.DataFrame(n_pop)                             # add the number of clusters per sample in the second column
    n_pop.columns = ['n_pop']                               # name the new column

    ''' Create a global dataframe to store cell counts per metacluster per sample '''
    df_cell_counts_mc = pd.DataFrame(np.zeros((len(df_cell_counts.index),int(2+np.max(metaclusters)))),index= df_cell_counts.index,columns=np.arange(-1,np.max(metaclusters)+1))

    ''' Create a global dataframe to store cell counts per functional population per sample '''
    PHYT_VAR = ['Pico-red','Pico-synecho','Nano-red','Nano-sws','Nano-crypto','Micro-red']
    df_cell_counts_mc_tax = pd.DataFrame(np.zeros((len(df_cell_counts.index),len(PHYT_VAR))),index=df_cell_counts.index, columns=PHYT_VAR)

    ''' Create a global dataframe to store biovolume per sample '''
    df_biovol_mc = pd.DataFrame(np.zeros((len(df_cell_counts.index),int(2+np.max(metaclusters)))),index= df_cell_counts.index,columns=np.arange(-1,np.max(metaclusters)+1))

    ''' Loop over the clustered datafiles to assign cells to metaclusters'''
    for index in df_cell_counts.index:
        df_median_cluster_val_subset = df_median_cluster_val.loc[df_median_cluster_val.loc[:,'Sample'] == index]        # open df with median values of the metaclusters 
        df = pd.read_csv(PATH + index[0:len(index)-4] + '_PhenoGraph_k=' + str(x) + '.csv', index_col=0, header=0)      # open 1 file
        df.loc[:,'MC'] = -1                                                                                             # create new column called MC with -1 values
    
        ''' Loop over all counts per cluster per sample to aggregate results into cell 
        counts per metacluster '''
        for pop_id in df_median_cluster_val_subset.index:                      # for each population in the df with metaclusters
            print(datetime.datetime.now())
            print('Assigning cells in sample:',index)   
            print('Processing cluster number:', pop_id)    
            pop = df_median_cluster_val.loc[pop_id,'Cluster']                  # determine its corresponding single-file cluster-number in this sample
            counts = df_cell_counts.loc[index[0:len(index)],pop]               # determine the number of cells of this metacluster in this sample
            mc = df_median_cluster_val_subset.loc[pop_id,'MC']                 # determine its metacluster number
            df_cell_counts_mc.loc[index,mc] += counts                          # store the samplename, metacluster number and cell count in the first global df
            idx_pop = df[df.loc[:,'Cluster'] == pop].index                     # find row-number of all cells of this cluster in the sample
            df.loc[idx_pop,'MC'] = mc                                          # add metacluster number to the sample datafile      
            df.to_csv('CSV_metaclustered/k1=' + str(x) + '-k2='+str(y)+ '/' +  # Store a new copy of the clustered datafile with the metaclusters added
                      index[0:len(index)-4] + '_PhenoGraph_k=' + str(x) + '.csv')        
        
            id_tax = df_median_cluster_val.loc[pop_id,'ID']                    # determine the functional group of this metacluster
            df_cell_counts_mc_tax.loc[index,id_tax] += counts                  # store the samplename, functional group and cell count in the second global df
            
        ''' Export the results of the metaclustering to the global dataframe '''
        sfcm_meta.loc[index,'Number of clusters'] = n_pop.loc[index,'n_pop']
        sfcm_meta.loc[index,'Number of metaclusters'] = len(df['MC'].unique())
        sfcm_meta.loc[index,'Fraction of cells clustered'] = df_cell_counts_mc.loc[index,:].sum() / N_cells
        
    ''' Drop unassigned clusters (denoted by -1) '''
    df_cell_counts_mc.drop(-1, axis=1, inplace=True)

    ''' Determine the relative cell counts using total cell concentration (hence corrected for) measured volume '''
    df_cell_counts_tot = df_cell_counts_mc.sum(axis=1)
    df_cell_counts_rel = df_cell_counts_mc.div(df_cell_counts_tot, axis=0)

    ''' Determine median detector values for each metacluster '''
    df_median_mc = calc_median_mc(df_cell_counts_rel, FEATURES, 'MC', PATHMC, FILENAMES)

    ''' Determine cell counts per functional cell population '''
    df_cell_counts_tot_tax = df_cell_counts_mc_tax.sum(axis=1)
    df_cell_counts_rel_tax = df_cell_counts_mc_tax.div(df_cell_counts_tot_tax, axis=0)
    
    ''' Write out results '''
    df_cell_counts_rel.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/REL_cell_counts.csv')                     # Contains the fraction of cells per metacluster in each sample
    df_cell_counts_rel_tax.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/REL_cell_counts_tax.csv')             # Contains the fraction of cells per taxonomic group in each sample
    df_median_cluster_val.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/PhenoGraph_pop_median_values.csv')     # Contains the median of all parameters for each metacluster
    df_median_mc.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/Median_Cluster_SFCM.csv')                       # Contains the global median of all parameters across all particles
    sfcm_meta.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/SFCM_overview.csv')                                # Contains a summary of the clustering and metaclustering steps
    print('\007')                                                                                                   # (on Windows) play an audio alert

'''Generate summaries of the diversity and Adjusted Rand Scores for each k-value attempted'''
ITERATIONS = [7,9,11,13,15,17,19,21,23,25]            

# Calculate and plot the species richness
df_cell_counts = pd.read_csv('Results/k=' + str(x) + '/PhenoGraph_pop_cell_counts.csv',index_col=0, header=0)       # open summary of clustering to obtain the number of clusters per sample
n_pop = df_cell_counts.count(axis=1)                                                                                # isolate the sample names in a single column df
n_pop = pd.DataFrame(n_pop)                                                                                         # add the number of clusters per sample in the second column
n_pop.columns = ['n_pop']                                                                                           # name the new column

for z in ITERATIONS:                                                                                                # for each k-value
    temp = pd.read_csv('Results/k=' + str(x) + '/k2=' + str(z) + '/REL_cell_counts.csv',index_col=0, header=0)      # open the relative cell count file
    n_pop['k=' + str(z)] = temp.select_dtypes(np.number).gt(0).sum(axis=1)                                          # determine the number of metaclusters

output = plt.figure()                                                                      # create a boxplot with the number of metaclusters for each k2-value
n_pop.boxplot(return_type='axes', rot=45, grid=False) 
plt.ylabel('# (meta)clusters', size=16)
output.savefig('Figures/Performance/Metaclustering_div.png', bbox_inches='tight', dpi=500)                 # Store the figure
n_pop.to_csv('Results/k=' + str(x) + '/Metaclustering.csv')                                # Store the diversity overview in a csv-file                                                

# Calculate and plot the Adjusted rand scores
rand = pd.DataFrame(columns=[ITERATIONS],index=FILENAMES)                                                           # create an empty dataframe to store the ARI

for i in range(1,len(ITERATIONS)):                                                                                  # for each k value minus the lowest
    print('Calculating ARI for k2=', ITERATIONS[i])     
    for z in range(0,len(FILENAMES)):                                                                                # for each file
       file = FILENAMES[z]                                                                                          # determine the filename
       y1 = ITERATIONS[i-1]                                                                                         # determine two adjacent k-values; store as y1 and y2
       y2 = ITERATIONS[i]                                                                                           # open the metaclustered data for these two k-values
       df1 = pd.read_csv('CSV_metaclustered/k1=' + str(x) + '-k2='+str(y1)+ '/' + file[0:(len(file)-4)] + '.csv', header=0) 
       df2 = pd.read_csv('CSV_metaclustered/k1=' + str(x) + '-k2='+str(y2)+ '/' + file[0:(len(file)-4)] + '.csv', header=0)  
       df1 = df1[['Particle ID','MC']]                                                                              # store the pairings of particle ID and metacluster ID for y1
       df2 = df2[['Particle ID','MC']]                                                                              # store the pairings of particle ID and metacluster ID for y2
       df = pd.merge(df1, df2, how='inner', on=['Particle ID'],copy=False)                                          # merge the pairings 
       rand.loc[FILENAMES[z],ITERATIONS[i]] = adjusted_rand_score(df['MC_x'],df['MC_y'])                            # Calculate the Adjusted Rand Score (part of scikit-learn)
rand.to_csv('Results/MetaRandScore.csv')                                                                            # Store the results in a csv-file

MetaRand = pd.read_csv('Results/MetaRandScore.csv',index_col=0, header=0)                  # Load the results without headers and multi-indexing
output = plt.figure()                                                                      # create a boxplot with the ARI for each k2-value
MetaRand.boxplot(return_type='axes', rot=45, grid=False) 
plt.xlabel('k =', size=16)
plt.ylabel('ARI', size=16)
output.savefig('Figures/Performance/Metaclustering_ARI.png', bbox_inches='tight', dpi=500)                 # Store the figure

# Calculate and plot the derivatives of species richness and the Adjusted rand scores
n_pop = pd.read_csv('Results/k=' + str(x) + '/Metaclustering.csv',index_col=0, header=0)                 
dARI = pd.DataFrame(columns=MetaRand.columns[2:len(MetaRand.columns)],index=MetaRand.index)    
dDiv = pd.DataFrame(columns=n_pop.columns[2:len(MetaRand.columns)],index=n_pop.index) 
dDiv.drop(columns=dDiv.columns[0], axis=1, inplace=True)

for i in range(2,len(ITERATIONS)):
    for z in range(0,len(dARI.index)):
        dARI.loc[dARI.index[z],MetaRand.columns[i]] = (MetaRand.loc[dARI.index[z],MetaRand.columns[i]] - MetaRand.loc[dARI.index[z],MetaRand.columns[i-1]]) / (int(MetaRand.columns[i])-int(MetaRand.columns[i-1]))
        dDiv.loc[dDiv.index[z],n_pop.columns[i]] = (n_pop.loc[dDiv.index[z],n_pop.columns[i]] - n_pop.loc[dDiv.index[z],n_pop.columns[i-1]]) / (int(n_pop.columns[i][2:len(n_pop.columns[i])])-int(n_pop.columns[i-1][2:len(n_pop.columns[i])]))

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
plt.savefig('Figures/Performance/Metaclustering_d(ARI).png', bbox_inches='tight', dpi=500)

fig, ax = plt.subplots(figsize=(6,6))
g = sns.boxplot(x = 'k', y='dDIV', data=dDivmelt)
ax.set_xlabel('k=', size=14)
ax.set_ylabel('d(Div)/d(k)', size=14)
plt.savefig('Figures/Performance/Metaclustering_d(Div).png', bbox_inches='tight', dpi=500)                                                                           # Store the results in a csv-file

# Visualise the metaclustering results using t-distributed Stochastic Neighbor Embedding (t-SNE)
if not os.path.exists('Figures/Performance/tSNE/'): os.makedirs('Figures/Performance/tSNE/')                           # Create a new folder if needed

tsne = TSNE(n_components=2, init='pca')                                                    # input parameters for the t-SNE function 
scaler = StandardScaler()                                                                  # picking the same scaler as used before 

for i in range(0,len(ITERATIONS)):                                                         # for each k2-value, plot the metaclusters on a t-SNE using their median parameter values
    results_median_pop = pd.read_csv('Results/k=' + str(x) + '/k2=' + str(ITERATIONS[i]) + '/PhenoGraph_pop_median_values.csv',index_col=0, header=0)
    results_median_pop_asinh = transform_asinh(results_median_pop, FEATURES)
    results_median_pop_stand = scaler.fit_transform(results_median_pop_asinh.loc[:,FEATURES])
    tsne_pop = tsne.fit_transform(results_median_pop_stand)
    tsne_pop = pd.DataFrame(tsne_pop, index=results_median_pop.index, columns=['t-SNE 1', 't-SNE 2'])
    results = pd.concat([tsne_pop,results_median_pop], axis=1)

    fig, ax = plt.subplots() 
    g = sns.scatterplot(x='t-SNE 1', y='t-SNE 2', data=results, palette='Spectral',
                            hue='MC', size='Proportion', legend=False, alpha=0.5)
    sns.despine(ax=ax)
    ax.set_xlabel("t-SNE 1", fontsize=18)
    ax.set_ylabel("t-SNE 2", fontsize=18)
    ax.set_title("k2= "+ str(ITERATIONS[i]), fontsize=18)       
    plt.setp(g.get_xticklabels(), size=14)
    plt.setp(g.get_yticklabels(), size=14)
    plt.savefig('Figures/Performance/tSNE/tSNE_MC_k2=' + str(ITERATIONS[i]) + '.png', bbox_inches='tight', dpi=500)
    plt.show()

