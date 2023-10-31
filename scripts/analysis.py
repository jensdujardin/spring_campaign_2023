# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 14:14:00 2023

@author: maartend
"""

''' Required modules'''
import pandas as pd
import numpy as np
import seaborn as sns

from functions import assign_tax, remove_feat, transform_asinh, plot_tsne
from functions import calc_D0, calc_D1, calc_D2

from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE

''' Enter the input parameters (k=) previously used '''
x = 60           # Insert the number of nearest neighbours (k) previously used for clustering
y = 15           # Insert the desired number of nearest neighbours (k) to be used for metaclustering


''' Load the metadata and the (meta)clustering results '''
results_cell_counts_rel = pd.read_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/REL_cell_counts.csv',index_col=0)
results_median_mp = pd.read_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/Median_Cluster_SFCM.csv',index_col=0)
results_median_pop = pd.read_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/PhenoGraph_pop_median_values.csv',index_col=0, header=0)
retained_cells = pd.read_csv('Results/retained_cells.csv', index_col = 0)
metadata = pd.read_excel('metadata.xlsx', index_col = 0)

''' Extrapolate clustering results (data subset) to the entire particle count'''

# Convert relative cell counts to absolute cell concentrations using the total number of cells (after data clean up) and the measured sample volume
results_cell_counts_abs = results_cell_counts_rel.copy()

for sample in results_cell_counts_abs.index:
    results_cell_counts_abs.loc[sample,] = 1000*results_cell_counts_abs.loc[sample,]*retained_cells.loc[sample,'Total']/metadata.loc[sample,'Volume (microl)']

# Convert cell concentrations to biovolume concentrations using the average of FWS & SWS; assuming a sphere
results_cell_counts_vol = results_cell_counts_abs.copy()

for mc in results_median_mp.index:
    results_cell_counts_vol.iloc[:,mc] = results_cell_counts_abs.iloc[:,mc]  * (4/3) * np.pi * (((results_median_mp.loc[mc,'FWS Length'] + results_median_mp.loc[mc,'SWS Length'])/2)/2)**3

# Add the total cell count and total biovolume to the metadata-file
metadata['TCC (cells/mL)'] = results_cell_counts_abs.sum(axis=1)
metadata['Biovolume (μm3/mL)'] = results_cell_counts_vol.sum(axis=1)
metadata.to_excel('metadata.xlsx')


# Determine the relative contribution of each MC to the total biovolume per sample 
cell_counts_rel_vol = results_cell_counts_vol.copy()

cell_counts_vol_tot = cell_counts_rel_vol.sum(axis=1)                               # determine total biovolume per sample
cell_counts_rel_vol = cell_counts_rel_vol.div(cell_counts_vol_tot, axis=0)          # calculate relative contribution of each taxa to biovolume

# Determine the cell concentration and biovolume per functional group'''
PHYT_VAR = ['Pico-red','Pico-synecho','Nano-red','Nano-sws','Nano-crypto','Micro-red']
results_median_mp = assign_tax(results_median_mp, 'FWS Total', 6000,               # assign taxa to each MC
                    'FWS Total',17500,'SWS Maximum', 2000,'O/R ratio', 0.5, 1)
ID = list(set(results_median_mp['ID']))                                            # create list of taxa found
counts_tax = pd.DataFrame(index=results_cell_counts_abs.index)                     # create empty df to store counts per taxa
vol_tax = pd.DataFrame(index=results_cell_counts_vol.index)                        # create empty df to store biovolume per taxa

for tax in ID:                                                                     # for each taxonomic group 
    clusters = results_median_mp[results_median_mp['ID']==tax].index.tolist()      # find the metaclusters belonging to this group 
    counts_tax.loc[:,tax] = results_cell_counts_abs.iloc[:,clusters].sum(axis=1)   # calculate the sum of their cell concentrations 
    vol_tax.loc[:,tax] = results_cell_counts_vol.iloc[:,clusters].sum(axis=1)      # calculate the sum of their biovolume concentrations 

vol_tax_rel = vol_tax.div(cell_counts_vol_tot, axis=0)                             # calculate relative contribution of taxa to biovolume

# Store the results
results_cell_counts_abs.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/ABS_cell_counts.csv')
results_cell_counts_vol.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/ABS_biovolume.csv')
cell_counts_rel_vol.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/REL_biovolume.csv')
counts_tax.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/ABS_cell_counts_tax.csv') 
vol_tax.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/ABS_biovolume_tax.csv')
vol_tax_rel.to_csv('Results/k=' + str(x) + '/k2=' + str(y) + '/REL_biovolume_tax.csv')

''' Perform dimensionality reduction using t-distributed Stochastic Neighbor Embedding (t-SNE) to visualise the metaclustering results'''
FEATURES = remove_feat(results_median_mp, ['Sample Length', 'Arrival Time','Sample','Proportion','Cluster','MC','ID']) 

# Transform median cluster values using f(x) = asinh(x) and standardize median values for t-sne to mimick the clustering process
scaler = StandardScaler()
results_median_pop_asinh = transform_asinh(results_median_pop, FEATURES)
results_median_pop_stand = scaler.fit_transform(results_median_pop_asinh.loc[:,FEATURES])

# Calculate the t-SNE components (2 dimensions)
tsne = TSNE(n_components=2, init='pca')
tsne_pop = tsne.fit_transform(results_median_pop_stand)
tsne_pop = pd.DataFrame(tsne_pop, index=results_median_pop.index, columns=['t-SNE 1', 't-SNE 2'])
results = pd.concat([tsne_pop,results_median_pop], axis=1)
results['Day'] = results['Sample'].str[0:16]

# Plot t-SNE for individual cluster and metacluster functional identity '''
palette = sns.color_palette('Spectral_r', as_cmap=True)
plot_tsne(results, 't-SNE_ID_pop', hue_var='Cluster', size_var='Proportion', legend =False)             # plot all clusters by cluster-ID
plot_tsne(results, 't-SNE_ID_mc', hue_var='MC', size_var='Proportion', legend=True, pal=palette)        # plot all clusters by metaclusters
plot_tsne(results, 't-SNE_Taxa', hue_var='ID', hue_order=PHYT_VAR, size_var='Proportion', legend=True)  # plot all clusters by taxa

''' Calculate and plot the cytometric alpha-diversity indices'''
d0 = calc_D0(results_cell_counts_rel)                          # Calculate Hill-index 0
d1 = calc_D1(results_cell_counts_rel)                          # Calculate Hill-index 1
d2 = calc_D2(results_cell_counts_rel)                          # Calculate Hill-index 2
div = pd.concat([d0,d1,d2], ignore_index=False, axis=1)        # Merge dataframes
div.to_csv('Results/HillDiversity.csv')

''' Determine Bray-Curtis dissimilarities between relative cell counts and PCA with MDS to identify possible communities'''
from sklearn.neighbors import DistanceMetric
from sklearn.manifold import MDS
import matplotlib.pyplot as plt

bc = DistanceMetric.get_metric('braycurtis')
cont_table_bc = bc.pairwise(results_cell_counts_rel)
mds = MDS(dissimilarity='precomputed', n_init=100, metric=True)
cont_table_mds = pd.DataFrame(mds.fit_transform(cont_table_bc), columns=['PCoA 1','PCoA 2'], index=results_cell_counts_rel.index)
cont_table_mds['Year-Month'] = cont_table_mds.index.str[0:7]

fig, ax = plt.subplots() 
g = sns.scatterplot(x='PCoA 1', y='PCoA 2', hue='Year-Month', data=cont_table_mds, legend='brief', palette = 'Spectral_r', alpha=1, ax=ax)
sns.despine(ax=ax)
ax.set_xlabel("PCoA 1", fontsize=18)
ax.set_ylabel("PCoA 2", fontsize=18)    
plt.setp(g.get_xticklabels(), size=14)
plt.setp(g.get_yticklabels(), size=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('Figures/PCoA_MC_trip_REL.png',bbox_inches='tight', dpi=500)
plt.show()