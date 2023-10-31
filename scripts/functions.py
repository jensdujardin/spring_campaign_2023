# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 09:39:33 2022

@authors: maartend & prubbens
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

''' Remove cells that are larger or smaller than threshold for a specific feature

Parameters: 
----------
df: pandas dataframe
feature: feature used to threshold
threshold: threshold value below or above which cells will be retained
maximum: boolean indicating whether the threshold is a maximum or minimum

returns subset of filtered cells
'''

def remove_cells(df, feature, threshold, maximum=True):        # Remove cells that are larger or smaller than threshold for a specific feature
    if maximum: df = df.loc[(df.loc[:,feature] < threshold)]   # Input: pandas dataframe (df); threshold feature (feature); threshold value (threshold)
    else: df = df.loc[(df.loc[:,feature] > threshold)]         # Output: return subsetted df
    return df   


''' Remove variables from dataframe that are in list features_to_remove

Parameters: 
----------
df: pandas dataframe
features_to_remove: list of features to remove

returns dataframe without removed features
'''
def remove_feat(df, features_to_remove):                        # Create function to select and store column names
    for feature in features_to_remove: 
        df = df.loc[:, ~df.columns.str.endswith(feature)]
    return df.columns.to_list()


''' Transform each variable in the dataframe df by f(x) = asinh(x) 

Parameters:
----------
df: dataframe
features: variables to transform
'''
def transform_asinh(df, features):                              # Transform each variable in the dataframe df by log(x)
    return df.loc[:,features].apply(np.arcsinh)


''' Assign taxonomic group to cluster

df: dataframe
fws_thresh_pico: threshold to distinguish pico-phytoplankton based on FWS Length
fws_thresh_nano: threshold to distinguish nano-phytoplankton based on FWS Length
fws_thresh_micro: threshold to distinguish micro-phytoplankton based on FWS Length

'''
def assign_tax(df, fws_var_pico, fws_thresh_pico, fws_var_micro, fws_thresh_micro,
               sws_var_cocco, sws_thresh_cocco, fl_var_synecho, fl_thresh_synecho, fl_thresh_crypto): 
    df.loc[:,'O/R ratio'] = df.loc[:,'FL Orange Total'] / df.loc[:,'FL Red Total']
    df.loc[:,'O/Y ratio'] = df.loc[:,'FL Orange Total'] / df.loc[:,'FL Yellow Total']
    df.loc[:,'ID'] = 'ND'

    for cluster in df.index:
        if df.loc[cluster,fws_var_pico] < fws_thresh_pico: 
            if (df.loc[cluster,fl_var_synecho] > fl_thresh_synecho):
                df.loc[cluster,'ID'] = 'Pico-synecho'
            else: 
                df.loc[cluster,'ID'] = 'Pico-red'
        if df.loc[cluster,fws_var_micro] >= fws_thresh_micro:
            df.loc[cluster,'ID'] = 'Micro-red'
        if df.loc[cluster,'ID'] == 'ND': 
            if df.loc[cluster,sws_var_cocco] > sws_thresh_cocco: 
                df.loc[cluster,'ID'] = 'Nano-sws' 
            elif df.loc[cluster,'O/R ratio'] < 1: 
                df.loc[cluster,'ID'] = 'Nano-red'
            else: 
                df.loc[cluster,'ID'] = 'Nano-crypto'
            
    return df


''' Determine median values for each metacluster '''
def calc_median_mc(df, var, var_mc, path, list_df_mc): 
    df_median = pd.DataFrame(index=df.columns, columns=var)
    for mc in df.columns: 
        df_mc = pd.DataFrame()
        for file in list_df_mc: 
            print('Metacluster',mc,'in',file)
            df = pd.read_csv(path + file, index_col=0, header=0)
            df = df[df.loc[:,var_mc] == mc]
            df_mc = pd.concat([df_mc,df], axis=0)
            df_mc = df_mc.select_dtypes('number')
        df_median.loc[mc,var] = df_mc.median()
    
    return df_median

''' Make visualization of t-SNE dimensionality reduction, laying over a third variable 
to visualize phenotypic (or other) properties

Parameters: 
-----------

df: dataframe
str_name: string identifier to  provide a name for the figure
pal: color palette used to visualize hue_var
hue_var: hue variable that is laid over the t-SNE map
size_var: variable to visualize size (i.e. proportion) of each data point
legend: boolean that indicates whether to add a legend to the figure or not

'''
def plot_tsne(df, str_name, pal='colorblind', hue_var=None, hue_order=None,
              size_var=None, legend=True):
    if legend ==  True:
        fig, ax = plt.subplots() 
        g = sns.scatterplot(x='t-SNE 1', y='t-SNE 2', data=df, palette=pal,
                            hue=hue_var, hue_order=hue_order, size=size_var, legend='auto', alpha=0.5)
        sns.despine(ax=ax)
        ax.set_xlabel("t-SNE 1", fontsize=18)
        ax.set_ylabel("t-SNE 2", fontsize=18)    
        plt.setp(g.get_xticklabels(), size=14)
        plt.setp(g.get_yticklabels(), size=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig('Figures/' + str_name + '.png', bbox_inches='tight', dpi=500)
        plt.show()
    else: 
        fig, ax = plt.subplots() 
        g = sns.scatterplot(x='t-SNE 1', y='t-SNE 2', data=df, palette=pal,
                            hue=hue_var, hue_order=hue_order, size=size_var, legend=False, alpha=0.5)
        sns.despine(ax=ax)
        ax.set_xlabel("t-SNE 1", fontsize=18)
        ax.set_ylabel("t-SNE 2", fontsize=18)    
        plt.setp(g.get_xticklabels(), size=14)
        plt.setp(g.get_yticklabels(), size=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig('Figures/' + str_name + '.png', bbox_inches='tight', dpi=500)
        plt.show()


''' Calculate Hill diversity for q=0 (richness), q=1 (exponential Shannon index) & q=2 (reciprocal Simpson index)

Parameters:
----------
df: dataframe containing abundance data

'''
def calc_D0(df):
    result = df.astype(bool).sum(axis=1)
    return pd.DataFrame(result, columns=[r'D0'])

def calc_D1(df): 
    return pd.DataFrame(np.exp(stats.entropy(df.T.values)), index=df.index, columns=[r'D1'])
 
def calc_D2(df): 
    df_sum = df.sum(axis=1)
    df_final = (df.div(df_sum, axis=0))**2
    return pd.DataFrame(1./df_final.sum(axis=1), columns=[r'D2'])

''' Plot data geographically using longitude and latitude coordinates

Parameters: 
-----------

df: dataframe
back_img: background image
bbox: bounding box containing minimum and maximum longitude and latitude 
[lang_min, lang_max, long_min, long_max]
hue_var: hue variable that is laid over spatiotemporal map
str_id: string identifier to save figure
str_pal: string to choose color palette

'''
def plot_cruise_data(df, back_img, bbox, hue_var, str_id, str_pal): 
    fig, ax = plt.subplots() 
    ax.imshow(back_img, zorder=0, extent = bbox, aspect= 'auto', alpha=0.8)
    g = sns.scatterplot(x='Longitude', y='Latitude', hue=hue_var, data=df, legend='brief', palette = str_pal, alpha=1, ax=ax)
    ax.set_xlabel("Longitude", fontsize=18)
    ax.set_ylabel("Latitude", fontsize=18)    
    ax.set_xlim(bbox[0],bbox[1])
    ax.set_ylim(bbox[2],bbox[3])
    plt.setp(g.get_xticklabels(), size=14)
    plt.setp(g.get_yticklabels(), size=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title=str_id,
               fontsize=12, title_fontsize=14)
    plt.savefig('Figures/' + hue_var + '_lat_long.png',bbox_inches='tight', dpi=500)
    plt.show()
    plt.close()