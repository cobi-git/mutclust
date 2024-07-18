from multiprocessing import Pool
import os
from os import  listdir,  walk
from os.path import isfile, join, isdir
import numpy as np
from scipy.stats import entropy
import pandas as pd
from pandas import read_csv, Series
import pickle

import matplotlib.pyplot as plt
import matplotlib.transforms
import matplotlib.patches

HSCORE = 'H-score'
POS = 'Position'
FREQ = 'Frequency'
PER = 'Percentage'
ENT = 'Entropy'
NUCLEOTIDE_ANNOTATION_PATH = '/data3/projects/2020_MUTCLUST/Data/Annotation/Nucleotide/covid_annotation.tsv'
CLUSTER_INFO_PATH = '/data3/projects/2020_MUTCLUST/Data/Output/GISAID/Mutclust/clusters_avr.txt'

def readPickle(filepath):
    with open(filepath, 'rb') as file:
        return pickle.load(file)
    
def mutation_filtering(mutclust_input_df):
    new_mutclust_input_df = mutclust_input_df.loc[ (mutclust_input_df[POS] >= 266) & (mutclust_input_df[POS] <= 29674)]
    new_mutclust_input_df.reset_index(drop=True, inplace=True)
    return new_mutclust_input_df

def get_GeneInfo_df():
    annotation_df = read_csv(NUCLEOTIDE_ANNOTATION_PATH, sep=' ')
    return annotation_df
    
def make_bedgraph(mutClustInput_df, cluster_df):
    color_list = ['maroon', 'r', 'coral', 'chocolate', 'orange', 'gold', 'olive', 'yellow', 'lawngreen', 'palegreen', 'forestgreen', 'lime', 'mediumaquamarine', 'aquamarine', 'teal', 'aqua', 'steelblue', 'slategrey', 'cornflowerblue', 'blue', 'slateblue', 'indigo', 'plum', 'magenta', 'deeppink', 'pink']

    col_list = ['H-score']
    gene_df = get_GeneInfo_df()
    #cluster_df = pd.read_csv(CLUSTER_INFO_PATH, sep='\t')
    input_df = mutClustInput_df[col_list]

    plt.figure(constrained_layout=True, figsize=(12, 8))
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(13.5, 3))
    axes = [ax] if len(col_list) == 1 else ax

    for ax, col in zip(axes, input_df):
        plot_df = input_df[col]
        plot_df.plot.bar(width=30, ax=ax, xticks=np.arange(0, 29903, 1000), ylabel=col)
        ax.tick_params(labelrotation=45, labelsize=5)
    
        for i, row in cluster_df.iterrows():
            ax.axvspan(row['left_position'], row['right_position'], facecolor='lightgrey', alpha=0.8)

    pre_y = 0
    for i, row in gene_df.iterrows():
        trans = matplotlib.transforms.blended_transform_factory(axes[-1].transData, fig.transFigure)
        r = matplotlib.patches.Rectangle((row['start'], 0.02), row['end'] - row['start'] + 1, 0.02, facecolor=color_list[i], transform=trans, edgecolor='black', lw=0.5)
        fig.add_artist(r)
        text_y = 0.01
        if ((row['start'] - row['end']) < 1000) and (pre_y == 0.01):
            text_y = 0.05
        pre_y = text_y
        fig.text((row['start'] + row['end']) / 2, text_y, row['gene'], ha='center', va='center', fontsize=5, fontweight='bold', zorder=10, transform=trans)
    
    #fig.suptitle('Clustering result')
    plt.tight_layout()
    plt.show()
    
if __name__ == '__main__':
    print()