import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import gpmap.src.plot.mpl as mplot

from os.path import join
from scripts.settings import DATA_DIR, FIGURES_DIR
from scripts.utils import (plot_seq_names, format_pos, get_subsequences,
                           get_edges_df, rm_alleles, calc_lims, format_subseq)


def plot_landscape(axes, ndf, path, int_label):
    ref = path[0]
    perc['d'] = [np.sum([x1 != x2 for x1, x2 in zip(ref, s2)])
                         for s2 in perc.index.values]
    xlim = calc_lims(perc['d'], p=0.2)
    ylim = calc_lims(perc['fit'], p=0.2)
    
    edf = get_edges_df(ndf.index.values)
    mplot.plot_visualization(axes, ndf, edges_df=edf, 
                             x='d', y='fit', nodes_color='black',
                             nodes_size=30,
                             edges_color='grey', edges_width=1,
                             edges_alpha=0.5)
    
    edf = get_edges_df(path)
    mplot.plot_visualization(axes, ndf.loc[path], edges_df=edf, 
                             x='d', y='fit', nodes_color='black',
                             nodes_size=30,
                             edges_color='black', edges_width=1.1,
                             edges_alpha=1)
    
    plot_seq_names(ndf, axes, x='d', y='fit')
    axes.set(xlabel='Hamming distance',
             ylabel='% viable {} genotypes'.format(int_label),
             xlim=xlim, ylim=ylim)



if __name__ == '__main__':
    label = 'model2'
    fpath = join(DATA_DIR, 'calibrated_predictions.pq')
    thresholds = [1, 1.25, 1.5, 1.75, 2, 2.25]
    nthresholds = len(thresholds)
    configs = [
               (['E:98', 'Im:33', 'Im:34'], ['E:83'], 'K', np.array(['R/DN', 'R/DV', 'M/DV', 'M/LV']), {0: 'V', 2: 'V', 3: 'D'}),
               (['E:97', 'Im:38'], [], '', np.array(['E/R', 'K/R', 'K/T']), {}),
               (['E:72', 'E:73'], ['E:77', 'E:83'], 'TK', np.array(['KA', 'NA', 'NP']), {}),
               (['E:72', 'E:73'], ['E:77', 'E:83'], 'SK', np.array(['KA', 'NA', 'NP']), {}),
               (['E:72', 'E:73'], ['E:77', 'E:83'], 'TY', np.array(['KA', 'NA', 'NP']), {}),
               (['E:72', 'E:73'], ['E:77', 'E:83'], 'SY', np.array(['KA', 'NA', 'NP']), {})
               ]
    npanels = len(configs)
    
    fig, subplots =  plt.subplots(npanels, nthresholds,
                                  figsize=(nthresholds * 3, npanels * 3))

    print('Reading data from {}'.format(fpath))
    data = pd.read_parquet(fpath)
    data.loc[data['y'] > 6, 'y'] = 6
    
    for i, config in enumerate(configs):
        int_pos, context_pos, seq, path, alleles_remove = config
        int_label = format_pos(int_pos)
        print('Interaction {}'.format(int_label))
        if int_label not in data.columns:
            data[int_label] = get_subsequences(int_pos, data.index.values) 
        
        if context_pos:
            context_label = format_pos(context_pos)    
            if context_label not in data.columns:
                data[context_label] = get_subsequences(context_pos, data.index.values)
            df = data.loc[data[context_label] == seq, ['y', int_label]].copy()
        else:
            df = data[['y', int_label]].copy()
        
        for j, x in enumerate(thresholds):
            print('\tWith Norm MIC score > {}'.format(x))
            axes = subplots[i, j]
            df['fit'] = (df['y'] > x).astype(float)
            perc = df.groupby([int_label])[['fit']].mean() * 100
            
            if alleles_remove:
                perc = rm_alleles(perc, alleles_remove)
                
            plot_landscape(axes, perc, path, int_label)
            
            if j > 0:
                axes.set_ylabel('')
            if i < (npanels - 1):
                axes.set_xlabel('')
            if i == 0:
                axes.set_title('Norm. MIC score > {}'.format(x))
            axes.set_xticks(np.arange(len(int_pos)+1))
            
            if context_pos:
                label = format_subseq(context_pos, seq)
                axes.text(0.97, 0.07, label, transform=axes.transAxes,
                          fontsize=10, ha='right', va='top')

    fig.tight_layout()
    fig.savefig(join(FIGURES_DIR, 'figS14.svg'), dpi=300, format='svg')
    fig.savefig(join(FIGURES_DIR, 'figS14.png'), dpi=300, format='png')
    