import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import gpmap.src.plot.ds as dplot
import gpmap.src.plot.mpl as mplot

from os.path import join
from scripts.utils import load_visualization, get_subsequences
from scripts.settings import FIGURES_DIR


if __name__ == '__main__':
    x, y = '1', '2'
    cmap = cm.get_cmap('coolwarm')
    palette = {'Other': cmap(0.1), 
               'E:Lys72': cmap(0.67),
               'E:Tyr83': cmap(0.9),
               'E:Lys72-Tyr83': cmap(0.82)}
    
    print('Loading visualization data')
    nodes, edges = load_visualization()
    
    print('Extracting subsequences of interest')
    nodes['E:72'] = get_subsequences(['E:72'], nodes.index.values)
    nodes['E:83'] = get_subsequences(['E:83'], nodes.index.values)
    nodes['v'] = 0
    nodes.loc[nodes['E:72'] == 'K', 'v'] = 0.95
    nodes.loc[nodes['E:83'] == 'Y', 'v'] = 0.9
    nodes.loc[(nodes['E:83'] == 'Y') & (nodes['E:72'] == 'K'), 'v'] = 1
    
    nodes['c'] = 'Other'
    nodes.loc[nodes['E:72'] == 'K', 'c'] = 'E:Lys72'
    nodes.loc[nodes['E:83'] == 'Y', 'c'] = 'E:Tyr83'
    nodes.loc[(nodes['E:83'] == 'Y') & (nodes['E:72'] == 'K'), 'c'] = 'E:Lys72-Tyr83'
    
    print('Rendering edges for both panels')
    edges_dsg = dplot.plot_edges(nodes, edges, x, y,
                                 shade=True, resolution=1200)
    
    print('Rendering nodes for normalized MIC score')
    dsg1 = dplot.plot_nodes(nodes, x=x, y=y, color='function',
                            cmap='coolwarm', resolution=800,
                            sort_by='function', sort_ascending=False)
    dsg = (edges_dsg * dsg1).opts(padding=0.05, aspect='square', title='',
                                  xlabel='Diffusion axis 1',
                                  ylabel='Diffusion axis 2')
    dplot.savefig(dsg, join(FIGURES_DIR, 'figS12a1'), dpi=300, figsize=(4, 4))
    
    print('Rendering nodes coloring by E:72-83 genotypes')
    fig = dplot.dsg_to_fig(edges_dsg.opts(padding=0.05, aspect='square'))
    fig.set_size_inches(4, 4)
    axes = fig.axes[0]
    # fig, axes = plt.subplots(1, 1)
    mplot.plot_nodes(axes, nodes, x=x, y=y, color='c', palette=palette,
                     sort_by='function', sort_ascending=True, 
                     size=0.5, legend=False)
    axes.set(xlabel='Diffusion axis 1', ylabel='Diffusion axis 2')
    
    fpath = join(FIGURES_DIR, 'figS12a2')
    mplot.savefig(fig, fpath, fmt='png', dpi=300)
    
    print('Plotting histograms')
    fig, axes = plt.subplots(1, 1, figsize=(4, 4))
    
    idx = np.logical_or(nodes['E:83'] == 'Y', nodes['E:72'] == 'K')
    bins = np.linspace(0, 6, 25)
    sns.histplot(nodes.loc[idx, 'function'].values, 
                 label='E:Lys72 or E:Tyr83', color=cmap(0.6), bins=bins)
    sns.histplot(nodes.loc[~idx, 'function'].values, 
                 label='Other genotypes', color=palette['Other'], bins=bins)
    axes.legend(loc=1)
    axes.set(xlabel='Normalized MIC score', ylabel='# genotypes')
    sns.despine(ax=axes)
    fig.tight_layout()
    fig.savefig(join(FIGURES_DIR, 'figS12b.png'), format='png', dpi=300)
    fig.savefig(join(FIGURES_DIR, 'figS12b.svg'), format='svg', dpi=300)
