import numpy as np
import pandas as pd
import seaborn as sns
import holoviews as hv
import networkx as nx

import gpmap.src.plot.ds as dplot
import gpmap.src.plot.mpl as mplot


from scripts.utils import load_visualization
from os.path import join
from scripts.settings import (FIGURES_DIR, COL2, COL9,
                              RECOMBINANTS, CHANGE_POS, E2, E9, IM2, IM9,
                              EXP_PATH)
from gpmap.src.genotypes import select_genotypes


def add_cbar_hist_inset(axes, values, cmap='viridis',
                        label='Function', fontsize=9, pos=(0.6, 0.5),
                        width=0.4, height=0.2, bins=50):
    vmin, vmax = values.min(), values.max()    
    ax1 = mplot.get_cbar_inset_axes(axes, pos=(pos[0], pos[1]), width=0.05,
                                    horizontal=True, height=width)
    mplot.draw_cbar(ax1, cmap=cmap, label=label, vmin=vmin, vmax=vmax, width=16,
              orientation='horizontal', fontsize=fontsize)
    ax2 = mplot.get_hist_inset_axes(axes, pos=(pos[0], pos[1] + height-0.06),
                                    width=width, height=height)
    mplot.plot_color_hist(ax2, values, cmap=cmap, bins=bins, fontsize=fontsize)


def get_orthogonal_pairs(nodes, p=0.25):
    fmins = nodes.loc[[COL2, COL9], 'function'].values
    fmaxs = nodes.loc[RECOMBINANTS, 'function'].values
    print(fmins, fmaxs)
    fmax = max(fmaxs)
    
    threshold1 = fmax + p
    threshold2 = fmins[1] - p
    threshold1, threshold2 = 1.65, 5.25
    print('Pair2', threshold1, threshold2)
    idx1 = (nodes['function'] >= threshold2) & (nodes['E2/Im'] <= threshold1) & (nodes['E/Im2'] <= threshold1)
    
    threshold2 = fmins[0] - p
    threshold1, threshold2 = 1.65, 1.70
    print('Pair9', threshold1, threshold2)    
    idx2 = (nodes['function'] >= threshold2) & (nodes['E9/Im'] <= threshold1) & (nodes['E/Im9'] <= threshold1)
    
    print('Col2 Orthogonal {} ({:.2f} %)'.format(idx1.sum(), idx1.mean() * 100))
    print('Col9 Orthogonal {} ({:.2f} %)'.format(idx2.sum(), idx2.mean() * 100))
    return(idx1, idx2)


def add_legend(axes, label):
    axes.text(0.82, 0.08, label,
              transform=axes.transAxes, fontsize=9,
              ha='center', va='center')
    axes.scatter(0.67, 0.08, c='black', s=5, transform=axes.transAxes)


def plot_path(axes, nodes, x, y, s=12, zorder=3, edges_width=0.8, color='col2_orthogonal'):
    path = nodes.loc[EXP_PATH, :]
    palette = {'Yes': 'black', 'No': 'darkgrey'}
    print(path[['function', color]])
    n = path.shape[0]
    path_edges = pd.DataFrame({'i': np.arange(n-1), 'j': np.arange(1, n)})
    path_edges['c'] = 'No'
    path['function'].values
    idx = np.logical_and(path[color].values[:-1] == 'Yes', path[color].values[1:] == 'Yes')
    path_edges.loc[idx, 'c'] = 'Yes'
    
    mplot.plot_visualization(axes, path, edges_df=path_edges,
                             x=x, y=y, nodes_color=color,
                             nodes_palette=palette,
                             edges_alpha=1,
                             nodes_zorder=zorder+1, edges_zorder=zorder,
                             nodes_cbar=False, nodes_size=s,
                             edges_color='c', edges_palete=palette,
                             edges_width=edges_width)
    axes.legend().set_visible(False)


def plot_orthogonality(nodes, edges, x, y, cmap, idx1, idx2):
    dsg = dplot.plot_visualization(nodes, x=x, y=y, edges_df=edges,
                                   nodes_cmap=cmap, nodes_resolution=1200,
                                   edges_resolution=1800, square=True)
    
    # Col2 orthogonal
    nodes['col2_orthogonal'] = 'No'
    nodes.loc[idx1, 'col2_orthogonal'] = 'Yes'
    ndf = nodes.loc[nodes['col2_orthogonal'] == 'Yes', :].copy()
    print(nodes.loc[EXP_PATH, ['function', 'E2/Im', 'E/Im2', 'col2_orthogonal']])
    
    ndf, edf = select_genotypes(nodes, ndf.index.values, edges=edges)
    g = nx.Graph()
    g.add_edges_from(edf.values)
    connected = nx.is_connected(g)
    print('Are sequences orthogonal to E2/Im2 connected? {}'.format(connected))
    n = len(list(nx.connected_components(g)))
    print('\tThere are {} connected components'.format(n))
    
    
    ndf['v'] = 1
    dsg1 = dsg * dplot.plot_visualization(ndf, x=x, y=y, nodes_color='v', nodes_resolution=800,
                                          nodes_vmin=0, nodes_vmax=1, nodes_cmap='binary')
    fig = dplot.dsg_to_fig(dsg1.opts(padding=0.05))
    fig.set_size_inches(3.5, 3.5)
    axes = fig.axes[0]
    # plot_path(axes, nodes, x, y, s=2, edges_width=0.4, color='col2_orthogonal')
    axes.set(title='') #, xlabel='', ylabel='')
    axes.set_xlabel('Diffusion axis 1', fontsize=12)
    axes.set_ylabel('Diffusion axis 2', fontsize=12)
    axes.text(-0.2, 1.0, r'$\bf{E}$', transform=axes.transAxes, fontsize=14)
    sns.despine(ax=axes, right=False, top=False)
    add_legend(axes, label='Orthogonal\nto E2/Im2')
    fig.tight_layout()
    fig.savefig(join(FIGURES_DIR, 'fig3e.png'), dpi=300)
    fig.savefig(join(FIGURES_DIR, 'fig3e.svg'), dpi=300, format='svg')
    
    # Col9 orthogonal
    nodes['col9_orthogonal'] = 'No'
    nodes.loc[idx2, 'col9_orthogonal'] = 'Yes'
    print(nodes.loc[EXP_PATH, ['function', 'E9/Im', 'E/Im9', 'col9_orthogonal']])
    
    ndf = nodes.loc[nodes['col9_orthogonal'] == 'Yes', :].copy()
    
    ndf, edf = select_genotypes(nodes, ndf.index.values, edges=edges)
    g = nx.Graph()
    g.add_edges_from(edf.values)
    connected = nx.is_connected(g)
    print('Are sequences orthogonal to E9/Im9 connected? {}'.format(connected))
    n = len(list(nx.connected_components(g)))
    print('\tThere are {} connected components'.format(n))
    
    ndf['v'] = 1
    dsg2 = dsg * dplot.plot_visualization(ndf, x=x, y=y, nodes_color='v', nodes_resolution=800,
                                          nodes_vmin=0, nodes_vmax=1, nodes_cmap='binary')
    fig = dplot.dsg_to_fig(dsg2.opts(padding=0.05))
    fig.set_size_inches(3.5, 3.5)
    axes = fig.axes[0]
    # plot_path(axes, nodes, x, y, s=5, edges_width=0.4, color='col9_orthogonal')
    axes.set_title('')
    axes.set_xlabel('Diffusion axis 1', fontsize=12)
    axes.set_ylabel('Diffusion axis 2', fontsize=12)
    sns.despine(ax=axes, right=False, top=False)
    axes.text(-0.2, 1.0, r'$\bf{F}$', transform=axes.transAxes, fontsize=14,
              ha='right', va='bottom')
    add_legend(axes, label='Orthogonal\nto E9/Im9')
    
    fig.tight_layout()
    fig.savefig(join(FIGURES_DIR, 'fig3f.png'), dpi=300)
    fig.savefig(join(FIGURES_DIR, 'fig3f.svg'), dpi=300, format='svg')


def main(nodes, edges):
    hv.extension('matplotlib')
    cmap = 'coolwarm'
    
    nodes['E2/Im'] = nodes.reindex([x[:CHANGE_POS] + E2 for x in nodes.index])['function'].values
    nodes['E/Im2'] = nodes.reindex([IM2 + x[CHANGE_POS:] for x in nodes.index])['function'].values
    nodes['E9/Im'] = nodes.reindex([x[:CHANGE_POS] + E9 for x in nodes.index])['function'].values
    nodes['E/Im9'] = nodes.reindex([IM9 + x[CHANGE_POS:] for x in nodes.index])['function'].values
    
    x, y = '1', '2'
    idx1, idx2 = get_orthogonal_pairs(nodes, p=0.5)
    plot_orthogonality(nodes, edges, x, y, cmap, idx1, idx2)


if __name__ == '__main__':
    nodes, edges = load_visualization()
    main(nodes, edges)
    