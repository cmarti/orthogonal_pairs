import seaborn as sns
import holoviews as hv
import gpmap.src.plot.ds as dplot

from os.path import join

from scripts.utils import load_visualization
from scripts.settings import (FIGURES_DIR, COL2, COL9,
                              RECOMBINANTS, CHANGE_POS, E2, E9, IM2, IM9)


def plot_landscape(nodes, edges, x, y, cmap):
    dsg = dplot.plot_visualization(nodes, x=x, y=y, edges_df=edges,
                                   nodes_cmap=cmap, nodes_resolution=1200,
                                   edges_resolution=1800, square=True)
    return(dsg)

def highlight_orthogonal_pairs(nodes, idx, x, y):
    nodes['orthogonal'] = 'No'
    nodes.loc[idx, 'orthogonal'] = 'Yes'
    ndf = nodes.loc[nodes['orthogonal'] == 'Yes', :].copy()
    ndf['v'] = 1
    dsg = dplot.plot_visualization(ndf, x=x, y=y, nodes_color='v', nodes_resolution=800,
                                   nodes_vmin=0, nodes_vmax=1, nodes_cmap='binary')
    return(dsg)


def to_fig(dsg, t1, t2):
    fig = dplot.dsg_to_fig(dsg.opts(padding=0.05))
    fig.set_size_inches(3.5, 3.5)
    axes = fig.axes[0]
    axes.set(title='Cognate > {} & non-cognate < {}'.format(t2, t1))
    axes.set_xlabel('Diffusion axis 1', fontsize=12)
    axes.set_ylabel('Diffusion axis 2', fontsize=12)
    sns.despine(ax=axes, right=False, top=False)
    fig.tight_layout()
    return(fig)


def get_pair2_orthogonal(nodes, threshold1, threshold2):
    idx1 = (nodes['function'] >= threshold2) & (nodes['E2/Im'] <= threshold1) & (nodes['E/Im2'] <= threshold1)
    return(idx1)


def get_pair9_orthogonal(nodes, threshold1, threshold2):
    idx2 = (nodes['function'] >= threshold2) & (nodes['E9/Im'] <= threshold1) & (nodes['E/Im9'] <= threshold1)
    return(idx2)


def main(nodes, edges):
    hv.extension('matplotlib')
    cmap = 'coolwarm'
    
    fmins = nodes.loc[[COL2, COL9], 'function'].values
    fmaxs = nodes.loc[RECOMBINANTS, 'function'].values
    print('Predicted normalized MIC scores in extant variants:')
    print('\tE2/Im2: {:.2f}'.format(fmins[0]))
    print('\tE2/Im2: {:.2f}'.format(fmins[1]))
    print('Predicted normalized MIC scores in recombinant variants:')
    print('\tE9/Im2: {:.2f}'.format(fmaxs[0]))
    print('\tE2/Im9: {:.2f}'.format(fmaxs[1]))
    pair9_thresholds = [(1, 1.5), (1.5, 2),   (2, 2),
                        (1, 2.5), (1.5, 2.5), (2, 2.5)]
    
    pair2_thresholds = [(1, 2.5), (1.5, 4.0), (2, 4.0),
                        (1, 5.5), (1.5, 5.5), (2, 5.5)]
    
    print('Annotating recombinants')
    nodes['E2/Im'] = nodes.reindex([x[:CHANGE_POS] + E2 for x in nodes.index])['function'].values
    nodes['E/Im2'] = nodes.reindex([IM2 + x[CHANGE_POS:] for x in nodes.index])['function'].values
    nodes['E9/Im'] = nodes.reindex([x[:CHANGE_POS] + E9 for x in nodes.index])['function'].values
    nodes['E/Im9'] = nodes.reindex([IM9 + x[CHANGE_POS:] for x in nodes.index])['function'].values
    
    x, y = '1', '2'
    
    # Main landscape
    print('Generating main visualization')
    dsg = plot_landscape(nodes, edges, x, y, cmap)
    
    # Pair2 orthogonal
    print('Overlaying pairs orthogonal to E2/Im2')
    for i, (t1, t2) in enumerate(pair2_thresholds):
        print('\t{}: Cognate > {} ; Non-cognate < {}'.format(i, t2, t1))
        idx = get_pair2_orthogonal(nodes, t1, t2)
        orth_dsg = highlight_orthogonal_pairs(nodes, idx, x, y)
        fig = to_fig(dsg * orth_dsg, t1, t2)
        
        fpath = join(FIGURES_DIR, 'figS16.{}.png'.format(i))
        fig.savefig(fpath, dpi=300)
        
    # Pair2 orthogonal
    print('Overlaying pairs orthogonal to E9/Im9')
    for i, (t1, t2) in enumerate(pair9_thresholds):
        print('\t{}: Cognate > {} ; Non-cognate < {}'.format(i, t2, t1))
        idx = get_pair9_orthogonal(nodes, t1, t2)
        orth_dsg = highlight_orthogonal_pairs(nodes, idx, x, y)
        fig = to_fig(dsg * orth_dsg, t1, t2)
        
        fpath = join(FIGURES_DIR, 'figS17.{}.png'.format(i))
        fig.savefig(fpath, dpi=300)
    

if __name__ == '__main__':
    nodes, edges = load_visualization()
    main(nodes, edges)
    