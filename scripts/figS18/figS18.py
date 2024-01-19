import seaborn as sns
import holoviews as hv
import gpmap.src.plot.ds as dplot

from os.path import join

from scripts.utils import load_visualization
from scripts.settings import CHANGE_POS, EXP_PATH, FIGURES_DIR, POS_LABELS


def to_fig(dsg, title=''):
    fig = dplot.dsg_to_fig(dsg.opts(padding=0.05))
    fig.set_size_inches(3.5, 3.5)
    axes = fig.axes[0]
    axes.set(title=title)
    axes.set_xlabel('Diffusion axis 1', fontsize=12)
    axes.set_ylabel('Diffusion axis 2', fontsize=12)
    sns.despine(ax=axes, right=False, top=False)
    fig.tight_layout()
    return(fig)


def get_molecules_paths():
    ims, es = [], []
    ims_mut, es_mut = [], []
    
    for seq in EXP_PATH:
        im, e = seq[:CHANGE_POS], seq[CHANGE_POS:]
        
        if not es or e != es[-1]:
            if not es:
                es_mut.append('')
            else:
                es_mut.append(';'.join([' (E:{}{}{})'.format(a1, POS_LABELS[CHANGE_POS:][i].split(':')[1], a2)
                                        for i, (a1, a2) in enumerate(zip(es[-1], e))
                                        if a1 != a2]))
            es.append(e)
            
        if not ims or im != ims[-1]:
            if not ims:
                ims_mut.append('')
            else:
                ims_mut.append(';'.join([' (Im:{}{}{})'.format(a1, POS_LABELS[:CHANGE_POS][i].split(':')[1], a2)
                                         for i, (a1, a2) in enumerate(zip(ims[-1], im))
                                         if a1 != a2]))
            ims.append(im)
            
    return(ims, es, ims_mut, es_mut)


def main(nodes, edges):
    hv.extension('matplotlib')
    cmap = 'coolwarm'
    
    x, y = '1', '2'
    
    # Main landscape
    print('Generating edges background')
    edges_dsg = dplot.plot_edges(nodes, x=x, y=y, edges_df=edges,
                                 resolution=1800, square=True)
    ims, es, ims_mut, es_mut = get_molecules_paths()
    viable = (nodes['function'] >= 1).astype(float).values
    

    print('Generating specificity maps for Im')
    for i, (im, mut) in enumerate(zip(ims, ims_mut)):
        print('\t{}: {}'.format(i, im + mut))
        nodes['f'] = nodes.loc[[im + x[CHANGE_POS:] for x in nodes.index], 'function'].values * viable
        nodes_dsg = dplot.plot_nodes(nodes, x=x, y=y, cmap=cmap, color='f',
                                     resolution=800, square=True,
                                     sort_by='function', sort_ascending=False)
        fig = to_fig(edges_dsg * nodes_dsg, title='Im:{}{}'.format(im, mut))
        fpath = join(FIGURES_DIR, 'figS18.Im_specificity.{}.png'.format(i))
        fig.savefig(fpath, dpi=300)
    
    print('Generating specificity maps for E')    
    for i, (e, mut) in enumerate(zip(es, es_mut)):
        print('\t{}: {}'.format(i, e + mut))
        nodes['f'] = nodes.loc[[x[:CHANGE_POS] +e for x in nodes.index], 'function'].values * viable
        nodes_dsg = dplot.plot_nodes(nodes, x=x, y=y, cmap=cmap, color='f',
                                     resolution=800, square=True,
                                     sort_by='function', sort_ascending=False)
        fig = to_fig(edges_dsg * nodes_dsg, title='E:{}{}'.format(e, mut))
        fpath = join(FIGURES_DIR, 'figS18.E_specificity.{}.png'.format(i))
        fig.savefig(fpath, dpi=300)
        

if __name__ == '__main__':
    nodes, edges = load_visualization()
    main(nodes, edges)
    