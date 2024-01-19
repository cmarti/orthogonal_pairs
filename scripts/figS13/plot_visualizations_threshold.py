import pandas as pd
import seaborn as sns
import gpmap.src.plot.ds as dplot

from os.path import join
from gpmap.src.utils import read_edges

from scripts.utils import annotate_wts, plot_path
from scripts.settings import VIZ_DIR, FIGURES_DIR


if __name__ == '__main__':
    cmap = 'coolwarm'
    x, y = '1', '2'
    
    edges_fpath = join(VIZ_DIR, 'edges.npz')
    edges = read_edges(edges_fpath)
    
    for i, t in enumerate([1, 1.25, 1.5, 1.75, 2, 2.25]):
        print('Plotting visualization with Normalized MIC score > {}'.format(t))
        nodes_fpath = join(VIZ_DIR, 'threshold_landscape.t{}.nodes.pq'.format(t))
        nodes = pd.read_parquet(nodes_fpath)
        nodes['1'] = -nodes['1']
        
        dsg = dplot.plot_visualization(nodes, x=x, y=y, edges_df=edges,
                                       nodes_cmap=cmap, nodes_resolution=1200,
                                       edges_resolution=1800, square=False)
        fig = dplot.dsg_to_fig(dsg.opts(padding=0.15))
        fig.set_size_inches(4, 4)
        axes = fig.axes[0]

        plot_path(axes, nodes, x, y)        
        annotate_wts(axes, nodes, x, y)
        axes.set(ylabel='Diffusion axis 2', xlim=(-4.5, 2.5), ylim=(-3, 3),
                 title='Normalized MIC score > {}'.format(t))
        sns.despine(ax=axes, right=False, top=False)
        
        fig.tight_layout()
        fig.savefig(join(FIGURES_DIR, 'figS13.{}.t{}.png'.format(i, t)), dpi=300)
