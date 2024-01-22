import pandas as pd
import seaborn as sns
import gpmap.src.plot.ds as dplot

from os.path import join
from gpmap.src.utils import read_edges

from scripts.utils import annotate_wts, plot_path
from scripts.settings import VIZ_DIR, FIGURES_DIR


if __name__ == '__main__':
    model = 'model2'
    x, y = '1', '2'
    
    edges_fpath = join(VIZ_DIR, 'visualization.edges.npz')
    edges = read_edges(edges_fpath)
    
    for i, sd in enumerate([0.5, 1, 1.5, 2, 2.5, 3]):
        print('Plotting visualization with noised Normalized MIC score SD={}'.format(sd))
        nodes_fpath = join(VIZ_DIR, 'noised_landscape.sd{}.nodes.pq'.format(sd))
        nodes = pd.read_parquet(nodes_fpath)
        nodes['1'] = -nodes['1']
        
        dsg = dplot.plot_visualization(nodes, x=x, y=y, edges_df=edges,
                                       nodes_cmap='coolwarm', nodes_resolution=1200,
                                       edges_resolution=1800, square=False)
        fig = dplot.dsg_to_fig(dsg.opts(padding=0.15))
        fig.set_size_inches(4, 4)
        axes = fig.axes[0]

        plot_path(axes, nodes, x, y)        
        annotate_wts(axes, nodes, x, y)
        axes.set(ylabel='Diffusion axis 2', xlim=(-3.5, 3.5), ylim=(-3, 3),
                 title=r'$\Delta\Delta G_{binding}$ ' + 'Noise SD={} (R.e.u.)'.format(sd))
        sns.despine(ax=axes, right=False, top=False)
        
        fig.tight_layout()
        fig.savefig(join(FIGURES_DIR, 'figS15.{}.sd{}.png'.format(i, sd)), dpi=300)
