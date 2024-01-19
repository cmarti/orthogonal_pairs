import numpy as np

import gpmap.src.plot.ds as dplot
import gpmap.src.plot.mpl as mplot


from scripts.utils import load_visualization, plot_path, annotate_wts
from os.path import join
from scripts.settings import FIGURES_DIR, CMAP


def add_cbar_hist_inset(axes, values, cmap='viridis',
                        label='Function', fontsize=9, pos=(0.6, 0.5),
                        width=0.4, height=0.2, bins=50):
    vmin, vmax = values.min(), values.max()    
    ax1 = mplot.get_cbar_inset_axes(axes, pos=(pos[0], pos[1]), width=0.05,
                                    horizontal=True, height=width)
    mplot.draw_cbar(ax1, cmap=cmap, label=label, vmin=vmin, vmax=vmax, width=16,
              orientation='horizontal', fontsize=fontsize)
    ax2 = mplot.get_hist_inset_axes(axes, pos=(pos[0], pos[1] + height-0.085),
                                    width=width, height=height)
    mplot.plot_color_hist(ax2, values, cmap=cmap, bins=bins, fontsize=fontsize)


def main(nodes, edges, ylim=None):
    x, y = '1', '2'
    xticks = [-2, 2]
    yticks = np.array([-2, 2])
    dsg = dplot.plot_visualization(nodes, x=x, y=y, edges_df=edges,
                                   nodes_cmap=CMAP, nodes_resolution=1200,
                                   edges_resolution=1800, square=False,
                                   sort_by='function', sort_ascending=False)
    dsg.opts(padding=0.05, aspect='square')
    fig = dplot.dsg_to_fig(dsg)
    fig.set_size_inches(8, 8)
    axes = fig.axes[0]
    mplot.set_centered_spines(axes, xlabel='Diffusion axis {}'.format(x),
                              ylabel='Diffusion axis {}'.format(y), alpha=0.8,
                              zorder=-1, fontsize=12, add_grid=False,
                              xlabel_pos=(1.15, 0.15),
                              ylabel_pos=(0.1, 0.98))
    add_cbar_hist_inset(axes, nodes['function'],
                        cmap=CMAP, label='Normalized MIC score',
                        pos=(-0.05, 0.80), fontsize=9,
                        width=0.25, height=0.125, bins=20)

    # Path legend
    ypos = 0.95
    axes.text(0.03, ypos, 'Experimental\n        path',
              transform=axes.transAxes, fontsize=9,
              ha='left', va='center')
    axes.scatter(0.01, ypos, c='black', s=12,
                 transform=axes.transAxes)
    axes.plot((0.0, 0.02), (ypos, ypos), c='black', lw=0.8,
                 transform=axes.transAxes)
    
    plot_path(axes, nodes, x, y)
    annotate_wts(axes, nodes, x, y)    
    
#     axes.text(-0.15, 1., r'$\bf{A}$', transform=axes.transAxes, fontsize=14)

    axes.set(xticks=xticks, yticks=yticks,
             title='', xlabel='', ylabel='',
             ylim=ylim)
    
    fig.savefig(join(FIGURES_DIR, 'fig3a.png'), dpi=300)
    fig.savefig(join(FIGURES_DIR, 'fig3a.svg'), dpi=300, format='svg')


if __name__ == '__main__':
    nodes, edges = load_visualization()
    main(nodes, edges)
