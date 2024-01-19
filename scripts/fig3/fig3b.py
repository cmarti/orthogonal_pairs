import matplotlib.pyplot as plt

from os.path import join

from scripts.settings import FIGURES_DIR
from scripts.utils import (load_visualization, get_marginal_landscape,
                           plot_marginal_landscape, get_subsequences,
                           format_pos, format_subseq, calc_lims)


if __name__ == '__main__':
    panel_label=r'$\bf{B}$'
    ymax = 1.8
    p = 0.15
    axis = '1'
    interaction = ['E:97', 'Im:38']
    path = ['E/R', 'K/R', 'K/T']
    
    # Get labels and positions
    interaction_label = format_pos(interaction)
    
    # Get sublandscape
    nodes = load_visualization(load_edges=False)
    nodes[interaction_label] = get_subsequences(interaction, nodes.index.values)
    x = get_marginal_landscape(nodes, path=path, interaction=interaction_label)
    ndf, edf, ndfp, edfp = x
    
    # Make plot
    fig, axes = plt.subplots(1, 1, figsize=(3.75, 3.5))
    plot_marginal_landscape(axes, ndf, edf, ndfp, edfp,
                            ylim=(0, ymax), xlim=calc_lims(ndf[axis], p=p),
                            seed=0, x=axis)
    axes.text(-0.25, 1.0, panel_label, transform=axes.transAxes, fontsize=14,
              ha='right', va='bottom')
    axes.set_title('')
    axes.set_yticks([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8])
    
    fig.tight_layout()
    fpath = join(FIGURES_DIR, 'fig3b')
    fig.savefig(fpath + '.svg', dpi=300)
    fig.savefig(fpath + '.png', dpi=300)
    