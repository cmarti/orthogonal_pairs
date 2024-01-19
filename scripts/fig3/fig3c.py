import matplotlib.pyplot as plt

from os.path import join

from scripts.settings import FIGURES_DIR
from scripts.utils import (load_visualization, get_marginal_landscape,
                           plot_marginal_landscape, get_subsequences,
                           format_pos, format_subseq, calc_lims)


if __name__ == '__main__':
    positions=['E:83']
    allele = 'K'
    panel_label=r'$\bf{C}$'
    p = 0.15
    axis = '1'
    interaction = ['E:98', 'Im:33', 'Im:34']
    path = ['R/DN', 'R/DV', 'M/DV', 'M/LV']
    alleles_rm = {0: 'V', 2: 'V', 3: 'D'}
    
    # Get labels and positions
    label = format_subseq(positions, allele)
    pos_label = format_pos(positions)
    interaction_label = format_pos(interaction)
    
    # Get sublandscape
    nodes = load_visualization(load_edges=False)
    nodes[pos_label] = get_subsequences(positions, nodes.index.values)
    ndf = nodes.loc[nodes[pos_label] == allele].copy()
    ndf[interaction_label] = get_subsequences(interaction, ndf.index.values)
    x = get_marginal_landscape(ndf, path=path, interaction=interaction_label,
                               alleles_rm=alleles_rm)
    ndf, edf, ndfp, edfp = x
    
    # Make plot
    fig, axes = plt.subplots(1, 1, figsize=(3.9, 3.5))
    plot_marginal_landscape(axes, ndf, edf, ndfp, edfp,
                            ylim=calc_lims(ndf['function'], p=p),
                            xlim=calc_lims(ndf[axis], p=p),
                            seed=0, x=axis)
    axes.text(-0.25, 1.0, panel_label, transform=axes.transAxes, fontsize=14,
              ha='right', va='bottom')
    axes.text(0.97, 0.07, label, transform=axes.transAxes, fontsize=10,
              ha='right', va='top')
    axes.set_title('')
    axes.set_yticks([0, 0.3, 0.6, 0.9, 1.2])
    
    fig.tight_layout()
    fpath = join(FIGURES_DIR, 'fig3c')
    fig.savefig(fpath + '.svg', dpi=300)
    fig.savefig(fpath + '.png', dpi=300)
    