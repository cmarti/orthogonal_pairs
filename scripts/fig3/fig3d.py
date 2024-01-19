import matplotlib as mpl
import seaborn as sns

from os.path import join
from scripts.settings import FIGURES_DIR
from scripts.utils import (load_visualization, get_marginal_landscape,
                           plot_marginal_landscape,
                           get_subsequences, calc_lims)


def make_plot(axes, nodes, axis, label, ylim):
    x = get_marginal_landscape(nodes, path=['KA', 'NA', 'NP'], 
                               interaction='E:72-73')
    ndf, edf, ndfp, edfp = x
    xlim = calc_lims(ndf[axis], p=0.2)
    
    # Make plot
    plot_marginal_landscape(axes, ndf, edf, ndfp, edfp,
                            xlim=xlim, ylim=ylim,
                            seed=0, x=axis)
    axes.set_title('')
    axes.text(0.97, 0.07, label, transform=axes.transAxes, fontsize=10,
              ha='right', va='top')


def main(nodes, axis='1', ylim=(-0.3, 3.2)):
    nodes['context'] = get_subsequences(['E:77', 'E:83'], nodes.index.values)
    nodes['E:72-73'] = get_subsequences(['E:72', 'E:73'], nodes.index.values)
    fig, subplots = mpl.pyplot.subplots(1, 3, figsize=(7.25, 3.5))
    
    axes = subplots[0]
    label = 'E:77T-83K'
    ndf = nodes.loc[nodes['context'] == 'TK', :]
    make_plot(axes, ndf, axis, label, ylim)
    axes.text(-0.3, 1.0, r'$\bf{D}$', transform=axes.transAxes, fontsize=14,
              ha='right', va='bottom')
    sns.despine(ax=axes, top=False, right=True)
    
    axes = subplots[1]
    label = 'E:77T-83Y'
    ndf = nodes.loc[nodes['context'] == 'TY', :]
    make_plot(axes, ndf, axis, label, ylim)
    axes.set(ylabel='', yticks=[])
    sns.despine(ax=axes, top=False, right=True, left=True)
    
    axes = subplots[2]
    label = 'E:77S-83Y'
    ndf = nodes.loc[nodes['context'] == 'SY', :]
    make_plot(axes, ndf, axis, label, ylim)
    axes.set(ylabel='', yticks=[])
    sns.despine(ax=axes, top=False, right=False, left=True)
    
    fig.tight_layout()
    fpath = join(FIGURES_DIR, 'fig3d.')
    fig.savefig(fpath + 'svg', dpi=300)
    fig.savefig(fpath + 'png', dpi=300)


if __name__ == '__main__':
    nodes = load_visualization(load_edges=False)
    main(nodes)
    