import numpy as np
import pandas as pd
import matplotlib as mpl
import gpmap.src.plot.mpl as mplot

from os.path import join

from gpmap.src.utils import read_dataframe, read_edges
from scripts.settings import (COL2, COL9, EXP_DATA_FPATH, VIZ_DIR, EXP_PATH,
                              EXP_EXTRA_DATA_FPATH, OLD_EXP_DATA_FPATH,
    POS_LABELS)
from gpmap.src.space import SequenceSpace
from gpmap.src.genotypes import select_genotypes
from itertools import combinations, product
from _collections import defaultdict

def load_exp_data(energies=None, new=True):
    if new:
        data = pd.read_csv(EXP_DATA_FPATH,
                           usecols=['mut_seq', 'd_tot', 'd_ddg', 'mic1', 'mic2',
                                    'background', 'mic_avg', 'viability_score'])
        data['basal_mic'] = data['mic_avg'] + data['viability_score']
        data['norm_mic1'] = data['basal_mic'] - data['mic1']
        data['norm_mic2'] = data['basal_mic'] - data['mic2']
        data['err'] = np.abs(data['mic1'] - data['mic2'])
        data['d2'] = [np.sum([x != y for x, y in zip(seq, COL2)]) for seq in data['mut_seq']]
        data['d9'] = [np.sum([x != y for x, y in zip(seq, COL9)]) for seq in data['mut_seq']]
        data['d'] = np.min(data[['d2', 'd9']], axis=1)
        data['y_std'] = data[['norm_mic1', 'norm_mic2']].std(1)
    else:
        data = pd.read_csv(OLD_EXP_DATA_FPATH,
                           usecols=['mut_seq', 'background', 'viability_score'])
    
    if energies is not None:
        data['col2_binding'] = energies.loc[data['mut_seq'], 'col2_binding'].values
        data['col9_binding'] = energies.loc[data['mut_seq'], 'col9_binding'].values
    return(data)


def load_extra_exp_data(energies):
    data = pd.read_csv(EXP_EXTRA_DATA_FPATH,
                       usecols=['muts_seq', 'viability_score']).set_index('muts_seq')
    
    if energies is not None:
        data['col2_binding'] = energies.loc[data.index, 'col2_binding'].values
        data['col9_binding'] = energies.loc[data.index, 'col9_binding'].values
        
        if 'phi' in energies.columns:
            data['phi'] = energies.loc[data.index, 'phi'].values
            data['y'] = energies.loc[data.index, 'y'].values

    return(data)


def load_visualization(model=None, m=None, flip_axis=['1', '2'], load_edges=True, test=False):
    if model is None:
        nodes_fpath = join(VIZ_DIR, 'visualization.nodes.pq')
    else:
        nodes_fpath = join(VIZ_DIR, '{}.{}.nodes.pq'.format(model, m))
    
    if test:
        nodes_fpath = join(VIZ_DIR, 'test.nodes.pq')
    nodes = read_dataframe(nodes_fpath)
    
    for axis in flip_axis:
        nodes[axis] = -nodes[axis]
        
    res = nodes
    if load_edges:
        edges_fpath = join(VIZ_DIR, 'edges.npz')
        if test:
            edges_fpath = join(VIZ_DIR, 'test.edges.npz')
        edges = read_edges(edges_fpath)
        res = nodes, edges
    
    return(res)


def set_plotting_params():
    mpl.rcParams['ytick.labelsize'] = 5
    mpl.rcParams['xtick.labelsize'] = 5
    mpl.rcParams['xtick.major.width'] = 0.5
    mpl.rcParams['ytick.major.width'] = 0.5
    mpl.rcParams['xtick.major.size'] = 3
    mpl.rcParams['ytick.major.size'] = 3
    mpl.rcParams['axes.labelsize'] = 5
    mpl.rcParams['axes.titlesize'] = 5
    mpl.rcParams['font.sans-serif'] = ['Helvetica']
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['axes.linewidth'] = 0.5


def plot_path(axes, nodes, x, y, s=12, zorder=3, edges_width=0.8, color='black',
              cmap='coolwarm', nodes_vmin=0, nodes_vmax=6):
    path = nodes.loc[EXP_PATH, :]
    n = path.shape[0]
    path_edges = pd.DataFrame({'i': np.arange(n-1), 'j': np.arange(1, n)})
    mplot.plot_visualization(axes, path, edges_df=path_edges,
                             x=x, y=y, nodes_color=color,
                             nodes_cmap=cmap, nodes_vmin=nodes_vmin,
                             nodes_vmax=nodes_vmax,
                             edges_alpha=1, nodes_lw=edges_width,
                             nodes_zorder=zorder+1, edges_zorder=zorder,
                             nodes_cbar=False, nodes_size=s,
                             edges_color='black', edges_width=edges_width)

    
def annotate_wts(axes, nodes, x, y):
    xy = (nodes.loc[COL2, [x, y]])  
    axes.annotate('E2/Im2', xy=xy,
                  xytext=[xy[0] - 0.1, xy[1]+0.1],
                  fontsize=10, ha='right', va='bottom',
                  arrowprops={'arrowstyle':'-|>'})
    xy = (nodes.loc[COL9, [x, y]])
    axes.annotate('E9/Im9', xy=xy,
                  xytext=[xy[0] + 0.1, xy[1] + 0.1],
                  fontsize=10, ha='left', va='bottom',
                  arrowprops={'arrowstyle':'-|>'})
    

def get_extant_combinations(positions=None):
    templates = [COL2, COL9]
    if positions is not None:
        templates = get_subsequences(positions, templates)
    seqs = np.array([''.join([templates[t][i] for i, t in enumerate(seq)]) 
                     for seq in product([0, 1], repeat=len(templates[0]))])
    return(seqs)

    
def get_subsequences(pos_labels, seqs):
    molecule_pos = {'E': [], 'Im': []}
    for p in sorted(pos_labels):
        molecule_pos[p.split(':')[0]].append(POS_LABELS.index(p))
    return(np.array(['{}/{}'.format(''.join([x[i] for i in molecule_pos['E']]), 
                                    ''.join([x[i] for i in molecule_pos['Im']])).strip('/')
                     for x in seqs]))

def format_pos(pos_labels):
    molecule_pos = defaultdict(list)
    for label in pos_labels:
        molecule, pos = label.split(':')
        molecule_pos[molecule].append(pos)
    molecule_pos = sorted(['{}:{}'.format(k, '-'.join(v)) 
                           for k, v in molecule_pos.items()
                           if len(v) > 0])
    return('/'.join(molecule_pos))


def format_subseq(pos_labels, subseq):
    molecule_pos = {'E': [], 'Im': []}
    for p, a in zip(pos_labels, subseq.replace('/', '')):
        molecule, pos = p.split(':')
        molecule_pos[molecule].append('{}{}'.format(pos, a))
    molecule_pos = sorted(['{}:{}'.format(k, '-'.join(v)) 
                           for k, v in molecule_pos.items()
                           if len(v) > 0])
    return('/'.join(molecule_pos))


def format_mut(pos, a1, a2):
    molecule, pos = pos.split(':')
    return('{}:{}{}{}'.format(molecule, a1, pos, a2))


def path_to_mut(path):
    muts = []
    for s1, s2 in zip(path, path[1:]):
        pos = [i for i, (a1, a2) in enumerate(zip(s1, s2)) if a1 != a2]
        if len(pos) > 1:
            msg = 'Sequences {} and {} differ at more than 1 position'.format(s1, s2)
            raise ValueError(msg)
        p = pos[0]
        pos_label = POS_LABELS[p]
        mut = format_mut(pos_label, s1[p], s2[p])
        muts.append(mut)
    return(muts)


def calc_lims(v, p=0.1):
    vmax, vmin = np.max(v), np.min(v)
    dv = vmax - vmin
    return(vmin - p * dv, vmax + p * dv)


def calc_mutational_effect(df, pos, a1, a2, y='y1'):
    i = POS_LABELS.index(pos)
    seqs1 = [x[:i] + a1 + x[i+1:] for x in df.index.values] 
    seqs2 = [x[:i] + a2 + x[i+1:] for x in df.index.values]
    return(df.loc[seqs2, y].values - df.loc[seqs1, y].values)


def annotate_landscape(nodes, sep='/'):
    nodes['E:98/Im:33-34'] = [x[16] + sep + x[6:8] for x in nodes.index]
    nodes['E:97/Im:38'] = [x[15] + sep + x[8] for x in nodes.index]
    nodes['E:83-98/Im:33-34'] = [x[14] + x[16] + sep + x[6:8] for x in nodes.index]
    
    nodes['E:98'] = [x[16] for x in nodes.index]
    nodes['Im:33'] = [x[6] for x in nodes.index]
    nodes['Im:34'] = [x[7] for x in nodes.index]
    
    nodes['E:72'] = [x[10] for x in nodes.index]
    nodes['E:73'] = [x[11] for x in nodes.index]
    nodes['E:77'] = [x[12] for x in nodes.index]
    
    # Secondary interactions
    nodes['E:83'] = [x[14] for x in nodes.index]
    nodes['E:78'] = [x[13] for x in nodes.index]
    nodes['E:83-98/Im:34'] = [x[14] + x[16] + sep + x[7] for x in nodes.index]
    nodes['E:83-98'] = [x[14] + x[16] for x in nodes.index]
    nodes['E:78-83-98'] = [x[13:15] + x[16] for x in nodes.index]
    nodes['E:73-83'] = [x[11] + x[14] for x in nodes.index]
    nodes['E:72-73-77'] = [x[10:13] for x in nodes.index]
    nodes['E:73-78-83'] = [x[11] + x[13:15] for x in nodes.index]
    nodes['E:78-97/Im:38'] = [x[13] + x[15] + sep + x[8] for x in nodes.index]
    nodes['E:72-77-83'] = [x[10] + x[12] + x[14] for x in nodes.index]
    
    nodes['E:72-73'] = [x[10:12] for x in nodes.index]
    
    # Contexts
    nodes['E:97-98/Im:33-34-38'] = [x[15:] + sep + x[6:9] for x in nodes.index]
    nodes['context2'] = [x[:10] + x[13:] for x in nodes.index]
    
    nodes['83K'] = nodes.loc[[x[:14] + 'K' + x[15:] for x in nodes.index], 'function'].values
    nodes['83Y'] = nodes.loc[[x[:14] + 'Y' + x[15:] for x in nodes.index], 'function'].values
    nodes['K83Y'] = nodes['83Y'] - nodes['83K']


def rm_alleles(ndf, rm_alleles):
    for pos, allele in rm_alleles.items():
        keep_idx = np.array([seq[pos] != allele for seq in ndf.index])
        ndf = ndf.loc[keep_idx, :]
    return(ndf)

def plot_seq_names(ndf, axes, x='1', y='function', fontsize=10):
    for x, y, label in zip(ndf[x].values, ndf[y].values, ndf.index.values):
        axes.text(x+np.random.normal(0, 0.01), y+np.random.normal(0, 0.01), label,
                  fontsize=fontsize)
        
        
def get_marginal_landscape(nodes, path, interaction, alleles_rm={}):
    columns = ['function', '1', '2', '3']
    ndf = nodes.groupby(interaction)[columns].mean().reset_index().set_index(interaction)
    if alleles_rm:
        ndf = rm_alleles(ndf, alleles_rm)
    edf = get_edges_df(ndf.index)
    ndfp, edfp = select_genotypes(ndf, genotypes=path, edges=edf)
    return(ndf, edf, ndfp, edfp)


def plot_marginal_landscape(axes, ndf, edf, ndfp, edfp,
                            xlim=None, ylim=None, 
                            seed=0, x='1'):
    np.random.seed(seed)
    mplot.plot_visualization(axes, ndf, edges_df=edf, x=x, y='function',
                            nodes_size=20, nodes_color='black',
                            edges_alpha=1, edges_width=0.5)
    plot_seq_names(ndf, axes, x=x)
    mplot.plot_visualization(axes, ndfp, edges_df=edfp, x=x, y='function',
                             nodes_size=10, nodes_color='black',
                             edges_alpha=1, edges_color='black', edges_width=1)
    axes.set_ylabel('Average normalized MIC score', fontsize=12)
    axes.set_xlabel('Diffusion axis {}'.format(x), fontsize=12)
    axes.set(xlim=xlim, ylim=ylim)


def get_edges_df(seqs):
    i, j = [], []
    idxs = np.arange(seqs.shape[0])
    for i1, i2 in combinations(idxs, 2):
        seq1, seq2 = seqs[i1], seqs[i2]
        d = np.sum([a1 != a2 for a1, a2 in zip(seq1, seq2)])
        if d == 1:
            i.append(i1)
            j.append(i2)
    return(pd.DataFrame({'i': i, 'j': j}))
