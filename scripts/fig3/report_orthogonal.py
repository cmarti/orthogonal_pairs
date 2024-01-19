import pandas as pd
from scripts.utils import load_visualization
from scripts.settings import (COL2, COL9, RECOMBINANTS, CHANGE_POS, E2, E9, IM2, IM9,
    DATA_DIR, EXP_PATH)
from os.path import join


def get_orthogonal_pairs(nodes, p):
    fmins = nodes.loc[[COL2, COL9], 'function'].values
    fmaxs = nodes.loc[RECOMBINANTS, 'function'].values
    print(fmins, fmaxs)
    fmax = max(fmaxs)
    
    threshold1 = fmax + p
    threshold2 = fmins[1] - p
    # threshold1 = 1.35
    # threshold2 = 3.75
    # threshold1, threshold2 = 1.65, 5.25
    print('Pair2', threshold1, threshold2)
    idx1 = (nodes['function'] >= threshold2) & (nodes['E2/Im'] <= threshold1) & (nodes['E/Im2'] <= threshold1)
    
    nodes['col2_orthogonal'] = 'No'
    nodes.loc[idx1, 'col2_orthogonal'] = 'Yes'
    print(nodes.loc[EXP_PATH, ['function', 'E2/Im', 'E/Im2', 'col2_orthogonal']])
    
    
    threshold2 = fmins[0] - p
    # threshold1, threshold2 = 1.65, 1.7
    print('Pair9', threshold1, threshold2)    
    idx2 = (nodes['function'] >= threshold2) & (nodes['E9/Im'] <= threshold1) & (nodes['E/Im9'] <= threshold1)
    nodes['col9_orthogonal'] = 'No'
    nodes.loc[idx2, 'col9_orthogonal'] = 'Yes'
    print(nodes.loc[EXP_PATH, ['function', 'E9/Im', 'E/Im9', 'col9_orthogonal']])
    
    print('Col2 Orthogonal {} ({:.2f} %)'.format(idx1.sum(), idx1.mean() * 100))
    print('Col9 Orthogonal {} ({:.2f} %)'.format(idx2.sum(), idx2.mean() * 100))
    return(idx1, idx2)


def main(nodes):
    nodes['E2/Im'] = nodes.reindex([x[:CHANGE_POS] + E2 for x in nodes.index])['function'].values
    nodes['E/Im2'] = nodes.reindex([IM2 + x[CHANGE_POS:] for x in nodes.index])['function'].values
    nodes['E9/Im'] = nodes.reindex([x[:CHANGE_POS] + E9 for x in nodes.index])['function'].values
    nodes['E/Im9'] = nodes.reindex([IM9 + x[CHANGE_POS:] for x in nodes.index])['function'].values
    get_orthogonal_pairs(nodes, p=0.5)


if __name__ == '__main__':
    data = pd.read_csv(join(DATA_DIR, 'data.calibrated_predictions.csv'), index_col=0).set_index('seq')
    data['function'] = data['y']
    print(data)
    main(data)
    exit()
    
    
    nodes = load_visualization(load_edges=False)
    nodes['>3'] = nodes['function'] > 3
    
    print('Neutrality')
    print(nodes['function'].mean())
    f1 = nodes['>3'].mean()
    print(f1)
    print('Selection')
    print((nodes['function'] * nodes['stationary_freq']).sum())
    f2 = (nodes['>3'] * nodes['stationary_freq']).sum()
    print(f2)
    print(f2 / f1)
    main(nodes)
    
    
    
    
    
    
    