import pandas as pd

from os.path import join
from scripts.settings import DATA_DIR, VIZ_DIR
from gpmap.src.space import SequenceSpace
from gpmap.src.randwalk import WMWalk


if __name__ == '__main__':
    fpath = join(DATA_DIR, 'calibrated_predictions.pq')
    print('Reading data from {}'.format(fpath))
    data = pd.read_parquet(fpath)
    
    odds_ratio = 10
    for x in [1, 1.25, 1.5, 1.75, 2, 2.25]:
        print('Calculating visualization for viability threshold {}'.format(x))
        data['fit'] = (data['y'] > x).astype(float)
        neutral_p = data['fit'].mean()
        neutral_odds = neutral_p / (1-neutral_p)
        selection_odds = neutral_odds * odds_ratio
        selection_p = selection_odds / (1 + selection_odds)
        
        print('\tTime spent in functional sequences in neutrality: {}'.format(neutral_p))
        print('\tTime spent in functional sequences under selection: {}'.format(selection_p))
        space = SequenceSpace(X=data.index.values, y=data.fit.values)
        rw = WMWalk(space)
        
        print('\tCalculating coordinates')
        rw.calc_visualization(mean_function=selection_p)
        
        fpath = join(VIZ_DIR, 'threshold_landscape.t{}'.format(x))
        print('\tSaving results at {}'.format(fpath))
        rw.write_tables(fpath, write_edges=False)
