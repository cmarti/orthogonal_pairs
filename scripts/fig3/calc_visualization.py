import pandas as pd

from os.path import join
from scripts.settings import DATA_DIR, VIZ_DIR
from gpmap.src.space import SequenceSpace
from gpmap.src.randwalk import WMWalk


if __name__ == '__main__':
    fpath = join(DATA_DIR, 'calibrated_predictions.pq')
    print('Loading calibrated predictions from {}'.format(fpath))
    data = pd.read_parquet(fpath)
    
    print('Creating sequence space object and building rate matrix')
    space = SequenceSpace(X=data.index.values, y=data.y.values)
    rw = WMWalk(space)
    
    print('Calculating visualization for y_mean=2')
    rw.calc_visualization(mean_function=2)
    
    fpath = join(VIZ_DIR, 'visualization')
    print('Writing output to {}.*'.format(fpath))
    rw.write_tables(fpath, write_edges=True)
