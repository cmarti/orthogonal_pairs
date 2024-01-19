import pandas as pd
import numpy as np

from os.path import join

from gpmap.src.space import SequenceSpace
from gpmap.src.randwalk import WMWalk

from scripts.settings import ENERGIES_FPATH, DATA_DIR ,VIZ_DIR
from scripts.model2.figS10.calibration import CalibrationModel


if __name__ == '__main__':
    np.random.seed(0)
    energies = pd.read_parquet(ENERGIES_FPATH)    

    print('Loading calibration model')
    model = CalibrationModel()
    model.read_params(join(DATA_DIR, 'calibration.params.csv'))
    data_pred = {'x1': energies.col2_binding.values,
                 'x2': energies.col9_binding.values}
    phi = model.calc_phi(data_pred)
    odds_ratio = 3.5
    
    for sd in [0.5, 1, 1.5, 2, 2.5, 3]:
        print('Prediction with noised energies with SD={}'.format(sd))
        y_pred = model.phi_to_yhat(np.random.normal(phi, sd))
        
        print('\tCalculating visualization')
        space = SequenceSpace(X=energies.index.values, y=y_pred)
        rw = WMWalk(space)
        
        neutral_p = y_pred.mean()
        neutral_odds = neutral_p / (6-neutral_p)
        selection_odds = neutral_odds * odds_ratio
        selection_p = 6 * selection_odds / (1 + selection_odds)
        print('\tMean function at stationarity of {}'.format(selection_p))
        
        rw.calc_visualization(mean_function=selection_p)
        rw.write_tables(join(VIZ_DIR, 'noised_landscape.sd{}'.format(sd)))
