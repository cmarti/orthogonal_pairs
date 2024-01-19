import pandas as pd
import numpy as np

from os.path import join

from scripts.settings import DATA_DIR
from scipy.stats.stats import spearmanr, pearsonr


if __name__ == '__main__':
    training = pd.read_csv(join(DATA_DIR, 'data.calibrated_predictions.csv'), index_col=0)
    cv = pd.read_csv(join(DATA_DIR, 'cv_table.csv'))
    pars = pd.read_csv(join(DATA_DIR, 'calibration.params.csv'), index_col=0).set_index('param')['estimate'].to_dict()
    pred = pd.read_parquet(join(DATA_DIR, 'calibrated_predictions.pq'))
    
    print('Total number of observations: {}'.format(cv.shape[0]))
    n = np.sum(np.logical_and(cv['obs'] > 5, cv['pred'] > 5))
    print('Genotypes with predicted and observed norm MIC > 5: {}'.format(n))
    n = np.sum(np.logical_and(cv['obs'] < 1, cv['pred'] < 1))
    print('Genotypes with predicted and observed norm MIC > 1: {}'.format(n))
    n = np.sum((cv['obs'] > 3) == (cv['pred'] > 3))
    print('Genotypes correctly predicted te be larger or smaller than 3: {} ({:.2f} %)'.format(n, n / cv.shape[0] * 100))
    print('Maximum Likelihood sigma = {:.2f}'.format(np.exp(pars['log_sigma'])))

    print('\tTraining set metrics')
    residuals = training['y'] - training['y_pred']
    calibration = np.mean((training['y_pred_lower'] <= training['y']) & (training['y'] <= training['y_pred_upper']))
    r, p = pearsonr(training['y_pred'], training['y'])
    s = spearmanr(training['y_pred'], training['y'])[0]
    rmse = np.sqrt(np.mean(residuals ** 2))
    mae = np.mean(np.abs(residuals))
    
    print('\tPearson r = {:.2f} ; p = {}'.format(r, p))
    print('\tSpearman r = {:.2f}'.format(s))
    print('\tMean of residuals = {:.2f}'.format(np.mean(residuals)))
    print('\tRMSE = {:.2f}'.format(rmse))
    print('\tMAE = {:.2f}'.format(mae))
    print('\tCalibration = {:.2f}'.format(calibration))
    
    print('\nCross-validation metrics')
    calibration = np.mean((cv['pred_lower'] <= cv['obs']) & (cv['obs'] <= cv['pred_upper']))
    r, p = pearsonr(cv['pred'], cv['obs'])
    s = spearmanr(cv['pred'], cv['obs'])[0]
    rmse = np.sqrt(np.mean((cv['pred'] - cv['obs']) ** 2))
    mae = np.mean(np.abs(cv['pred'] - cv['obs']))
    logL = cv['cv_ll'].sum()
    
    print('\tPearson r = {:.2f} ; p = {}'.format(r, p))
    print('\tSpearman r = {:.2f}'.format(s))
    print('\tMean of residuals = {:.2f}'.format(np.mean(cv['obs'] - cv['pred'])))
    print('\tRMSE = {:.2f}'.format(rmse))
    print('\tMAE = {:.2f}'.format(mae))
    print('\tlogL = {:.2f}'.format(logL))
    print('\tCalibration = {:.2f}'.format(calibration))

    print('\nSimpsons paradox evaluation')
    x, y = cv['pred'].values, cv['obs'].values
    
    idx = cv['obs'].values > 3
    r, p = pearsonr(x[idx], y[idx])
    s = spearmanr(x[idx], y[idx])[0]
    print('\tNorm MIC > 3: Pearson r = {:.2f}, p = {}; Spearman r = {:.2f}'.format(r, p, s))
    
    idx = cv['obs'].values <= 3
    r, p = pearsonr(x[idx], y[idx])
    s = spearmanr(x[idx], y[idx])[0]
    print('\tNorm MIC <= 3: Pearson r = {:.2f}, p = {}; Spearman r = {:.2f}'.format(r, p, s))

    print('Complete dataset')
    print('\tMedian Normalized MIC score: {:.2f}'.format(pred['y'].median()))