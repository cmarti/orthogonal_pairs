import pandas as pd
import numpy as np

from os.path import join
from tqdm import tqdm
from scipy.optimize import minimize
from scipy.stats import norm

from scripts.settings import ENERGIES_FPATH, DATA_DIR, CHANGE_POS, IM9
from scripts.utils import load_extra_exp_data
from scipy.stats.stats import spearmanr, pearsonr


class CalibrationModel():
    def __init__(self, train_gamma=True):
        params = ['alpha', 'beta', 'log_sigma']
        if train_gamma: 
            params.append('gamma')
        else:
            self.gamma = 0
            
        self.train_gamma = train_gamma
        self.n_params = len(params)
        self.param_names = params
    
    def init_param(self, param_name, value, params, extra_params):
        if value is None:
            params.append(param_name)
        else:
            extra_params[param_name] = value
    
    def calc_phi(self, data, params=None):
        pars = self.get_params(params)
        phi1 = data['x1']
        phi2 = data['x2'] + pars['gamma']
        phi = np.vstack([phi1, phi2])
        return(phi)
    
    def phi_to_yhat(self, phi, params=None):
        pars = self.get_params(params)
        w = np.exp(pars['alpha'] + pars['beta'] * phi.min(0))
        yhat = np.log(1 + w)
        return(yhat)
    
    def get_params(self, params=None):
        if params is None:
            params = self.params
        pars = {k: v for k, v in zip(self.param_names, params)}
        if not self.train_gamma:
            pars['gamma'] = 0
        return(pars)
    
    def likelihood(self, y, yhat, sigma, y_max=6):
        p = norm.pdf(y, loc=yhat, scale=sigma)

        # Change for upper bound
        upper_sat = y == y_max
        p[upper_sat] = 1 - norm.cdf(x=6, loc=yhat[upper_sat], scale=sigma)

        # Update for lower bound
        lower_sat = y == 0
        p[lower_sat] = norm.cdf(x=0, loc=yhat[lower_sat], scale=sigma)
        return(p)
    
    def fit(self, data, niter=10000):
        self.data = data
        x0 = np.zeros(self.n_params)
        self.res = minimize(self.neg_log_likelihood, x0=x0, options={'maxiter': niter})
        self.ll = -self.res.fun
        self.params = self.res.x
        return(self.params)
    
    def x_to_yhat(self, parameters):
        phi = self.calc_phi(self.data, parameters)
        yhat = self.phi_to_yhat(phi, parameters)
        return(yhat)
    
    def neg_log_likelihood(self, parameters):
        yhat = self.x_to_yhat(parameters)
        sigma = np.exp(self.get_params(parameters)['log_sigma'])
        
        with np.errstate(divide='ignore'):
            neg_ll = -np.log(self.likelihood(self.data['y'], yhat, sigma)).sum()
            if np.isinf(neg_ll):
                neg_ll = 1e16
        return(neg_ll)
    
    def predict(self, data, params=None):
        if params is None:
            params = self.params
            
        phi = self.calc_phi(data, params)
        y = self.phi_to_yhat(phi, params)
        
        if 'y_max' not in data:
            y_max = np.ones(y.shape[0]) * 6
        else:
            y_max = data['y_max']
            
        saturated = y > y_max
        y[saturated] = y_max[saturated]
        return(phi, y)
    
    def yhat_expectation(self, yhat, y_min=0, y_max=6):
        sigma = np.exp(self.get_params()['log_sigma'])
        alpha = (y_min - yhat) / sigma
        beta = (y_max - yhat) / sigma
        # Expectation conditioning on being within the [0, 6] interval
        exp_yhat = yhat - sigma * (norm.pdf(beta) - norm.pdf(alpha)) / (norm.cdf(beta) - norm.cdf(alpha))
        
        p1 = norm.cdf(y_min, loc=yhat, scale=sigma)
        p3 = 1-norm.cdf(y_max, loc=yhat, scale=sigma)
        p2 = 1 - p1 - p3
        return(y_min * p1 + exp_yhat * p2 + y_max * p3)
    
    def y_pred_ci(self, yhat, q, y_min=0, y_max=6):
        sigma = np.exp(self.get_params()['log_sigma'])
        lower = norm.ppf(q, loc=yhat, scale=sigma)
        upper = norm.ppf(1-q, loc=yhat, scale=sigma)
        lower[lower < y_min] = y_min
        lower[lower > y_max] = y_max
        upper[upper < y_min] = y_min
        upper[upper > y_max] = y_max
        return(lower, upper)
    
    def write_params(self, fpath):
        pars_df = pd.DataFrame({'param': pars.keys(),
                                'estimate': pars.values()})
        pars_df.to_csv(fpath)
    
    def read_params(self, fpath):
        pars_df = pd.read_csv(fpath).set_index('param')
        self.params = pars_df.loc[self.param_names, 'estimate'].values
        
    
def get_predictions(data, model):
    df = data.copy()
    model_data = {'x1': df.col2_binding.values, 
                  'x2': df.col9_binding.values,
                  'y': df.y.values}
    params = model.fit(model_data)
    phi, y = model.predict(model_data, params=params)
    if len(phi.shape) == 1:
        df['phi_pred'] = phi
    else:
        df['phi_pred'] = phi.min(0)
    
    df['y_pred'] = y    
    df['y_pred_lower'], df['y_pred_upper'] = model.y_pred_ci(y, q=0.05)
    return(df, params)


def calc_cv_table(data, model):
    flags = np.full(data.shape[0], False)
    cv_data = []
    for i in tqdm(np.arange(data.shape[0])):
        sel = flags.copy()
        sel[i] = True

        training, validation = data.loc[~sel, :], data.loc[sel, :]
        train_data = {'x1': training.col2_binding.values,
                      'x2': training.col9_binding.values,
                      'y': training.y.values}
        model.fit(train_data)
        sigma = np.exp(model.get_params()['log_sigma'])
        
        val_data = {'x1': validation.col2_binding.values,
                    'x2': validation.col9_binding.values,
                    'y': validation.y.values}
        val_pred = model.predict(val_data)[1]
        lower, upper = model.y_pred_ci(val_pred, q=0.05)
        cv_ll = np.log(model.likelihood(val_data['y'], val_pred, sigma))[0]
        
        if hasattr(val_pred, 'values'):
            val_pred = val_pred.values

        cv_data.append({'ll': model.ll, 'cv_ll': cv_ll,
                        'pred': min(val_pred[0], 6),
                        'pred_upper': upper[0],
                        'pred_lower': lower[0],
                        'obs': val_data['y'][0]})
    cv_data = pd.DataFrame(cv_data)
    return(cv_data)


if __name__ == '__main__':
    data = pd.read_csv(join(DATA_DIR, 'merged_exp_data.csv'))
    data = data.loc[data['backbone'] == 'E2/Im2', :].dropna()
    data['norm_mic1'] = data['E_naive_mic'] - data['mic1']
    data['norm_mic2'] = data['E_naive_mic'] - data['mic2']
    
    # Fit with the experimental data
    print('Model fitting')
    model = CalibrationModel(train_gamma=True)
    d, params = get_predictions(data, model)
    pars = model.get_params(model.params)
    for param, value in pars.items():
        print('\t{} = {:.2f}'.format(param, value))
    calibration = np.mean((d['y_pred_lower'] <= d['y']) & (d['y'] <= d['y_pred_upper']))
    print('\tCalibration = {:.2f}'.format(calibration))
    print('\tlogL = {}'.format(model.ll))
    model.write_params(join(DATA_DIR, 'calibration.params.csv'))
    d.to_csv(join(DATA_DIR, 'data.calibrated_predictions.csv'))
    
    
    # Making predictions for the whole range of phi and y
    print('Making predictions for the whole range of phi')
    phi = np.linspace(d.phi_pred.min()-1, d.phi_pred.max()+1, 100)
    y = model.phi_to_yhat(np.array([phi]), model.params)
    pred = {'phi': phi, 'y': y, 'y_mean': model.yhat_expectation(y)}
    for q in np.linspace(0.01, 0.50, 50):
        lower, upper = model.y_pred_ci(y, q=q)
        pred['lower_{}'.format(q)] = lower
        pred['upper_{}'.format(q)] = upper
    pd.DataFrame(pred).to_csv(join(DATA_DIR, 'phi_to_y_pred.csv'))
    
    print('Making predictions for the whole range of y')
    y = np.linspace(0, 6, 100)
    pred = {'y': y, 'y_mean': model.yhat_expectation(y)}
    for q in np.linspace(0.01, 0.50, 50):
        lower, upper = model.y_pred_ci(y, q=q)
        pred['lower_{}'.format(q)] = lower
        pred['upper_{}'.format(q)] = upper
    pd.DataFrame(pred).to_csv(join(DATA_DIR, 'yhat_to_y_pred.csv'))
    
    # Predict in the full landscape
    print('Making predictions in the complete landscape')
    energies = pd.read_parquet(ENERGIES_FPATH)
    data_pred = {'x1': energies.col2_binding.values,
                 'x2': energies.col9_binding.values,
                 'x3': np.zeros(energies.shape[0]),
                 'x4': np.zeros(energies.shape[0])}
    phi_pred, ypred = model.predict(data_pred)
    energies['phi'] = phi_pred.min(0)
    energies['phi1'] = phi_pred[0]
    energies['phi2'] = phi_pred[1]
    energies['y'] = ypred
    energies['y1'] = model.phi_to_yhat(np.array([phi_pred[0]]))
    energies['y2'] = model.phi_to_yhat(np.array([phi_pred[1]]))
    energies.to_parquet(join(DATA_DIR, 'calibrated_predictions.pq'))
    
    # Predict in the extra test data
    print('Making predictions in the test data')
    extra = load_extra_exp_data(energies)
    extra_data = {'x1': extra.col2_binding.values,
                  'x2': extra.col9_binding.values,
                  'x3': np.zeros(extra.shape[0]),
                  'x4': np.zeros(extra.shape[0])}
    extra['y'] = model.predict(extra_data)[1]
    extra.to_csv(join(DATA_DIR, 'test_pred.csv'))
    
    # Cross-validated correlation
    print('Running cross-validation analysis')
    cv = calc_cv_table(data, model)
    cv.to_csv(join(DATA_DIR, 'cv_table.csv'))
