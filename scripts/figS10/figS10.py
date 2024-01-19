import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from os.path import join
from scipy.stats import pearsonr
from scripts.settings import FIGURES_DIR, DATA_DIR


def plot_uncertainty(df, x, shade=False, alpha=0.12, color='black', q=0.05,
                     zorder=0, add_expectation=False):
    lower, upper = df['lower_{}'.format(q)], df['upper_{}'.format(q)]
    
    if add_expectation:
        axes.plot(df[x], df['y'], lw=1, c=color, #label='Theoretical',
              zorder=zorder, linestyle='--')
        axes.plot(df[x], df['y_mean'], lw=1, c=color, #label='Predicted',
                  zorder=zorder)
    else:
        axes.plot(df[x], df['y'], lw=1, c=color,
                  zorder=zorder, linestyle='--')
    
    if shade:
        axes.fill_between(df[x], y1=lower, y2=upper, alpha=alpha,
                          color=color, lw=0, zorder=zorder)
    else:
        for v in [lower, upper]:
            axes.plot(df[x], v, linestyle='--', lw=1, color=color,
                      alpha=0.5, zorder=zorder)


if __name__ == '__main__':
    print('Reading data tables for plotting')
    data = pd.read_csv(join(DATA_DIR, 'data.calibrated_predictions.csv'))
    phi_ypred = pd.read_csv(join(DATA_DIR, 'phi_to_y_pred.csv'))
    cv_data = pd.read_csv(join(DATA_DIR, 'cv_table.csv'))
    y_ypred = pd.read_csv(join(DATA_DIR, 'yhat_to_y_pred.csv'))
    full_data = pd.read_parquet(join(DATA_DIR, 'calibrated_predictions.pq'))
    
    fig, subplots = plt.subplots(1, 4, figsize=(12, 3))

    print('\tPanel A: replicates comparison')
    axes = subplots[0]
    axes.plot((-1, 7), (-1, 7), lw=0.3, c='grey', alpha=0.5)
    sns.scatterplot(x='norm_mic1', y='norm_mic2', data=data, linewidth=0, alpha=0.2, ax=axes,
                    color='purple')
    axes.set(xlabel='Normalized MIC score (rep 1)', ylabel='Normalized MIC score (rep 2)',
             xlim=(-1, 7), ylim=(-1, 7), title='Experimental replicates')
    r = pearsonr(data['norm_mic1'], data['norm_mic2'])[0]
    axes.text(0.95, 0.05, r'r=' + '{:.2f}'.format(r), transform=axes.transAxes,
              ha='right', va='bottom')
    axes.text(-0.25, 1.05, r'$\bf{A}$', transform=axes.transAxes, fontsize=14)

    print('\tPanel B: training set predictions')
    calibration = np.mean((data['y_pred_lower'] <= data['y']) & (data['y'] <= data['y_pred_upper']))
    print('\t\tCalibration = {:.2f}%'.format(calibration * 100))
    axes = subplots[1]
    plot_uncertainty(phi_ypred, x='phi', shade=True, alpha=0.1, color='black',
                     q=0.05,add_expectation=True)
    sns.scatterplot(x='phi_pred', y='viability_score', data=data, linewidth=0,
                    alpha=0.2, ax=axes, color='purple')
    axes.set(xlabel=r'Predicted $\Delta\Delta G_{binding}$ [R.e.u.]', ylabel=r'Obs. normalized MIC score',
             xlim=(phi_ypred['phi'].values[0], phi_ypred['phi'].values[-1]),
             ylim=(-1, 7), title='Model fitting')
    axes.text(-0.25, 1.05, r'$\bf{B}$', transform=axes.transAxes, fontsize=14)

    print('\tPanel C: cross-validation predictions')
    axes = subplots[2]
    axes.plot((-1, 7), (-1, 7), lw=0.3, c='grey', alpha=0.5)
    plot_uncertainty(y_ypred, x='y', shade=True, alpha=0.1, color='black',
                     q=0.05,add_expectation=False)
    sns.scatterplot(x='pred', y='obs', data=cv_data, linewidth=0,
                    alpha=0.2, ax=axes, color='purple')
    axes.set(xlabel='Pred. normalized MIC score', ylabel='Obs. normalized MIC score',
             xlim=(-1, 7), ylim=(-1, 7), title='Leave-one-out cross-validation')
    
    r = pearsonr(cv_data['pred'], cv_data['obs'])[0]
    logL = cv_data['cv_ll'].sum()
    calibration = np.mean((cv_data['pred_lower'] <= cv_data['obs']) & (cv_data['obs'] <= cv_data['pred_upper']))
    print('\t\tr={:.2f}; logL={:.2f}; calibration={:.2f}%'.format(r, logL, calibration * 100))
    axes.text(0.95, 0.05, r'Pearson r=' + '{:.2f}'.format(r, logL), transform=axes.transAxes,
              ha='right', va='bottom')
    axes.text(-0.25, 1.05, r'$\bf{C}$', transform=axes.transAxes, fontsize=14)

    # Show predictions in the full dataset and calibration range
    print('\tPanel D: calibration range')
    axes = subplots[3]
    sns.histplot(full_data['phi'], color='purple', ax=axes, bins=30, stat='probability')
    ymax = axes.get_ylim()[1]
    axes.fill_between((data.phi_pred.min(), data.phi_pred.max()),
                      y1=(0, 0), y2=(ymax, ymax),
                      color='purple', alpha=0.1, lw=0)
    axes.set(xlabel=r'Predicted $\Delta\Delta G_{binding}$ [R.e.u.]',
             ylabel='Fraction Im/E pairs', ylim=(0, ymax), xlim=(-10, 25),
             title='Full dataset')
    axes.text(-0.25, 1.05, r'$\bf{D}$', transform=axes.transAxes, fontsize=14)
    
    fig.tight_layout(w_pad=2.5)
    fig.savefig(join(FIGURES_DIR, 'figS10.png'), dpi=300)
    fig.savefig(join(FIGURES_DIR, 'figS10.pdf'), dpi=300, format='pdf')

    