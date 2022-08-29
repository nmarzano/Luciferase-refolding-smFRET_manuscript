import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pandas as pd
import matplotlib
import os
import math
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize, signal
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.integrate as integrate
from lmfit import models

input_folder = 'Experiment_1-description/python_results'
output_folder = f"{input_folder}/GaussianFits"  ### modify for each experiment
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

filename = f'{input_folder}/Cleaned_FRET_histogram_data.csv'
compiled_df = pd.read_csv(filename, header="infer")

############ 
############
############
############   Hsp70 low-FRET peak not constrained here
def fit_gauss_dif_constrained_nativespont(df, treatment, mu_1, sigma_1, amplitude_1, gamma_1, mu_2, sigma_2, amplitude_2, mu_3, sigma_3, amplitude_3, gamma_3):
    filt_df = df[df['treatment_name'] == treatment]
    bins = np.arange(-0.21, 1.1, 0.025) 
    inds = np.digitize(filt_df['FRET'].astype(float), bins)
    xdata, ydata = np.unique(inds, return_counts=True)
    ydata = ydata[1:-1] #### trim off outside range bins at the end
    xdata = [np.mean(bins[x : x + 2]) for x in range(len(bins)- 1)]  ##### convert bin edges to bin centres, therefore end up with one less bin
    sns.lineplot(xdata, ydata)

    model_1 = models.SkewedGaussianModel(prefix='m1_')
    model_2 = models.GaussianModel(prefix='m2_')
    model_3 = models.SkewedGaussianModel(prefix='m3_')
    model = model_1 + model_2 + model_3 
   
    model_1.set_param_hint('m1_center', vary=False)

    model_2.set_param_hint('m2_sigma', vary=False)

    model_2.set_param_hint('m2_center', vary=False)
    model_3.set_param_hint('m3_gamma', vary=False)
    model_3.set_param_hint('m3_sigma', vary=False)
    model_3.set_param_hint('m3_center', vary=False)


    params_1 = model_1.make_params(center = mu_1, sigma = sigma_1, amplitude = amplitude_1, gamma = gamma_1, min = 0)
    params_2 = model_2.make_params(center = mu_2, sigma = sigma_2, amplitude = amplitude_2, min = 0)
    params_3 = model_3.make_params(center = mu_3, sigma = sigma_3, amplitude = amplitude_3, gamma = gamma_3, min = 0)
    params = params_1.update(params_2)
    params = params.update(params_3)

    output = model.fit((ydata/np.max(ydata)), params, x=xdata)
    fig = output.plot(data_kws={'markersize': 3})

    paramaters = {name:output.params[name].value for name in output.params.keys()}
    fitx = np.arange(-0.2, 1.2, 0.025)

    fit1 = model_1.eval(x = fitx, center = paramaters['m1_center'], amplitude = abs(paramaters['m1_amplitude']), sigma = paramaters['m1_sigma'], gamma = paramaters['m1_gamma'])
    fit2 = model_2.eval(x = fitx, center = paramaters['m2_center'], amplitude = abs(paramaters['m2_amplitude']), sigma = paramaters['m2_sigma'], fwhm = paramaters['m2_fwhm'])
    fit3 = model_3.eval(x = fitx, center = paramaters['m3_center'], amplitude = abs(paramaters['m3_amplitude']), sigma = paramaters['m3_sigma'], gamma = paramaters['m3_gamma'])

    sns.lineplot(fitx, fit1)
    sns.lineplot(fitx, fit2)
    sns.lineplot(fitx, fit3)
    plt.show()

    # Calculate area under the curve for each gaussian
    aoc_m1 = paramaters['m1_amplitude']
    aoc_m2 = paramaters['m2_amplitude']
    aoc_m3 = paramaters['m3_amplitude']

    # aoc_m1 = (paramaters['m1_amplitude']*paramaters['m1_sigma'])/0.3989
    # aoc_m2 = (paramaters['m2_amplitude']*paramaters['m2_sigma'])/0.3989
    # aoc_m3 = (paramaters['m3_amplitude']*paramaters['m3_sigma'])/0.3989

    sum_aoc = aoc_m1 + aoc_m2 + aoc_m3 

    aoc_m1_percent_of_total = (aoc_m1/sum_aoc)*100
    aoc_m2_percent_of_total = (aoc_m2/sum_aoc)*100
    aoc_m3_percent_of_total = (aoc_m3/sum_aoc)*100
    list_of_gaus_proportion = [aoc_m1_percent_of_total, aoc_m2_percent_of_total, aoc_m3_percent_of_total]
    labels_of_gaus_proportion = ['m1', 'm2', 'm3']
    proportion_df = pd.DataFrame([labels_of_gaus_proportion, list_of_gaus_proportion])
    proportion_df.columns = proportion_df.iloc[0]
    proportion_df = proportion_df.drop(0)
    proportion_df['treatment'] = treatment
    proportion_df.to_csv(f'{output_folder}/gaussian_proportions_for_{treatment}.csv')
    return proportion_df

gaussian_kj_skew_con_nat = fit_gauss_dif_constrained_nativespont(compiled_df, 'KJ', 0.05, .1, 1, 10, .63, .1, .05, .95, .22, .03, -2.7)
gaussian_high_skew_con_na2 = fit_gauss_dif_constrained_nativespont(compiled_df, 'high', 0, .1, .1, 10, .63, .13, .5, .95, .2, 1, -2.7)
gaussian_medium_skew_con2 = fit_gauss_dif_constrained_nativespont(compiled_df, 'medium', 0.00, .1, .5, 10, .63, .13, 1, .95, .2, .5, -2.7)
gaussian_low_skew_con2 = fit_gauss_dif_constrained_nativespont(compiled_df, 'low', 0.02, .1, .5, 10, .63, .1, 1, .95, .2, .5, -2.7)

collated = pd.concat([gaussian_kj_skew_con_nat,gaussian_low_skew_con2, gaussian_medium_skew_con2, gaussian_high_skew_con_na2 ])
collated.to_csv(f'{output_folder}/histogram_proportions.csv', index = False)


###### - Do not change - these conditions are now set (if you want to play around just duplicate the function)

def fit_gauss_dif_constrained_allpeaks(df, treatment, mu_1, sigma_1, amplitude_1, gamma_1, mu_2, sigma_2, amplitude_2, mu_3, sigma_3, amplitude_3, gamma_3):
    """Set paramaters and fit histogram data to a 3-gaussian model. 

    Args:
        df (dataframe): dataframe containing cleaned FRET values used to plot histogram
        treatment (str): determines what treatment you want to look at within the dataset
        mu_1 (float): set the mean of the first gaussian
        sigma_1 (float): set the value of the width of the first gaussian
        amplitude_1 (float): estimate for the height of the first gaussian
        gamma_1 (float): sets the skew parameter - positive values result in skew to right and negative values result in skew to the left
        mu_2 (float): set the mean of the second gaussian
        sigma_2 (float): estimate for the width of the second gaussian
        amplitude_2 (float): estimate for the height of the second gaussian
        mu_3 (float): estimate for the mean of the third gaussian
        sigma_3 (float): set the width of the third gaussian
        amplitude_3 (float): estimate for the height of the third gaussian
        gamma_3 (float): set the skew parameter - positive values result in skew to right and negative values result in skew to the left

    Returns:
        dataframe, plots: returns the proportional area of each gausssian relative to the sum of all three gaussians. Also shows what the fit of each gaussian looks like.
    """
    filt_df = df[df['treatment_name'] == treatment]
    bins = np.arange(-0.21, 1.1, 0.025) 
    inds = np.digitize(filt_df['FRET'].astype(float), bins)
    xdata, ydata = np.unique(inds, return_counts=True)
    ydata = ydata[1:-1] #### trim off outside range bins at the end
    xdata = [np.mean(bins[x : x + 2]) for x in range(len(bins)- 1)]  ##### convert bin edges to bin centres, therefore end up with one less bin
    sns.lineplot(xdata, ydata)

    model_1 = models.SkewedGaussianModel(prefix='m1_')
    model_2 = models.GaussianModel(prefix='m2_')
    model_3 = models.SkewedGaussianModel(prefix='m3_')
    model = model_1 + model_2 + model_3 

    model_1.set_param_hint('m1_gamma', vary=False)
    model_1.set_param_hint('m1_sigma', vary=False)
    model_1.set_param_hint('m1_center', vary=False)

    model_2.set_param_hint('m2_sigma', vary=False)
    model_2.set_param_hint('m2_center', vary=False)
    model_3.set_param_hint('m3_gamma', vary=False)
    model_3.set_param_hint('m3_sigma', vary=False)


    params_1 = model_1.make_params(center = mu_1, sigma = sigma_1, amplitude = amplitude_1, gamma = gamma_1, min = 0)
    params_2 = model_2.make_params(center = mu_2, sigma = sigma_2, amplitude = amplitude_2, min = 0)
    params_3 = model_3.make_params(center = mu_3, sigma = sigma_3, amplitude = amplitude_3, gamma = gamma_3, min = 0)
    params = params_1.update(params_2)
    params = params.update(params_3)

    output = model.fit((ydata/np.max(ydata)), params, x=xdata)
    fig = sns.set_style('darkgrid')
    fig = output.plot(data_kws={'markersize': 3})

    paramaters = {name:output.params[name].value for name in output.params.keys()}
    fitx = np.arange(-0.2, 1.2, 0.025)

    fit1 = model_1.eval(x = fitx, center = paramaters['m1_center'], amplitude = abs(paramaters['m1_amplitude']), sigma = paramaters['m1_sigma'], gamma = paramaters['m1_gamma'])
    fit2 = model_2.eval(x = fitx, center = paramaters['m2_center'], amplitude = abs(paramaters['m2_amplitude']), sigma = paramaters['m2_sigma'], fwhm = paramaters['m2_fwhm'])
    fit3 = model_3.eval(x = fitx, center = paramaters['m3_center'], amplitude = abs(paramaters['m3_amplitude']), sigma = paramaters['m3_sigma'], gamma = paramaters['m3_gamma'])

    sns.lineplot(fitx, fit1)
    sns.lineplot(fitx, fit2)
    sns.lineplot(fitx, fit3)
    fig.savefig(f'{output_folder}/{treatment}_gaussfit.svg', dpi = 600)
    plt.show()

    # Calculate area under the curve for each gaussian
    aoc_m1 = paramaters['m1_amplitude']
    aoc_m2 = paramaters['m2_amplitude']
    aoc_m3 = paramaters['m3_amplitude']

    # aoc_m1 = (paramaters['m1_amplitude']*paramaters['m1_sigma'])/0.3989
    # aoc_m2 = (paramaters['m2_amplitude']*paramaters['m2_sigma'])/0.3989
    # aoc_m3 = (paramaters['m3_amplitude']*paramaters['m3_sigma'])/0.3989

    sum_aoc = aoc_m1 + aoc_m2 + aoc_m3 

    aoc_m1_percent_of_total = (aoc_m1/sum_aoc)*100
    aoc_m2_percent_of_total = (aoc_m2/sum_aoc)*100
    aoc_m3_percent_of_total = (aoc_m3/sum_aoc)*100
    list_of_gaus_proportion = [aoc_m1_percent_of_total, aoc_m2_percent_of_total, aoc_m3_percent_of_total]
    labels_of_gaus_proportion = ['m1', 'm2', 'm3']
    proportion_df = pd.DataFrame([labels_of_gaus_proportion, list_of_gaus_proportion])
    proportion_df.columns = proportion_df.iloc[0]
    proportion_df = proportion_df.drop(0)
    proportion_df['treatment'] = treatment
    proportion_df.to_csv(f'{output_folder}/gaussian_proportions_for_{treatment}.csv')
    return proportion_df

gaussian_medium_constrainall = fit_gauss_dif_constrained_allpeaks(compiled_df, 'medium', 0.02, .1, .5, 3, .63, .12, 1, .9, .22, .5, -2.7)
gaussian_col_high = fit_gauss_dif_constrained_allpeaks(compiled_df, 'high', 0.0, .20, .1, 10, .64, .1, .6, .9, .18, 1, -2.7)
gaussian_low_constrainall = fit_gauss_dif_constrained_allpeaks(compiled_df, 'low', 0.00, .1, .8, 2, .63, .11, .5, .95, .22, .4, -2.7)

collated_allconstrained = pd.concat([gaussian_kj_skew_con_nat,gaussian_low_constrainall, gaussian_medium_constrainall, gaussian_col_high ])
collated_allconstrained.to_csv(f'{output_folder}/histogram_proportions_constrained.csv', index = False)

