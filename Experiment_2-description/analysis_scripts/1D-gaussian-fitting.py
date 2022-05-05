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
from lmfit import Model, Parameter, report_fit
import math 
import numpy as np


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





##########
########## Theoretical demonstration of skew
##########

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

x_values = np.linspace(0, 1, 100)
for mu, sig in [(.8, .1)]:
    plt.plot(x_values, gaussian(x_values, mu, sig))

plt.show()


pair_distance =  np.linspace(0, 1, 100)
forster_radius = 0.51
FRET = 1/(1 + np.power((pair_distance/forster_radius), 6.))

sns.scatterplot(x = pair_distance, y = FRET)
plt.show()

Forster_values = np.linspace(0.5, 1, 6)

for mu, sig in [(.404, .05)]:
    plt.plot(FRET, gaussian(pair_distance, mu, sig))
    plt.xlabel('FRET')

FRET_peaks = np.linspace(0, 1, 9)
FRET_peaks_df = pd.DataFrame(FRET_peaks)
FRET_peaks_df.columns = ['FRET_peak']



for value in FRET_peaks:
    FRET_2 = 1/(1 + np.power((pair_distance/forster_radius), 6.))
    for mu, sig in [(value, .1)]:
        plt.plot(FRET_2, gaussian(pair_distance, mu, sig))
        plt.xlabel('FRET')
plt.show()


plot1 = plt.figure()
plt.rcParams['svg.fonttype'] = 'none'

for mu, sig in [(.8, .05)]:
    plt.plot(x_values, gaussian(x_values, mu, sig))

for mu, sig in [(.404, .05)]:
    plt.plot(FRET, gaussian(pair_distance, mu, sig))


sns.scatterplot(x = FRET, y = pair_distance)
# plot1.legend(['Pair-distance', 'FRET transformed pair-distance', 'Pair-distance/FRET relationship'])
plot1.savefig(f'{output_folder}/comparison-with-fret-on-xaxis.svg', dpi = 600)
plt.show


pair_distance =  np.linspace(0, 1, 100)
forster_radius = 0.51
FRET_peaks = np.linspace(0, 1, 9)
FRET_peaks_df = pd.DataFrame(FRET_peaks)
FRET_peaks_df.columns = ['FRET_peak']

def plot_skew(df, sigma, pair_distance = pair_distance, forster_radius = 0.51):
    y_data = []
    for value, dfs in df.groupby('FRET_peak'):
        FRET_2 = 1/(1 + np.power((pair_distance/forster_radius), 6.))
        for mu, sig in [(value, sigma)]:
            y = gaussian(pair_distance, mu, sig)
        y_df = pd.DataFrame(y)
        y_df['mu'] = mu
        y_df['xdata'] = FRET_2
        y_data.append(y_df)
        concat_data = pd.concat(y_data)
        concat_data.columns = ['ydata', 'mu', 'xdata']
    return concat_data


skew_plot = plot_skew(FRET_peaks_df, 0.05)

plt.figure()
plt.rcParams['svg.fonttype'] = 'none'
sns.set_style('darkgrid') ## alternatively, 'ticks'
sns.lineplot(data = skew_plot, y = 'ydata', x = 'xdata', hue = 'mu', palette = 'mako', legend = 'full')
plt.legend(loc = "upper left", bbox_to_anchor = (1,1), ncol =1)
plt.ylabel('Probability')
plt.xlabel('FRET')
plt.savefig(f'{output_folder}/Different-mu.svg', dpi = 600)
plt.show()



####################
#################### Code to fit experimental data based on the inherent skew arising due to non-linear relationship
#################### between dye-pair distances and FRET
####################



def gaussian_amp(x, mu, sig, height):
    return (height)* np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def FRET_distance_gen(x, mu, sig, height, forster = 0.51): ### mu is dye pair distance
    return gaussian_amp((((1/x)-1)**(1/6))*forster, mu, sig, height)
   

def FRET_gaussian_fitting(df, treatment, mu_1, sig_1, height_1, mu_2, sig_2, height_2, mu_3, sig_3, height_3):
    combined_mod = Model(FRET_distance_gen, independent_vars=['x'], prefix='m1_') + Model(FRET_distance_gen, independent_vars=['x'], prefix='m2_')+ Model(FRET_distance_gen, independent_vars=['x'], prefix='m3_')
    pars = combined_mod.make_params()
    pars.add('m1_mu', value = mu_1, min=0, max=1) ## 0.833
    pars.add('m2_mu', value = mu_2, min=0, max=1)  ## 0.4
    pars.add('m3_mu', value = mu_3, min=0, max=1) ## 0.466
    pars.add('m1_sig', value = sig_1, min=0, max=.5) ## 0.0285
    pars.add('m2_sig', value = sig_2, min=0, max=.5) ## 0.03
    pars.add('m3_sig', value = sig_3, min=0, max=.5)  ## 0.05
    pars.add('m1_height', value = height_1, min=0, max=1) ## 0.55
    pars.add('m2_height', value = height_2, min=0, max=1) ## 0.6
    pars.add('m3_height', value = height_3, min=0, max=1) ## 1

    filt_med = df[df['treatment_name'] == treatment]
    bins = np.arange(0, 1, 0.025) 
    inds = np.digitize(filt_med['FRET'].astype(float), bins)
    xdata, ydata = np.unique(inds, return_counts=True)
    ydata = ydata[1:-1] #### trim off outside range bins at the end
    xdata = [np.mean(bins[x : x + 2]) for x in range(len(bins)- 1)]  ##### convert bin edges to bin centres, therefore end up with one less bin
    sns.lineplot(xdata, ydata)

    output = combined_mod.fit(ydata/ydata.max(), pars, x=np.array(xdata))
    output.plot()
    fig = sns.set_style('darkgrid')
    fig = output.plot(data_kws={'markersize': 3})


    mod_1 = Model(FRET_distance_gen, independent_vars=['x'], prefix='m1_')
    mod_2 = Model(FRET_distance_gen, independent_vars=['x'], prefix='m2_')
    mod_3 = Model(FRET_distance_gen, independent_vars=['x'], prefix='m3_')
    test_comb = mod_1 + mod_2 + mod_3

    params_1 = mod_1.make_params(mu = mu_1, sig = sig_1, height = height_1, min = 0, max = 1)
    params_2 = mod_2.make_params(mu = mu_2, sig = sig_2, height = height_2, min = 0, max = 1)
    params_3 = mod_3.make_params(mu = mu_3, sig = sig_3, height = height_3, min = 0, max = 1)
    params = params_1.update(params_2)
    params = params.update(params_3)

    paramaters = {name:output.params[name].value for name in output.params.keys()}


    fit1 = mod_1.eval(x = FRET, mu = paramaters['m1_mu'], height = abs(paramaters['m1_height']), sig = paramaters['m1_sig'])
    fit2 = mod_2.eval(x = FRET, mu = paramaters['m2_mu'], height = abs(paramaters['m2_height']), sig = paramaters['m2_sig'])
    fit3 = mod_3.eval(x = FRET, mu = paramaters['m3_mu'], height = abs(paramaters['m3_height']), sig = paramaters['m3_sig'])
    sns.lineplot(FRET, fit1)
    sns.lineplot(FRET, fit2)
    sns.lineplot(FRET, fit3)

    fig.savefig(f'{output_folder}/{treatment}_gaussfit_no-skew.svg', dpi = 600)

    # Calculate area under the curve for each gaussian
    aoc_m1 = paramaters['m1_height']
    aoc_m2 = paramaters['m2_height']
    aoc_m3 = paramaters['m3_height']

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
    return proportion_df, paramaters, output

proportion, parameters, output = FRET_gaussian_fitting(compiled_df, 'medium', 0.833, 0.0285, 0.55, 0.40, 0.03, 0.6, 0.466, 0.05, 1)