from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
from matplotlib import rc
import matplotlib as mpl
from matplotlib import cm 
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import glob as glob

from Utilities.Data_analysis import file_reader

output_folder = 'Experiment_1-PPR-with-1uM-atpH/python_results'
exposure = 0.200  # in seconds
frame_rate = 1/exposure

data_paths = {
    'treatment':('treatment_name', 'Experiment_1-description/raw_data/treatment_to_plot'),
}

def remove_outliers(compiled, plot_type, data_type = "raw"):
    """[removes outliers from dataframe]

    Args:
        compiled ([dataframe]): [raw dataframe containing outliers to be removed]
        plot_type ([str]): [string can either be 'hist' for histogram data or 'TDP' for TDP data]
        data_type (str, optional): [removes either raw FRET values or 'idealized' FRET values]. Defaults to "raw".

    Returns:
        [dataframe]: [returns cleaned data without outliers]
    """
    if plot_type == 'hist':
        if data_type == "raw":
            return compiled[(compiled['FRET'] > -0.5) & (compiled['FRET'] < 1.5)].copy()
        if data_type == "idealized":
            return compiled[(compiled[4] > -0.5) & (compiled[4] < 1.5)].copy()
    elif plot_type == 'TDP':
        outliers = compiled[(compiled["FRET before transition"] < -0.5)|(compiled["FRET before transition"] > 1.5)|(compiled["FRET after transition"] < -0.5) | (compiled["FRET after transition"] > 1.5)].index
        compiled.drop(outliers, inplace = True)
        return compiled
    else:
        print('invalid plot type, please set plot_type as "hist" or "TDP" - you idiot')

def plot_heatmap(df, gridsize, bins_hex):
    for treatment, dfs in df.groupby('treatment_name'):
        plt.rcParams['svg.fonttype'] = 'none'
        sns.set(style="whitegrid")
        g = sns.JointGrid(data = dfs, x='time', y='FRET', xlim = (0,100), ylim = (0, 1))
        g.plot_joint(plt.hexbin, gridsize=(gridsize, gridsize), cmap='ocean_r', mincnt=0, bins=bins_hex)
        g.plot_marginals(sns.histplot, kde=True, bins=20)
        identifier = treatment
        plt.savefig(f'{output_folder}/Heatmap_{identifier}.svg', dpi = 600)
        plt.show()
    return

compiled_df = []
for data_name, (label, data_path) in data_paths.items():
    imported_data = file_reader(data_path, 'heatmap', frame_rate)
    cleaned_raw = remove_outliers(imported_data, 'hist')    #### add "idealized" after imported_data to get idealized histograms
    cleaned_raw["treatment_name"] = data_name
    compiled_df.append(cleaned_raw)
compiled_df = pd.concat(compiled_df)   #### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed
compiled_df.columns = ["frames", "donor", "acceptor", "FRET", "idealized FRET", 'molecule_number', 'time', "treatment_name"]
plot_heatmap(compiled_df, 100, 80)

