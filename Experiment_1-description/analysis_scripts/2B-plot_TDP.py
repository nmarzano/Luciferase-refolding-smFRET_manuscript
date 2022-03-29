import pandas as pd
from scipy.signal import find_peaks
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

input_folder = 'Experiment_X-description/python_results'
output_folder = f'{input_folder}/TDP_plots'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#### the first item in the tuple will be the name that goes into the graph legend ######## CAN COPY FROM 2A-INITIAL_CLEANUP_tdp
data_paths = {
    "treatment_label_1":"data_directory_1",
    "treatment_label_2":"data_directory_1",
    "treatment_label_3":"data_directory_1",
    "treatment_label_4":"data_directory_1",
    "treatment_label_5":'data_directory_1'
}

############### import cleaned data
filename = f'{input_folder}/TDP_cleaned.csv'
TDP = pd.read_csv(filename, header="infer")


def plot(treatment):
    plt.rcParams['svg.fonttype'] = 'none'
    plot1 = plt.figure(figsize = (6, 6))
    plot1 = sns.JointGrid(data = treatment, x = treatment["FRET_before"], y = treatment["FRET_after"], xlim = (0,1.2), ylim = (0, 1.2))
    plot1.plot_joint(sns.kdeplot, cmap="PuBuGn", shade=bool, cbar=False, cbar_kws={'format': '%.0f%%', 'ticks': [0, 100]}, thresh=0.07, gridsize = 200)
    plot1.ax_joint.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plot1.ax_joint.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plot1.plot_marginals(sns.histplot, kde=True, bins=10)
    plot1.ax_joint.set_facecolor("#404040") # to change the colour of the background
    sns.set(style = 'whitegrid', font_scale = 1.5)
    return plot1


for data_name, data_path in data_paths.items():
    treatment = TDP[TDP["treatment_name"] == data_name]
    plot(treatment).savefig(f"{output_folder}/TDP_plot_{data_name}.svg", dpi = 600)

