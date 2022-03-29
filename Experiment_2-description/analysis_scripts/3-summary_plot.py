# library
import seaborn as sns
import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import Utilities.Data_analysis as uda

FRET_thresh = 0.5

input_folder = 'Experiment_X-description/python_results'
output_folder = f"{input_folder}/Heatmaps"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def concat_meandwell(input_folder):
    filenames = glob.glob(input_folder + "/*.csv")
    mean_concat = []
    for filename in filenames:
        mean_concat.append(pd.read_csv(filename))
    test = pd.concat(mean_concat, ignore_index=True)
    test_dfs = pd.DataFrame(test)
    return test_dfs

hist_data = uda.file_reader(f'{input_folder}/Cleaned_FRET_histogram_data.csv', 'other')
transition_frequency = uda.file_reader(f'{input_folder}/Transition_frequency.csv', 'other')
mean_dwell = concat_meandwell(f'{input_folder}/Mean_dwell')

from Utilities.Data_analysis import float_generator, heatmap_prep, mean_dwell_prep

def plot_heatmap(df_to_plot, arrow_details, mean):
    fig, ax = plt.subplots(figsize=(6,7))
    plt.rcParams['svg.fonttype'] = 'none'
    sns.set(font_scale = 1.5)
    sns.heatmap(df_to_plot, cmap = "PuBuGn", annot=True, linewidths=5, linecolor= "black", cbar_kws={"label":"fraction of total time",'use_gridspec': False, 'location': 'top'}, vmin = 0, vmax = 1)
    ax.arrow(.8,.225, .3,0, width=arrow_details[0], color="darkorange", head_width=arrow_details[0]*4, head_length=0.1)
    ax.arrow(1.2,.775, -.3,0, width=arrow_details[1], color="darkorange", head_width=arrow_details[1]*4, head_length=0.1)
    ax.arrow(1.85,.75,0 ,0.1, width=arrow_details[2], color="darkorange", head_width=arrow_details[2]*2, head_length=0.1)
    ax.arrow(1.65,.95,0 ,-0.1, width=arrow_details[2], color="darkorange", head_width=arrow_details[2]*2, head_length=0.1)
    ax.arrow(0.15,.05,0 ,0.1, width=arrow_details[3], color="darkorange", head_width=arrow_details[3]*2, head_length=0.1)
    ax.arrow(0.35,.25,0,-0.1, width=arrow_details[3], color="darkorange", head_width=arrow_details[3]*2, head_length=0.1)
    plt.text(.69, .3, str(f"{float(mean[0]):.3}" + " s"))
    plt.text(.10, .3, str(f"{float(mean[3]):.3}" + " s"))
    plt.text(1.05, .725, str(f"{float(mean[1]):.3}" + " s"))
    plt.text(1.6, .725, str(f"{float(mean[2]):.3}" + " s"))
    return fig


treatment_list = mean_dwell["sample"].to_list()  ### add .unique() if you were to use the hist_data dataframe

for data_name in treatment_list:
    arrow_info = float_generator(transition_frequency, data_name, FRET_thresh)
    heatmap_info = heatmap_prep(hist_data, data_name, FRET_thresh)
    mean = mean_dwell_prep(mean_dwell, data_name, FRET_thresh)
    plot_heatmap(heatmap_info, arrow_info, mean).savefig(f"{output_folder}/Heatmap_{data_name}.svg", dpi = 600)


