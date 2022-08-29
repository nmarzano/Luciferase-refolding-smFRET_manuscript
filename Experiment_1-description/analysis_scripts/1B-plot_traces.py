from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
from scipy.signal import savgol_filter
import numpy as np
import os

output_folder = 'Experiment_1-description/python_results/Traces'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#### the first item in the tuple will be the name that goes into the graph legend
data_paths = {
    "example1":"Figures/Figure_5/raw_data/200715_Fluc_FoldandUnfold/05_PATH_D111020_T2254.dat",
    "example2":"Figures/Figure_5/raw_data/200715_Fluc_FoldandUnfold/92_PATH_D111020_T2254.dat",
    "example3":"Figures/Figure_5/raw_data/200715_Fluc_FoldandUnfold/80_PATH_D111020_T2254.dat"
}

exposure = 0.2  ### exposure in seconds


##############
############## Code to import data and concatenate all molecule data for sequential plotting
def load_data(filename):
    trace_df = pd.DataFrame(np.loadtxt(filename))
    trace_df.columns = ["frames", "donor", "acceptor", "FRET", "idealized FRET"]
    trace_df["Time"] = trace_df["frames"]/(1/exposure)
    trace_df["smoothed_FRET"] = savgol_filter(trace_df["FRET"], 5, 2)
    return trace_df

compiled_df = []
for data_name, data_path in data_paths.items():
    imported_data = load_data(data_path)
    imported_data["treatment_name"] = data_name
    compiled_df.append(imported_data)
compiled_df = pd.concat(compiled_df)   #### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed



############
############ Code to plot FRET and/or intensity traces
def plot_intensity(treatment):
    plt.rcParams['svg.fonttype'] = 'none'
    sns.set_style("whitegrid", {'grid.linestyle':'--'})
    plot1, ax = plt.subplots(figsize = (5, 2))
    plt.xlim(0, 150, 10)
    plt.ylim(0, 4000, 0.2)
    sns.lineplot(x = treatment["Time"], y = treatment["donor"], color = 'green')
    sns.lineplot(x = treatment["Time"], y = treatment["acceptor"], color = 'purple')
    [x.set_linewidth(2) for x in ax.spines.values()]
    [x.set_color('black') for x in ax.spines.values()]
    plt.xlabel("Time (s)")
    plt.ylabel("FRET")
    plt.show()
    return plot1

def plot_FRET(treatment):
    plt.rcParams['svg.fonttype'] = 'none'
    sns.set_style("whitegrid", {'grid.linestyle':'--'})
    plot2, ax = plt.subplots(figsize = (5, 2))
    plt.xlim(0, 200, 10)
    plt.ylim(0, 1.1, 0.2)
    sns.lineplot(x = treatment["Time"], y = treatment["smoothed_FRET"], color = 'black')
    sns.lineplot(x = treatment["Time"], y = treatment["idealized FRET"], color = 'darkorange')
    [x.set_linewidth(2) for x in ax.spines.values()]
    [x.set_color('black') for x in ax.spines.values()]
    plt.xlabel("Time (s)")
    plt.ylabel("FRET")
    plt.show()
    return plot2

for data_name, data_path in data_paths.items():
    treatment = compiled_df[compiled_df["treatment_name"] == data_name]
    mol_ident = data_path.split('/')[-1]
    plot_FRET(treatment).savefig(f'{output_folder}/{data_name}_Trace_{mol_ident}.svg', dpi = 600)
    plot_intensity(treatment).savefig(f'{output_folder}/{data_name}_Trace_{mol_ident}_intensity.svg', dpi = 600)


