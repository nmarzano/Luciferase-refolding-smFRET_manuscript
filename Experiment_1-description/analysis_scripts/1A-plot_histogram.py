import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pandas as pd
import matplotlib
from scipy import optimize, signal
import numpy as np
from scipy import optimize
import scipy.integrate as integrate
from lmfit import models

output_folder = "Experiment_X-description/python_results"  ### modify for each experiment

#### the first item in the tuple will be the name that goes into the graph legend
data_paths = {
    "treatment_label_1":("figure_label", "data_directory_1/"),
    "treatment_label_2":("figure_label", "data_directory_2/"),
    "treatment_label_3":("figure_label", "data_directory_3/"),
    "treatment_label_4":("figure_label", "data_directory_4/"),
    "treatment_label_5":("figure_label", "data_directory_5/")
}


##########
########## Data from all data sets in the dict will be imported and concatenated into a single dataframe. Outliers wil be removed.
from Utilities.Data_analysis import remove_outliers, file_reader

compiled_df = []
for data_name, (label, data_path) in data_paths.items():
    imported_data = file_reader(data_path, 'hist')
    cleaned_raw = remove_outliers(imported_data, 'hist')    #### add "idealized" after imported_data to get idealized histograms
    cleaned_raw["treatment_name"] = data_name
    compiled_df.append(cleaned_raw)
compiled_df = pd.concat(compiled_df)   #### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed
compiled_df.columns = ["frames", "donor", "acceptor", "FRET", "idealized FRET", 'molecule_number', "treatment_name"]
compiled_df.to_csv(f'{output_folder}/Cleaned_FRET_histogram_data.csv', index = False)


###### FOR RIDGELINE PLOTS
###################################################
###################################################

labels = {data_name:label for data_name, (label, data_path) in data_paths.items()}
order = ["5","4", "3", "2", "1"] ### sets order for histogram
font = {'weight' : 'normal', 'size'   : 12 }
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

####### Sets order for ridgeline plot
data_paths_ridgeline = {
    "1":"label",
    '2':'label',
    "3":"label",
    "4":"label", 
    "5":"label", 
    "6":"label",
    "7":"label", 
    "8":"label", 
    "9":"label", 
    "10":"label", 
    "11":"label", 
    "12":"label",
}

n_colors = len(data_paths_ridgeline)
pal = sns.color_palette(palette='mako', n_colors=n_colors)
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
g = sns.FacetGrid(compiled_df, row='treatment_name', hue='treatment_name', aspect=10, height=5, palette=pal)
# then we add the densities kdeplots for each condition
g.map(sns.kdeplot, 'FRET',
      bw_adjust=1, clip_on=False,
      fill=True, alpha=1, linewidth=10)
# here we add a white line that represents the contour of each kdeplot
g.map(sns.kdeplot, 'FRET', 
      bw_adjust=1, clip_on=False, 
      color="white", linewidth=6)
# here we add a horizontal line for each plot
g.map(plt.axhline, y=0,
      lw=2, clip_on=False)
# we loop over the FacetGrid figure axes (g.axes.flat) and add the condition as text with the right color
# notice how ax.lines[-1].get_color() enables you to access the last line's color in each matplotlib.Axes
for i, ax in enumerate(g.axes.flat):
    ax.text(0, .5, data_paths_ridgeline[f'{i+1}'],
            fontweight='bold', fontsize=100,
            color=ax.lines[-1].get_color())
    ax.set_facecolor((0, 0, 0, 0))  ### removes background so KDE plots are not obstructed when stacked
# we use matplotlib.Figure.subplots_adjust() function to get the subplots to overlap
g.fig.subplots_adjust(hspace=-0.7)
# eventually we remove axes titles, yticks and spines
g.set_titles("")
g.set(yticks=([]))
g.despine(bottom=True, left=True)
plt.setp(ax.get_xticklabels(), fontsize=100, fontweight='bold')
plt.xlabel('FRET', fontweight='bold', fontsize=100)
plt.xlim(-.4, 1.2)
g.savefig(f'{output_folder}/Histogram-ridgeline.svg', dpi = 600)
plt.show()

#################
################# Code to plot regular histogram
#################
colors = {
    "Native": "black", 
    "Spontaneous":"darkorange", 
    "treatment_3":"skyblue", 
    "treatment_4":"royalblue", 
    "treatment_5":"darkviolet"
}

def plot_hist_type(df, order, kind = 'kde'):
    plot_hist, ax = plt.subplots()
    sns.set_style("whitegrid",{'figure.figsize':(8,5.5), 'grid.linestyle':'--', 'font_scale':1.5} )
    if kind == 'kde':
        sns.kdeplot(
            data = df, 
            palette = 'mako', 
            x = "FRET",
            hue="treatment_name",
            hue_order = order,
            common_norm=False, 
            fill = True, 
            linewidth = 1.5, 
            alpha = 0.25)
    if kind == 'bar':
        sns.histplot(
            data = df, 
            palette = 'mako', 
            x = "FRET",
            hue="treatment_name",
            hue_order = order,
            common_norm=False, 
            stat = 'density',
            binwidth = 0.025,
            fill = True, 
            kde = True,
            linewidth = 1.5, 
            alpha = 0.25)
    plt.xlim(0, 1, 10)
    plt.xlabel("FRET")
    plt.legend([labels[treatment] for treatment in reversed(order)],loc='upper left', fontsize = 12)
    [x.set_linewidth(2) for x in ax.spines.values()]
    [x.set_color('black') for x in ax.spines.values()]
    plot_hist.savefig(f'{output_folder}/Histogram_{kind}.svg', dpi = 600)
    plt.show()


plot_hist_type(compiled_df, order, 'kde')

