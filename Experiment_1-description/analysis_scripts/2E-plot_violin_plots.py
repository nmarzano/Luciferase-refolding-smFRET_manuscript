from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import Utilities.Data_analysis as uda

output_folder = "Experiment_X-description/python_results"
FRET_thresh = 0.5 #### FRET value at which to filter data above or below. 
order = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

data_paths_violin = {
    "treatment_label_1":"data_directory_1",
    "treatment_label_2":"data_directory_2",
    "treatment_label_3":"data_directory_3",
    "treatment_label_4":"data_directory_4",
    "treatment_label_5":'data_directory_5',
}
#### Directory for above should come from the Dwell_times folder in python_results

def compiled(df, data_name):
    """Will filter transitions dependent on a threshold defined above as FRET_thresh to calculate residenc time for each transition class

    Args:
        df (dataframe): dataset containing the residence times  for each treatment
        data_name (_type_): treatment name  

    Returns:
        dataframe: compiles all transition classes (with residence times) from all treatments together
    """
    violin_data_lowtolow = pd.DataFrame(df[f"< {FRET_thresh} to < {FRET_thresh}"])
    violin_data_lowtolow.columns = ["y_axis"]
    violin_data_lowtolow["transition_type"] = f"< {FRET_thresh} to < {FRET_thresh}"
    violin_data_lowtolow["treatment"] = data_name

    violin_data_lowtohigh = pd.DataFrame(df[f"< {FRET_thresh} to > {FRET_thresh}"])
    violin_data_lowtohigh.columns = ["y_axis"]
    violin_data_lowtohigh["transition_type"] = f"< {FRET_thresh} to > {FRET_thresh}"
    violin_data_lowtohigh["treatment"] = data_name

    violin_data_hightohigh = pd.DataFrame(df[f"> {FRET_thresh} to > {FRET_thresh}"])
    violin_data_hightohigh.columns = ["y_axis"]
    violin_data_hightohigh["transition_type"] = f"> {FRET_thresh} to > {FRET_thresh}"
    violin_data_hightohigh["treatment"] = data_name

    violin_data_hightolow = pd.DataFrame(df[f"> {FRET_thresh} to < {FRET_thresh}"])
    violin_data_hightolow.columns = ["y_axis"]
    violin_data_hightolow["transition_type"] = f"> {FRET_thresh} to < {FRET_thresh}"
    violin_data_hightolow["treatment"] = data_name
    return pd.concat(
        [
            violin_data_lowtolow,
            violin_data_lowtohigh,
            violin_data_hightohigh,
            violin_data_hightolow,
        ]
    )


test = []
for data_name, data_path in data_paths_violin.items():
    data = uda.file_reader(data_path, 'other')
    compiled_data = compiled(data, data_name)
    test.append(compiled_data)
final = pd.concat(test)
final["y_axis_log10"] = np.log10(final['y_axis']) ## if need to plot in log scale

############# Code to plot data
############# NEED TO MODIFY BASED ON DATASET but can COPY FROM 1A-PLOT-HISTOGRAM
colors_violin = {
    "Native": "black", 
    "Spontaneous":"darkorange", 
    "treatment_3":"skyblue", 
    "treatment_4":"royalblue", 
    "treatment_5":"darkviolet"
}


def plot_violin(data, scale = "y_axis"):
    if scale == "y_axis":
        plt.rcParams['svg.fonttype'] = 'none'
        plot1 = plt.figure(figsize = (12, 6))
        sns.set(style = "darkgrid", font_scale = 1.5)
        sns.violinplot(
            data=data, 
            x= "transition_type", 
            y = "y_axis",
            palette= colors_violin, 
            hue = "treatment", 
            log_scale = True,
            cut = 0)
        plt.ylabel("Residence time (s)")
        plt.xlabel("Transition class")
        plt.legend(loc = "upper right")
        plot1.savefig(f"{output_folder}/Violin_plot_normal.svg", dpi = 600)
        plt.show()
    if scale == "y_axis_log10":
        plt.rcParams['svg.fonttype'] = 'none'
        plot2 = plt.figure(figsize = (12, 6))
        sns.set(style = "darkgrid", font_scale = 1.5)
        sns.violinplot(
            data=data, 
            x= "transition_type", 
            y = "y_axis_log10",
            palette= colors_violin, 
            hue = "treatment", 
            log_scale = True)
        plt.ylabel("Log residence time (s)")
        plt.xlabel("Transition class")
        plt.legend(loc = "upper left", bbox_to_anchor = (1,1), ncol =1)
        plot2.savefig(f"{output_folder}/Violin_plot_log.svg", dpi = 600)
        plt.show()
    if scale == 'split':     
        f, (ax_top, ax_bottom) = plt.subplots(ncols=1, nrows=2, sharex=True, gridspec_kw={'hspace':0.05})
        sns.violinplot(x="transition_type", y="y_axis", hue="treatment",data=final, ax=ax_top, palette = 'mako')
        sns.violinplot(x="transition_type", y="y_axis", hue="treatment",data=final, ax=ax_bottom, cut = 0, palette = 'mako')
        ax_top.set_ylim(bottom=40)   # those limits are fake
        ax_bottom.set_ylim(0,40)
        sns.despine(ax=ax_bottom)
        sns.despine(ax=ax_top, bottom=True)
        ax = ax_top
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax2 = ax_bottom
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        #remove one of the legend
        ax_bottom.legend_.remove()
        f.savefig(f"{output_folder}/Violin_plot_splitaxis.svg", dpi = 600)
        plt.show()

plot_violin(final, 'split')


############# Code to generate and collate summary statistics of dwell times
mean = final.groupby(['treatment', 'transition_type']).mean()
sem =  final.groupby(['treatment', 'transition_type']).sem()
N = final.groupby(['treatment', 'transition_type']).count()
collated = pd.concat([mean,sem, N], axis = 1)
collated.drop([col for col in collated.columns.tolist() if 'y_axis_log10' in col], axis = 1, inplace = True)
collated.columns = ['mean_residence_time', 'sem', 'n']
collated.reset_index(inplace = True)
collated.to_csv(f"{output_folder}/summary.csv", index = False)


### code to plot with SEM
def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def plot_parameters_scattbar(df, x_col, x_order, y_col, hue_col, hue_order, show_bars=False, ax=False, alpha=1, err_type='std', colors=False, scat_alpha = 0):
    """Master plotting function to generate improved 'scatterbar' plots, with bar versus mean line option.

    Parameters
    ----------
    df : DataFrame
        Long-form dataset containing quantitative information to be plotted
    x_col : str
        column containing categories to be plotted along x-axis
    x_order : list[str]
        order in which x_col categories should be plotted
    y_col : str
        column containing quantitative values to be plotted as y-values
    hue_col : str
        column containing categories defining the per-x category split
    hue_order : list[str]
        order in which hue_col categories should be plotted
    show_bars : bool, optional
        If true, shows barplot drawn to mean of each category, by default False
    ax : bool, optional
        If ax object is passed, plots are added to existing axes object. Otherwise, by default False which generates new ax object
    alpha : float, optional
        If both alpha and show_bars, controls the transparency of the barplot patches, by default 0.3
    err_type : str, optional
        What type of error bar to plot - options are standard deviation (std) and standard error of mean (sem), by default 'std'
    colors : bool, optional
        To use custom colours, pass a dictionary mapping each hue category to colour name, by default False and random colours are assigned to each hue category.

    Returns
    -------
    matplotlib.axes
        axes object containing plotted elements. From here, additional axes parameters can readily be adjusted, including labels and title etc.
    """    
    if not show_bars:
        alpha = 0
    if not ax:
        fig, ax = plt.subplots()
    if not colors:
        colors = dict(zip(hue_order, get_cmap(len(hue_order), name='hsv')))

    # Generate figures
    scat = sns.stripplot(data=df, x=x_col, y=y_col, hue=hue_col,
                         hue_order=hue_order, palette=colors, order=x_order, dodge=True, ax=ax, alpha = scat_alpha)

    box = sns.boxplot(data=df.groupby([x_col, hue_col]).mean().reset_index(
    ), x=x_col, y=y_col, hue=hue_col, hue_order=hue_order, palette=colors, order=x_order, ax=ax)

    bar = sns.barplot(data=df, x=x_col, y=y_col, hue=hue_col, hue_order=hue_order,
                      palette=colors, order=x_order, alpha=alpha, errwidth=0, ax=ax,)

    # To generate custom error bars
    bars = bar.patches  # has format [legend*len(hue), hue*xcol]
    num_legend = len(hue_order)

    widths = [bar.get_width() for bar in bars][num_legend:]
    xpos = [bar.get_x() for bar in bars][num_legend:]
    locations = [x + width/2 for x, width in zip(xpos, widths)]

    sample_list = [(hue, x) for hue in hue_order for x in x_order]
    number_groups = len(hue_order)

    # collect mean, sd for each bar
    errvals = dict(df.groupby([hue_col, x_col]).aggregate(
        ['mean', 'std', 'sem'])[y_col].T)
    for sample in sample_list:
        if sample in errvals:
            errvals[sample] = errvals[sample].tolist()
        else:
            errvals[sample] = [np.nan, np.nan, np.nan]
    errvals = pd.DataFrame(errvals, index=['mean', 'std', 'sem'])
    errvals = errvals[sample_list].T

    # add location info
    errvals['xpos'] = locations

    (_, caps, _) = ax.errorbar(x=errvals['xpos'], y=errvals['mean'],
                               yerr=errvals[err_type], capsize=2, elinewidth=1.25, ecolor="black", linewidth=0)
    for cap in caps:
        cap.set_markeredgewidth(2)

    # To only label once in legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[:number_groups],
        labels[:number_groups],
        bbox_to_anchor=(1.0, 1.00),
        title=hue_col,
    )


    return ax



### to plot with SEM
x_order = final['transition_type'].unique().tolist()
sample_sum = len(order)
palette = dict(zip(order, sns.color_palette('mako', sample_sum))) ### adjust number to be the same as the number of treatments you are plotting
font = {'weight' : 'normal','size'   : 12 }
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

fig, ax = plt.subplots()
sns.set_style('whitegrid', {'font_scale':1.5, 'grid.linestyle':'--'})
plot_parameters_scattbar(final, 'transition_type', x_order, 'y_axis', 'treatment', order, show_bars = True, err_type = 'sem', colors = palette, ax = ax)
plt.ylim(0, 50)
plt.ylabel('Transition class')
plt.xticks(rotation = 45)
[x.set_linewidth(2) for x in ax.spines.values()]
[x.set_color('black') for x in ax.spines.values()]
fig.savefig(f'{output_folder}/mean_residence_withSEM.svg', dpi = 600)
plt.show()
