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
import os as os
import math 

from Utilities.Data_analysis import file_reader

input_folder = 'Experiment_1-description/python_results'
output_folder = f'{input_folder}/Heatmaps-and-first-transitions'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
exposure = 0.200  # in seconds
frame_rate = 1/exposure
FRET_thresh = 0.5 #### Used to filt the threshold for transitions
transition_type = 'high_to_low'   ##### select either 'low_to_high' or 'high_to_low' for plotting first specified transition
time_thresh  = 15 #### in seconds

data_paths = {
    'treatment':('treatment_name', 'Experiment_2-description/raw_data/treatment_to_plot'),
}

################
################ Plot heatmaps
################

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

def plot_heatmap(df, gridsize, bins_hex, plot_type = 'hex'):
    for treatment, dfs in df.groupby('treatment_name'):
        if plot_type == 'hex':
            plt.rcParams['svg.fonttype'] = 'none'
            sns.set(style="whitegrid")
            g = sns.JointGrid(data = dfs, x='time', y='FRET', xlim = (0,100), ylim = (0, 1))
            g.plot_joint(plt.hexbin, gridsize=(gridsize, gridsize), cmap='ocean_r', mincnt=0, bins=bins_hex)
            g.plot_marginals(sns.histplot, kde=True, bins=20)
            plt.savefig(f'{output_folder}/Heatmap_{treatment}_{plot_type}.svg', dpi = 600)
        if plot_type == 'kde':
            plt.rcParams['svg.fonttype'] = 'none'
            sns.set_style("whitegrid", {'grid.linestyle':'--', 'axes.linewidth':20, 'axes.color':'black', 'axes.edgecolor': 'black', 'font.size':10})
            fig = sns.jointplot(data = dfs, x='time', y='FRET', xlim = (0,200), ylim = (0, 1), alpha = 0.05, color = '#2AA6CF', marginal_kws = dict(bins = 20, kde = True))   
            fig.plot_joint(sns.kdeplot, cmap = 'mako') 
            fig.ax_joint.spines['top'].set_visible(True)
            fig.ax_joint.spines['right'].set_visible(True)
            fig.ax_marg_x.spines['left'].set_visible(True)
            fig.ax_marg_x.spines['top'].set_visible(True)
            fig.ax_marg_x.spines['right'].set_visible(True)
            fig.ax_marg_y.spines['top'].set_visible(True)
            fig.ax_marg_y.spines['right'].set_visible(True)
            fig.ax_marg_y.spines['bottom'].set_visible(True)
            plt.savefig(f'{output_folder}/Heatmap_{treatment}_{plot_type}.svg', dpi = 600)
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

plot_heatmap(compiled_df, 100, 80, 'kde')


no_filt = compiled_df['treatment_name'].unique().tolist()

def plot_average_FRET_over_time(df, filt, ci = 'sd', x_axis = 'time', subplot = False):
    sns.set(style = "ticks", font_scale = 1.5)
    filt_data = df[df['treatment_name'].isin(filt)]
    if x_axis == 'time':
        fig, ax = plt.subplots()
        sns.lineplot(data=filt_data, x = x_axis, y = 'FRET', ci = ci, hue = 'treatment_name')
        plt.xlim(0, 200, 10)
        ax.axvline(x = 10, linestyle = '-', color = 'grey')
        ax.axvspan(0, 10, facecolor='grey', alpha = .2)   
        plt.ylim(0, 1, 10)
        plt.xlabel('Time (s)')
        plt.legend(fontsize = 'small')
        plt.savefig(f'{output_folder}/Average_Heatmap_{filt}.svg', dpi = 600)
    if (x_axis == 'normalised_to_event' and subplot == True):
        nrow = math.ceil(len(filt)/2)
        fig, axes = plt.subplots(nrow, 2, figsize = (8, 3*nrow), sharex = True, sharey = True)
        axes = axes.flatten()
        for i, treatment in enumerate(filt):
            sns.lineplot(data=df[df['treatment_name'] == treatment], x = x_axis, y = 'FRET', ci = ci, hue = 'treatment_name', ax = axes[i])
            axes[i].axvline(x = 0, linestyle = '--', color = 'grey')
            plt.xlim(-30, 30, 10)
            plt.ylim(0, 1, 10)
            axes[i].set_xlabel('Time (s)')
            axes[i].set_title(f"{treatment} (n = {str(normalised_data[normalised_data['treatment_name'] == treatment]['molecule_number'].nunique())})")
            axes[i].get_legend().remove()
        if len(filt) < len(axes):
            axes[-1].remove()
            axes[3].set_xlabel('Time (s)')
        plt.savefig(f'{output_folder}/Traces_normalised_to_first_trans.svg', dpi = 600)
        plt.show()
    elif x_axis == 'normalised_to_event':
        for treatment, df in filt_data.groupby('treatment_name'):
            sns.lineplot(data=df[df['treatment_name'] == treatment], x = x_axis, y = 'FRET', ci = ci, hue = 'treatment_name')
            plt.xlim(-30, 30, 10)
            plt.axvline(x = 0, linestyle = '--', color = 'grey')
            plt.ylim(0, 1, 10)
            plt.xlabel('Time (s)')
            plt.legend(fontsize = 'small', loc = 'upper left')
            plt.savefig(f'{output_folder}/Traces_normalised_to_first_trans_{treatment}.svg', dpi = 600)
            plt.show()
    plt.show()
    return
###############
############### Generate transitions, filter for the first transition that meets criteria and plot
###############

def calculate_dwell_time(df):
    """Function to convert raw idealized data to a form in which the duration of each idealized state is calculated

    Args:
        df (dataframe): dataframe containing each molecule and the idealized fret 

    Returns:
        dataframe: returns dataframe containing the duration of each FRET state from each molecule
    """
    df_test2 = []
    for Molecule, dfs in df.groupby('molecule_number'):
        frame_length = len(dfs)
        test = dfs.groupby([dfs['idealized FRET'].ne(dfs['idealized FRET'].shift()).cumsum(), 'idealized FRET']).size()
        test = test.reset_index(level = 1, drop = False)
        test['Molecule'] = Molecule
        test['number_of_frames'] = frame_length
        df_test2.append(test)
    df_test3 = pd.concat(df_test2)
    df_test3.columns = ['FRET_state', 'Time', 'Molecule', 'number_of_frames']
    df_test3 = df_test3.reset_index().drop('idealized FRET', axis = 1)
    return df_test3[df_test3.groupby('Molecule').Molecule.transform('count') > 1]

def generate_transitions(df):
    """Converts the duration of each FRET state into a transition, whereby the FRET state before, the FRET state after
    and the duration of the FRET state before a transition is given in a single line. Each line represents a single transition.

    Args:
        df (dataframe): dataframe generated following 'calculate_dwell_time' function in which the duration of a certain
        FRET state is given for each FRET state for all molecules

    Returns:
        dataframe: returns a dataframe in which each row represents a transition, with FRET before transition, FRET after transition
        and duration of FRET state before transition (given in number of frames in column Time) provided
    """
    df_toconcat = []
    for molecule, dfs in df.groupby('Molecule'):
        thing1 = dfs.assign(FRET_after = dfs.FRET_state.shift(-1)).dropna()
        df_toconcat.append(thing1)
    compiled_df = pd.concat(df_toconcat).reset_index(drop = True)
    compiled_final = compiled_df[['Molecule', 'FRET_state', 'FRET_after', 'Time', 'number_of_frames']]
    compiled_final.columns = ['Molecule', 'FRET_before', 'FRET_after', 'Time', 'number_of_frames']
    return compiled_final

def remove_outliers(compiled_TDP):
    outliers = compiled_TDP[(compiled_TDP["FRET_before"] < -0.5)|(compiled_TDP["FRET_before"] > 1.5)|(compiled_TDP["FRET_after"] < -0.5) | (compiled_TDP["FRET_after"] > 1.5)].index
    compiled_TDP.drop(outliers, inplace = True)
    return compiled_TDP

def filter_FRET_trans_if(dfs, thresh, trans_type = 'low_to_high'):
    """This function has several roles. Firstly, for each molecule in the dataframe it will add a column with the cumulative
    sum of all residence times that will be used in later functions. Secondly, depending on what kind of transitions
    you are interested in, it will filter the dataset to include only those transitions (e.g., low to high)

    Args:
        dfs (dataframe): Dataframe containing the molecules, fret before, fret after, idealized fret, fret and time
        thresh (value): The FRET value at which to set the transition threshold. Will only find those that are
        lower than the thresh going to above the thresh (or vice versa)
        trans_type (str, optional): _description_. Defaults to 'low_to_high'. Dictates if youu want to look at high-to-low
        or low-to-high transitions. This is set as variable at the top of the script.

    Returns:
        dataframe: Will return a filtered dataframe with the transitions of interest as well as the cumulative sum of time
        at which each transition occurs (essentially how long into imaging does the transition appear)
    """
    comb = []
    for Molecule, df in dfs.groupby('Molecule'):
        df['cum_sum'] = df['time (s)'].cumsum()
        comb.append(df)
    combined = pd.concat(comb)
    if trans_type == 'high_to_low':
        filt_data = combined[(combined['FRET_before'] > thresh) & (combined['FRET_after'] < thresh)]
    elif trans_type == 'low_to_high':
        filt_data = combined[(combined['FRET_before'] < thresh) & (combined['FRET_after'] > thresh)]
    return filt_data

def select_first_transition(dfs, time_thresh):
    """Will find the first transition for each molecule. Important to note that this function should be run after the 
    'filter_FRET_trans_if' function, which filters for only those transitions that meet a criteria. This function 
    will essentially then find the first transition for a molecule that meets a defined criteria (e.g., first low-to-high 
    transition)

    Args:
        dfs (dataframw): Filtered dataframe containing only transitions of interest. Same as that returned after 
        executing 'filter_FRET_trans_if'

    Returns:
        dataframe: Dataframe containing the first transition of each molecule and the cumulative time that this occurs.
    """
    first_trans = []
    for molecule, df in dfs.groupby('Molecule'):
        if df['cum_sum'].min() < time_thresh:
            continue
        first_trans_above_timethresh = df[df['cum_sum'] == df['cum_sum'].min()]
        first_trans.append(first_trans_above_timethresh) 
    return pd.concat(first_trans)



def plot_first_specified_transition(df, trans_type):
    plot1 = plt.figure()
    sns.set_style('ticks')
    sns.violinplot(data = df, y = 'cum_sum', x = 'treatment', scale = 'width')
    sns.stripplot(data = df, y = 'cum_sum', x = 'treatment', color ='black', alpha = 0.5)
    plt.ylabel('Time for RNA binding (s)')
    plot1.savefig(f'{output_folder}/time_until_first_{trans_type}_transition.svg', dpi = 600)


def normalise_to_event(df1, df2):
    """ This function uses two dataframes to normalise the x-axis for each molecule so that the first transition
    that meets a criteria (filtered for and identified using the 'filter_FRET_trans_if' and 'select_first_transition'
    functions) is set to 0. This should allow the first transition between molecules to be synchronised to potentially
    observe changes in FRET that occur immediately prior to or after the transition that is normally hidden by
    the asynchronous timing of transitions between molecules.

    Args:
        df1 (dataframe): Dataframe that has not been filtered. Contains all transitions for all molecules and treatments
        df2 (dataframe): Filtered dataframe. Contains only the first transition for each molecule that meets the criteria.
        Also contains the cumulative sum of time at which that transition occurs, which is then subtracted from all
        timepoints for the corresponding molecule in df1 

    Returns:
        dataframe: returns a dataframe that contains all molecules that contain the transition of interest with an extra 
        column containing the normalised time (time for molecule minus time of first transition)
    """
    collated = []
    for (treatment, mol), df in df2.groupby(['treatment', 'Molecule']):
        testies = df1[(df1['treatment_name'] == treatment) & (df1['molecule_number'] == mol)]
        testies['normalised_to_event'] = testies['time']-(float(df[(df['Molecule'] == mol) & (df['treatment'] == treatment)]['cum_sum']))
        if transition_type == 'low_to_high' and testies['idealized FRET'].iloc[0] <= FRET_thresh:
            collated.append(testies)
        elif transition_type == 'high_to_low' and testies['idealized FRET'].iloc[0] >= FRET_thresh:
            collated.append(testies)
    return pd.concat(collated)


compiled_filt = []
for treatment, df in compiled_df.groupby('treatment_name'):
    treatment_df = compiled_df[compiled_df['treatment_name'] == treatment]
    treatment_df2 = treatment_df.filter(items = ['idealized FRET','molecule_number'])
    treatment_df3 = calculate_dwell_time(treatment_df2)
    treatment_transitions = generate_transitions(treatment_df3)
    treatment_cleaned_transitions = remove_outliers(treatment_transitions)
    treatment_cleaned_transitions['time (s)'] = treatment_cleaned_transitions['Time'] * exposure
    treatment_cumsum = filter_FRET_trans_if(treatment_cleaned_transitions, FRET_thresh, transition_type) ##### add 'high_to_low' to look at how long it takes for high-low transitions occur
    treatment_first_transition = select_first_transition(treatment_cumsum, time_thresh)
    treatment_first_transition['treatment'] = treatment
    compiled_filt.append(treatment_first_transition)
col = pd.concat(compiled_filt)
normalised_data = normalise_to_event(compiled_df, col)

plot_first_specified_transition(col, transition_type)

print(col.groupby('treatment')['cum_sum'].mean())
print(col.groupby('treatment')['cum_sum'].sem())

##########
########## Plot datasets
##########

to_filt = ['TREATMENT'] ######### select if you want what datasets to include in average heatmap


plot_first_specified_transition(col, transition_type)
plot_average_FRET_over_time(normalised_data, to_filt, 'sd', 'time') ##### change to 'to_filt' to include only datasets that were mentioned above
plot_average_FRET_over_time(normalised_data, to_filt, 'sd', 'normalised_to_event', False)  #### add True to plot subplots
