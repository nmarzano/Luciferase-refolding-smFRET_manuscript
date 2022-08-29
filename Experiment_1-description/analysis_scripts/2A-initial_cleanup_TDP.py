from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob
import os
from Utilities.Data_analysis import filter_TDP, file_reader, count_filt_mol

output_folder = 'Experiment_1-description/python_results'
compiled_data = pd.read_csv(f'{output_folder}/Cleaned_FRET_histogram_data.csv')
FRET_value = 0.5 #### the proportion of molecules that travel below this threshold will be counted

################
################ Remove outliers from TDP data, calculate the dwell/residence times of each FRET state from HMM analysis and identify transitions. Concatenatea all datasets into a single dataframe.
################


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

compiled_filt = []
for treatment, df in compiled_data.groupby('treatment_name'):
    treatment_df = compiled_data[compiled_data['treatment_name'] == treatment]
    treatment_df2 = treatment_df.filter(items = ['idealized FRET','molecule_number'])
    treatment_dwell = calculate_dwell_time(treatment_df2)
    treatment_trans = generate_transitions(treatment_dwell)
    treatment_cleaned = remove_outliers(treatment_trans)
    treatment_cleaned["treatment_name"] = treatment
    compiled_filt.append(treatment_cleaned)
compiled_TDP = pd.concat(compiled_filt)   
compiled_TDP.to_csv(f'{output_folder}/TDP_cleaned.csv', index = False)


#############
############# Code to generate a dataset whereby only molecules that go below a defined threshold are included
#############

filtered_data = filter_TDP(compiled_TDP, 0.3)  ##### number refers to the FRET threshold to filter data. Will only include molecules that go below this set threshold. Will default to 0.3 
filtered_data.to_csv(f'{output_folder}/TDP_cleaned_filt.csv', index = False)

############### 
############### calculate the proportion of molecules that go below a defined threshold
###############

def filter_tdp(dfs, FRET_value):
    filt_tdp = []
    for treatment, df in dfs.groupby('treatment_name'):
        mol_list = df[(df["FRET_before"] <= FRET_value)|(df["FRET_after"] <= FRET_value)].Molecule.unique().tolist()
        filtered = df[df["Molecule"].isin(mol_list)]
        filt_tdp.append(filtered)
    return pd.concat(filt_tdp)


def count_filt_mol(dfs, thresh):
    """Will count the number of molecules in which the idealized FRET will go below a defined threshold at some point before photobleaching

    Args:
        df (dataframe): Contains all the data required to plot TDP and identify transitions below threshold (i.e., FRET, idealized FRET, molecule)
        thresh (float): Threshold to set. If set to 0.5, function will count the number of molecules that go below 0.5 at some point
        dataname (dict): Dictionary containing keys for each treatment - used to find mol count for each treatment

    Returns:
        dataframe: Will return dataframe with raw mol count and also corrected mol count. Corrected mol count is calculated as the Raw mol count subtracted
        by the molcount of another treatment. The treatment to subtract is defined by 'order', which is the index of the treatment you want to subtract
    """
    filtered_data = filter_tdp(dfs, thresh)
    percent_mol_concat = []
    for treatment, df in dfs.groupby('treatment_name'):
        total_mol = len(df.Molecule.unique())
        filt_mol = len(filtered_data[filtered_data['treatment_name']==treatment].Molecule.unique())
        df['proportion_mol_below_thresh'] = (filt_mol/total_mol)*100
        percent_mol_concat.append(df)
    return pd.concat(percent_mol_concat)

proportion_mol_below_thresh = pd.DataFrame(count_filt_mol(compiled_TDP, FRET_value).groupby('treatment_name')['proportion_mol_below_thresh'].mean())
proportion_mol_below_thresh['norm_percent_mol'] = proportion_mol_below_thresh['proportion_mol_below_thresh'] - proportion_mol_below_thresh['proportion_mol_below_thresh'].iloc[0]
proportion_mol_below_thresh.reset_index()
proportion_mol_below_thresh.to_csv(f'{output_folder}/mol_below_{FRET_value}.csv', index = None)

