from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob
import os

output_folder = 'Experiment_1-description/python_results'
data_paths = {
    "key":"TDPdata",
    'key':'TDPdata',
    "key":"TDPdata",
    'key':'TDPdata',
}



################
################ Remove outliers from TDP data, calculate the dwell/residence times of each FRET state from HMM analysis and identify transitions. Concatenatea all datasets into a single dataframe.
################

def remove_outliers(compiled_TDP):
    outliers = compiled_TDP[(compiled_TDP["FRET_before"] < -0.5)|(compiled_TDP["FRET_before"] > 1.5)|(compiled_TDP["FRET_after"] < -0.5) | (compiled_TDP["FRET_after"] > 1.5)].index
    compiled_TDP.drop(outliers, inplace = True)
    return compiled_TDP

def calculate_dwell_time(df):
    """Function to convert raw idealized data to a form in which the duration of each idealized state is calculated

    Args:
        df (dataframe): dataframe containing each molecule and the idealized fret 

    Returns:
        dataframe: returns dataframe containing the duration of each FRET state from each molecule
    """
    df_test2 = []
    for Molecule, dfs in df.groupby('Molecule'):
        frame_length = len(dfs)
        test = dfs.groupby([dfs['Idealized_FRET'].ne(dfs['Idealized_FRET'].shift()).cumsum(), 'Idealized_FRET']).size()
        test = test.reset_index(level = 1, drop = False)
        test['Molecule'] = Molecule
        test['number_of_frames'] = frame_length
        df_test2.append(test)
    df_test3 = pd.concat(df_test2)
    df_test3.columns = ['FRET_state', 'Time', 'Molecule', 'number_of_frames']
    df_test3 = df_test3.reset_index().drop('Idealized_FRET', axis = 1)
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

from Utilities.Data_analysis import filter_TDP, file_reader, count_filt_mol


compiled_TDP = []
for data_name, data_path in data_paths.items():
    imported_TDPdata = file_reader(data_path, 'TDP')
    df1 = calculate_dwell_time(imported_TDPdata)
    df2 = generate_transitions(df1)
    cleaned_TDPdata = remove_outliers(df2)
    cleaned_TDPdata["treatment_name"] = data_name
    compiled_TDP.append(cleaned_TDPdata)
compiled_TDP = pd.concat(compiled_TDP)   
compiled_TDP.to_csv(f'{output_folder}/TDP_cleaned.csv', index = False)


#############
############# Code to generate a dataset whereby only molecules that go below a defined threshold are included
#############

filtered_data = filter_TDP(compiled_TDP, 0.3)  ##### number refers to the FRET threshold to filter data. Will only include molecules that go below this set threshold. Will default to 0.3 
filtered_data.to_csv(f'{output_folder}/TDP_cleaned_filt.csv', index = False)

############### 
############### calculate the proportion of molecules that go below a defined threshold
###############

FRET_value = 0.6
molcount = count_filt_mol(compiled_TDP, FRET_value, data_paths, 1)  ## last number refers to the index of the treatment you want to subtract from others
molcount.to_csv(f'{output_folder}/mol_below_{FRET_value}.csv', index = None)

