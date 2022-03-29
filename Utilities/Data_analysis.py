import pandas as pd
import numpy as np
import glob as glob
import seaborn as sns
import matplotlib.pyplot as plt

def filter_TDP(data_frame, thresh = 0.3):
    """[optional function that removes molecules that do not transition below threshold at some point]

    Args:
        data_frame ([dataframe]): [dataframe that has been cleaned to remove outliers by 'remove_outliers' function]
        thresh (float, optional): [FRET value to threshold - will filter molecules and keep those that go below the thresh]. Defaults to 0.3.

    Returns:
        [dataframe]: [contains only molecules of interest]
    """
    filtered_mol = []
    for treatment, df in data_frame.groupby("treatment_name"):
        mol_list = df[(df["FRET_before"] <= thresh)|(df["FRET_after"] <= thresh)].Molecule.unique().tolist()
        filtered = df[df["Molecule"].isin(mol_list)]
        filtered_mol.append(filtered)
    filtered_mol = pd.concat(filtered_mol)
    return filtered_mol

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
            rawFRET = compiled[(compiled[3] > -0.5) & (compiled[3] < 1.5)].copy()
            return rawFRET
        if data_type == "idealized":
            idealizedFRET = compiled[(compiled[4] > -0.5) & (compiled[4] < 1.5)].copy()
            return idealizedFRET
    elif plot_type == 'TDP':
        outliers = compiled[(compiled["FRET before transition"] < -0.5)|(compiled["FRET before transition"] > 1.5)|(compiled["FRET after transition"] < -0.5) | (compiled["FRET after transition"] > 1.5)].index
        compiled.drop(outliers, inplace = True)
        return compiled
    else:
        print('invalid plot type, please set plot_type as "hist" or "TDP" - you idiot')

def cleanup_dwell(data, fps, thresh, first_dwell = "delete"):
    """[Will convert the data frome frame number to unit of time (seconds) and then delete all dwell times
    that are smaller than the set threshold (defined previously) in seconds. Will also delete the first dwell 
    state from each molecule]

    Args:
        data ([dataframe]): [raw data]
        first_dwell (str, optional): [Set to 'keep' to keep the first dwell state from each molecule otherwise 
        Will delete the first dwell state from each molecule by default]. Defaults to "delete".

    Returns:
        [dataframe]: [Data is now cleaned and ready for subsequent processing]
    """
    if first_dwell == "delete":
        filtered = []
        for molecule, df in data.groupby("Molecule"):
            filtered.append(df.iloc[1:])
        filtered = pd.concat(filtered)  #####filtered = pd.concat([df.iloc[1:] for molecule, df in A.groupby("Molecule")]) ##code here is the same as the for loop but in a list comprehension format
        filtered["Time (s)"] = filtered["Time"]/fps
        filtered = filtered[filtered["Time (s)"] >= thresh]
        return filtered
    if first_dwell == "keep":
        data["Time (s)"] = data["Time"]/fps
        data = data[data["Time (s)"] >= thresh]
        return data

def filter_dwell(df, FRET_thresh, headers):
    """[Will take the cleaned TDP data and will filter it using a threshold (defined by FRET_thresh)
    into seperate types of transitions (e.g., < 0.5 to > 0.5 FRET if FRET_thresh is = 0.5 is one example
    of a type of transition).

    Args:
        df ([dataframe]): [contains cleaned data that has been processed using the 'cleanup_dwell' function]

    Returns:
        [dataframe]: [contains dwell time that has been categorized into each transition class]
    """
    filtered_lowtohigh = df[(df["FRET_before"] < FRET_thresh) & (df["FRET_after"] > FRET_thresh)].copy()
    filtered_lowtolow = df[(df["FRET_before"] < FRET_thresh) & (df["FRET_after"] < FRET_thresh)].copy()
    filtered_hightolow = df[(df["FRET_before"] > FRET_thresh) & (df["FRET_after"] < FRET_thresh)].copy()
    filtered_hightohigh = df[(df["FRET_before"] > FRET_thresh) & (df["FRET_after"] > FRET_thresh)].copy()
    dataf = [filtered_lowtolow["Time (s)"], filtered_lowtohigh["Time (s)"], filtered_hightohigh["Time (s)"], filtered_hightolow["Time (s)"]]
    df_col = pd.concat(dataf, axis = 1, keys = headers)
    df_col = df_col.apply(lambda x:pd.Series(x.dropna().values))  ## removes NaN values from each column in df_col
    return df_col

def transition_frequency(filt):
    """calculates the transition frequency (i.e., the number of transitions per transition class divided
    by the total number of transitions observed). For example if there are 40 transitions total, and a 
    < 0.5 to > 0.5 transition occurs 10 times, then the transition probability for that transition type is 
    0.25 or 25%.

    Args:
        filt (dataframe): contains the dataframe with filtered data (cleaned data has been filtered by
        'filter_dwell' function)

    Returns:
        dataframe: returns a dataframe containing the percentage for each transition type
    """
    count_df_col = pd.DataFrame(filt.count(axis = 0)).transpose()
    count_df_col["sum"] = count_df_col.sum(axis = 1)
    dwell_frequency = pd.DataFrame([(count_df_col[column]/count_df_col["sum"])*100 for column in count_df_col]).transpose()
    print(dwell_frequency)
    return dwell_frequency

def calculate_mean(filtered_data, treatment_name):
    """calculates the mean dwell time of each type of transition class

    Args:
        filtered_data (dataframe): dataframe generated after the 'cleanup_dwell' and 'filter_dwell' functions 
        have been run
        treatment_name (str): not required, only present to receive input from for loop. set to treatment_name
    Returns:
        [dataframe]: returns dataframe containing the mean of each transition class
    """
    mean_dwell = pd.DataFrame([filtered_data.iloc[0:].mean()])
    mean_dwell["sample"] = treatment_name
    return mean_dwell

def float_generator(data_frame, treatment, FRET_thresh):
    """Will generate the float values used to scale the size of the arrows when plotting the summary heatmap    

    Args:
        data_frame (dataframe): Takes dataframe containing the transition frequencies generated using the 
    'transition_frequency' function
        treatment (str): will take 'treatment' from for loop. Leave as treatment.

    Returns:
        [dataframe]: returns values used to scale arrow - 
    """
    transition_frequency_arrow = data_frame[data_frame["sample"]== treatment]
    normalised_number_lowtohigh = float(np.array(transition_frequency_arrow[f"< {FRET_thresh} to > {FRET_thresh}"])/1000)
    normalised_number_hightolow = float(np.array(transition_frequency_arrow[f"> {FRET_thresh} to < {FRET_thresh}"])/1000)
    normalised_number_hightohigh = float(np.array(transition_frequency_arrow[f"> {FRET_thresh} to > {FRET_thresh}"])/1000)
    normalised_number_lowtolow = float(np.array(transition_frequency_arrow[f"< {FRET_thresh} to < {FRET_thresh}"])/1000)
    arrow_list = [normalised_number_lowtohigh,normalised_number_hightolow,normalised_number_hightohigh,normalised_number_lowtolow]
    return arrow_list

def heatmap_prep(histogram_data, treatment, FRET_thresh):
    """Takes the data and calculates the number of data points total (total), how much data points are below
    a threshold (time below) and above a threshold (time above) and will use these values to calculate what proportion
    of time the FRET is below or above thresh. Will then feed into the heatmap when plotting

    Args:
        histogram_data (dataframe): data used to plot the histogram - will include all the cleaned FRET and
    idealized FRET values (time not a factor here, just number of frames)
        treatment (str): only used to be fed from for loop. do not change

    Returns:
        df: contains the data required to plot the heatmap. will be as a proportion of time spent below or above
        threshold
    """
    subset_data = histogram_data[histogram_data["treatment_name"]==treatment]
    total = len(subset_data[(subset_data["FRET"] < FRET_thresh) | (subset_data["FRET"] > FRET_thresh)])
    subset_data_largerthanthresh = len(subset_data[subset_data["FRET"] > FRET_thresh])
    subset_data_lessthanthresh = len(subset_data[subset_data["FRET"] < FRET_thresh])
    time_below = (subset_data_lessthanthresh/total)
    time_above = (subset_data_largerthanthresh/total)
    thresh_dicts = {f"< {FRET_thresh}":[time_below], f"> {FRET_thresh}":[time_above] }
    thresh_dfs = pd.DataFrame(thresh_dicts)
    #thresh_dfs["treatment"] = treatment
    return thresh_dfs

def mean_dwell_prep(mean_dwell_data, treatment, FRET_thresh):
    """calculates the mean values for each transition class and converts to float values for plotting in
    summary heatmap 

    Args:
        mean_dwell_data (dataframe): dataframe containing mean values
        treatment (str): used for a for loop. do not change.

    Returns:
        [list]: contains list of float values
    """
    subset_mean_dwell = mean_dwell_data[mean_dwell_data["sample"]== treatment]
    meandwell_lowtohigh = float(np.array(subset_mean_dwell[f"< {FRET_thresh} to > {FRET_thresh}"]))
    meandwell_hightolow = float(np.array(subset_mean_dwell[f"> {FRET_thresh} to < {FRET_thresh}"]))
    meandwell_hightohigh = float(np.array(subset_mean_dwell[f"> {FRET_thresh} to > {FRET_thresh}"]))
    meandwell_lowtolow = float(np.array(subset_mean_dwell[f"< {FRET_thresh} to < {FRET_thresh}"]))
    mean_list = [meandwell_lowtohigh,meandwell_hightolow,meandwell_hightohigh,meandwell_lowtolow]
    return mean_list

def file_reader(input_folder, data_type, frame_rate = False, column_names = False): 
    """will import data

    Args:
        input_folder (directory): where data is stored
        data_type (str): what data will be used to plot, needs to be either 'hist', 'TDP', 'transition_frequency'
        or 'other'. 

    Returns:
        dataframe: dataframe with data to be used in subseqeunt codes
    """
    if data_type == 'hist':
        filenames = glob.glob(input_folder + "/*.dat")
        dfs = []
        for filename in filenames:
            molecule_number = filename.split('\\')[1].split('_')[0]
            hist_data = pd.read_table(filename, sep="\s+", header=None)
            hist_data['molecule number'] = molecule_number
            dfs.append(hist_data) ### will error if forward slash (e.g. "/s+")
        test = pd.concat(dfs)
        test_dfs = pd.DataFrame(test)
        return test_dfs
    elif data_type == 'TDP':
        filename = input_folder
        A = pd.read_table(filename, header = None, sep="\s+")
        A.columns = ['Molecule', 'Idealized_FRET']
        return A
    elif data_type == 'transition_frequency':
        if not column_names:
            print('no column_names found, specify list to use')
            return
        filenames = glob.glob(input_folder + "/*.csv")
        dfs = []
        for filename in filenames:
            dfs.append(pd.read_csv(filename, header=None))
        test = pd.concat(dfs, ignore_index=True)
        test_dfs = pd.DataFrame(test)
        test_dfs.columns = column_names
        return test_dfs
    elif data_type == 'heatmap':
        filenames = glob.glob(input_folder + "/*.dat")
        dfs = []
        for filename in filenames:
            name = filename.split('\\')[1].split('_')[0]
            mol_data = pd.read_table(filename, sep="\s+", header=None)
            mol_data.columns = ['frame', 'donor', 'acceptor', 'FRET', 'idealized_FRET']
            mol_data['molecule'] = name
            mol_data['time'] = mol_data['frame']/frame_rate
            dfs.append(mol_data) ### will error if forward slash (e.g. "/s+")
        test = pd.concat(dfs)
        test_dfs = pd.DataFrame(test)
        return test_dfs
    elif data_type == 'other':
        dfs = pd.read_csv(input_folder)
        return dfs
    else:
        print('invalid data_type, please set data_type as "hist", "TDP","transition_frequency" or "other" if using for violin or heatmap plots')

def count_filt_mol(df, thresh, dataname, order):
    """Will count the number of molecules in which the idealized FRET will go below a defined threshold at some point before photobleaching

    Args:
        df (dataframe): Contains all the data required to plot TDP and identify transitions below threshold (i.e., FRET, idealized FRET, molecule)
        thresh (float): Threshold to set. If set to 0.5, function will count the number of molecules that go below 0.5 at some point
        dataname (dict): Dictionary containing keys for each treatment - used to find mol count for each treatment

    Returns:
        dataframe: Will return dataframe with raw mol count and also corrected mol count. Corrected mol count is calculated as the Raw mol count subtracted
        by the molcount of another treatment. The treatment to subtract is defined by 'order', which is the index of the treatment you want to subtract
    """
    filtered_data = filter_TDP(df, thresh)
    data_paths = dataname
    percent_mol_concat = {}
    for data_name, data_path in data_paths.items():
        total_mol = len(df[df['treatment_name']==data_name].Molecule.unique())
        filt_mol = len(filtered_data[filtered_data['treatment_name']==data_name].Molecule.unique())
        percent_mol = (filt_mol/total_mol)*100
        percent_mol_concat[data_name] = percent_mol
    percent_mol_concat = pd.DataFrame(percent_mol_concat, index = ['percent_mol']).T.reset_index().rename(columns={'index':'treatment'})
    percent_mol_concat['norm_percent_mol'] = percent_mol_concat['percent_mol'] - percent_mol_concat.iloc[order,1]
    return percent_mol_concat

def fret_state_trans(dfs, thresh, fps_clean, thresh_clean, state = 'after'):
    """Prepares a dataset in which 
    Will plot a violin plot of all the FRET states immediately prior to a transition to another FRET state that is below a defined threshold

    Args:
        dfs (dataframe): Contains all the data required to plot TDP and identify transitions below threshold (i.e., FRET, idealized FRET, molecule)
        thresh (float): Threshold that defines the FRET state that you want to look at. For example, if you want to look at the FRET state immediately priort
        to a transition below 0.3 FRET then you will set 'thresh' as 0.3 
        fps_clean ([type]): Required for cleanup_dwell function. Needs this to convert frames to seconds
        thresh_clean ([type]): Required for cleanup_dwell function. Specifies the minimum residence time. All residence times less than thresh_clean will be deleted

    Returns:
        dataframe: Dataframe containing all the transitions in which the 'FRET after transition' is below 'thresh'
    """
    cleaned_df = []
    for treatment_name, df in dfs.groupby("treatment_name"):
        initial_data = df[df["treatment_name"] == treatment_name]    
        cleaned = cleanup_dwell(initial_data, fps_clean, thresh_clean)
        cleaned_df.append(cleaned)
    cleaned_concat = pd.concat(cleaned_df)
    filt = []
    for treatment_name, df in cleaned_concat.groupby("treatment_name"):
        if state == 'after':
            filt.append(df[df['FRET_before'] <= thresh])
        elif state == 'before':
            filt.append(df[df['FRET_after'] <= thresh])
    filtered_f_state = pd.concat(filt)
    return filtered_f_state





