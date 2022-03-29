from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob
import Utilities.Data_analysis as uda

output_folder = "Experiment_X-description/python_results" #### NEED TO CHANGE EACH TIME TO THE PYTHON RESULTS LEVEL

FRET_thresh = 0.5 #### FRET value at which to filter data above or below. IF CHANGED, WILL NEED TO CHANGE ALL 0.5 VALUES (E.G. BELOW IN HEADERS) TO THE NEW VALUE
headers_withsum =  [f"< {FRET_thresh} to < {FRET_thresh}", f"< {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to < {FRET_thresh}", "sum", "sample"]
order = [4, 1, 2, 3, 0]  # organise data using its index 

imported_df = uda.file_reader(f'{output_folder}/Dwell_frequency', 'transition_frequency', False, headers_withsum) ##### NEED TO ALTER EACH TIME
imported_df.drop(columns = ["sum"], inplace = True)
imported_df = imported_df.reindex(order) ### need to do this to plot in order that you would like
imported_df.to_csv(f'{output_folder}/Transition_frequency.csv')


######################## code to prepare for plotting
total = [i+j+k+m for i,j,k, m in zip(imported_df[f'< {FRET_thresh} to < {FRET_thresh}'], imported_df[f'< {FRET_thresh} to > {FRET_thresh}'], imported_df[f'> {FRET_thresh} to > {FRET_thresh}'], imported_df[f"> {FRET_thresh} to < {FRET_thresh}"])]
lowtolow = [i / j * 100 for i,j in zip(imported_df[f'< {FRET_thresh} to < {FRET_thresh}'], total)]
lowtohigh = [i / j * 100 for i,j in zip(imported_df[f'< {FRET_thresh} to > {FRET_thresh}'], total)]
hightohigh = [i / j * 100 for i,j in zip(imported_df[f'> {FRET_thresh} to > {FRET_thresh}'], total)]
hightolow = [i / j * 100 for i,j in zip(imported_df[f'> {FRET_thresh} to < {FRET_thresh}'], total)]

################### plotting code
labels = imported_df["sample"].to_list()
barWidth = 0.85
plt.rcParams['svg.fonttype'] = 'none'
sns.set(style = "darkgrid", color_codes = 'pastel')
# Create green Bars
plot1 = plt.figure()
plt.bar(labels, lowtolow, color='skyblue', edgecolor='white', width=barWidth, label = f"< {FRET_thresh} to < {FRET_thresh}" )
# Create orange Bars
plt.bar(labels, lowtohigh, bottom=lowtolow, color='royalblue', edgecolor='white', width=barWidth, label = f"< {FRET_thresh} to > {FRET_thresh}")
# Create blue Bars
plt.bar(labels, hightohigh, bottom=[i+j for i,j in zip(lowtolow, lowtohigh)], color='mediumpurple', edgecolor='white', width=barWidth, label = f"> {FRET_thresh} to > {FRET_thresh}")
# Create blue Bars
plt.bar(labels, hightolow, bottom=[i+j+k for i,j,k in zip(lowtolow, lowtohigh, hightohigh)], color='indigo', edgecolor='white', width=barWidth,label = f"> {FRET_thresh} to < {FRET_thresh}")
plt.legend(loc = "upper left", bbox_to_anchor = (1,1), ncol =1)
plt.ylabel("Transition probability (%)")
plt.xticks(rotation=45)
plot1.savefig(f'{output_folder}/Transition_frequency_plot.svg', dpi = 600)



