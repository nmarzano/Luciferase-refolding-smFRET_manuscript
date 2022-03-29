# Manuscript v1.01: Code used to analyze data published in 'Real-time single-molecule observation of chaperone-assisted protein folding'

This code was written to analyse single-molecule FRET (smFRET) data generated as part of the manuscript entitled ‘Real-time single-molecule observation of chaperone-assisted protein folding'. Generally, the code imports the raw FRET data, removes outliers and analyses the kinetic details originating from Hidden Markov Model (HMM) fitting of the data using the vbFRET MATLAB program. Furthermore, analysed data can be plotted in a variety of different formats (e.g., FRET histograms, transition density plots [TDPs], FRET-intensity traces, violin plots etc.) for easy data visualization.

Prerequisites: 

Two kinds of data are required for this analysis. 

    (1)	The idealized FRET data for individual molecules following HMM analysis using vbFRET. These are generated following fitting and saving the data as ‘Idealized Traces’ with ‘Individual Text Files’ (in .dat format). Data for individual molecules is exported from vbFRET and stored within a folder.
    (2)	The idealized FRET data stored as a concatenated file. These are generated following fitting and saving the data as ‘Idealized Traces’ with ‘Concatenated Text Files’ (in .dat format).  This dataset is for generating TDP data, which includes the initial FRET state [prior to a transition] and the final FRET state [FRET state after the transition] and the corresponding number of frames that the molecule was in the initial FRET state prior to a transition. 

Workflow: 

To reproduce analyses presented in the manuscript, it is recommended to use the following order: Numbered items in order (0 - 3) and within each numbered analysis scripts to follow the lettered order (A- Z). This will ensure any source data is generated as needed (although in some cases this is not essential).

    0-import_data - imports data into Visual Studio Code [vscode] environment from computer directory). 

    1A-plot-histogram - imports data from vscode environment, removes major outliers and plots histograms/ridgeline plots. 

    1B-plot_traces - used to plot any individual molecule of interest. Just need to use copy and paste in the relative file path of the specific molecule you wish to plot and ensure the correct treatment ‘key’ for your dataset of interest is provided.

    1C-heatmap - used to plot a FRET heatmap that concatenates all molecules from a particular dataset and plots the FRET intensity over time. Can plot multiple datasets simultaneously, just define the ‘key’ and provide the directory to the folder with the individual files for each molecule.

    1D-gaussian-fitting – used to fit histogram data with a Gaussian function. Data will be fit to a 3-gaussian model and will output plots demonstrating the fit and the relative proportion of each Gaussian.

    2A-initial_cleanup_TDP - cleans TDP data and removes outliers, finds the proportion of molecules from each treatment that goes below a defined threshold. 

    2B-plot_TDP - plots TDP plots. 

    2C-dwell_time_analysis - cleans all TDP data, removes outliers and removes the last residence time for each molecule. Can be used to generate violin plots of the FRET state before a transition to below a defined threshold, set a threshold and calculate the transition frequency of specific classes of transitions, calculate the number of times a molecule transitions above and below a defined threshold and the proportion of transitions per treatment that is larger than a defined change in FRET. 

    2D-pot_transition_frequency - plots transition frequency graphs. 

    2E-plot_violin_plots - plots violin plots in a variety of options/formations. Also plots the mean and standard error of the mean as a bar plot. 

    3-summary_plot - plots summary heatmaps, which includes information on the proportion of time below a threshold, residence times and transition probabilities.
