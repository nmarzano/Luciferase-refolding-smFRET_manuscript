from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import shutil
import os

output_folder = "Experiment_X-description/raw_data/"  ### Change for each experiment

input_folders = ["directory_1/",
"directory_2/",
"directory_3/",
"directory_4/"
]

def move_folders(input_folders, filetype, output_folder):
    for folder in input_folders:
        new_folder = folder.split("/")[-2]
        if not os.path.exists(f"{output_folder}{new_folder}/"):
            os.makedirs(f"{output_folder}{new_folder}/")
        filelist = [filename for filename in os.listdir(folder) if filetype in filename] # create list of files
        for filename in filelist:
            shutil.copyfile(f"{folder}{filename}", f"{output_folder}{new_folder}/{filename}")        # if need to change filenames could have a function defined above that changes it
move_folders(input_folders, ".dat", output_folder)



############# Data for TDP and subsequent analysis
containing_folders = ["directory_1/",
"directory_2/",
"directory_3/",
"directory_4/"
]

def move_file(containing_folders, filetype, output_folder):
    for files in containing_folders:
        name = files.split("/")[-2]
        filelist = [filename for filename in os.listdir(files) if filetype in filename] # create list of files
        for filename in filelist:
            shutil.copyfile(f"{files}{filename}", f"{output_folder}{name}")        # if need to change filenames could have a function defined above that changes it
move_file(containing_folders, ".dat", output_folder)