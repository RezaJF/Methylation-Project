#!/usr/bin/env python3

import os
import pandas as pd
import argparse

# parser=argparse.ArgumentParser(description="05_post_processing_aux.py")
# parser.add_argument("-i", "--input", type=str, required=True, help="Specifiy the path to input DMR directory containing DMR files")
# # parser.add_argument("-o", "--output", type=str, required=True, help="Specifiy the path where the results should be written")
# args=parser.parse_args()

# Set the directory where the DMR files are located
dir_path = '/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/04_DMW_variance_noise/DMR'
# dir_path = os.getcwd()


# dir_path = args.input
# dir_out = args.output


# Create an empty list to store the dataframes
dataframes = []

# Iterate through all text files in the DMR folder
for file_name in sorted(os.listdir(dir_path)):
    if file_name.endswith('.txt'):
        file_path = os.path.join(dir_path, file_name)
        df = pd.read_csv(file_path, delimiter='\t')
        dataframes.append(df)

# Create a list to store the variances & DMWs
variances = []
DMWs = []

# Iterate through the dataframes and calculate the variance of the 'diff' column
for df in dataframes:
    variance = df['diff'].var()
    variances.append(variance)
    
    # Calculate percentage of differentially methylated windows (Bonferroni significant level)
    DMW = len(df.loc[df['p'] <= (0.05/len(df))])/len(df)*100
    DMWs.append(DMW)

# Write the output to a tab-delimited text file
filenames = []
filenames.extend(sorted(os.listdir(dir_path)))

with open('/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/04_DMW_variance_noise/IR-TCells_variance_DMWs.txt', 'w') as f:
    for file_name, variance, DMWs in zip(filenames, variances, DMWs):
        f.write(f'{file_name}\t{variance}\t{DMWs}\n')

