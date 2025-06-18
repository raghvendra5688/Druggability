# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: ProtSol
#     language: python
#     name: protsol
# ---

# ## Script to process the results from AlphaFold domain docking of human proteins with FDA approved drugs using QuickVina2-GPU-1

#Import libraries
import os
import sys
import pandas as pd
import numpy as np
import pprint
import glob
import gc

#Get the name of all files with AF ids for which docking has been done
files = glob.glob('/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/General_Data/Chirag/QuickVina2_Results/**/**/*.csv', recursive=True)
files = list(set(files))


# ## Take a docking file and get docking scores

# +
#Process the files quickly and add as dataframes in a list
def process_file(file):
    temp_df = pd.read_csv(file)
    temp_df["ProteinId"] = temp_df["ProteinId"].str.replace("AF-","").str.replace("-v4","")
    return(temp_df)


frames = [process_file(file) for file in files]

#Combine the list of dataframes
combined_df = pd.concat(frames)

#Remove the list of dataframes
del frames
gc.collect()
# -

# ## Pivot the table with rows being drugs and columns being proteins

# +
#Convert the table from long table to wide table
final_df = combined_df.pivot_table(values="Dock_Score", index="DrugId", columns="ProteinId",dropna=False)
final_df.head()

##Save the table as a matrix
#final_df.to_csv("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/General_Data/Chirag/QuickVina2_Results/tmp_dock.csv",sep="\t",index=True)

#Remove the temporary dataframes
del combined_df
gc.collect()
# -

# ## Get all domains for one protein and select the one with best binding affinity

# +
import re
from collections import defaultdict

all_protein_domains = final_df.columns.tolist()
rev_protein_domains = [re.sub("-[0-9]+","",s) for s in all_protein_domains]

duplicates = defaultdict(list)

# iterate over positions and domains simultaneously
for i, protein_domain in enumerate(rev_protein_domains):
    # accumulate positions to the same number
    duplicates[protein_domain].append(i)

#Get the unique protein domain and their positions in column ids
result = {key: value for key, value in duplicates.items() if len(value) > 0}
print(result)


# +
#Create and empty dataframe and append a column to it
unique_protein_domains = list(result.keys())

frames = []
for protein_domain in unique_protein_domains:
    col_ids = result[protein_domain]
    temp_df = final_df.iloc[:,col_ids].copy()
    if (len(col_ids)>1):
        frames.append(temp_df.min(axis=1))
    else:
        frames.append(temp_df)

new_df = pd.concat(frames, axis=1)
new_df.columns = unique_protein_domains
new_df.head()

#Save the table as a matrix
new_df.to_csv("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/General_Data/Chirag/QuickVina2_Results/tmp_dock.csv",sep="\t",index=True)
# -


