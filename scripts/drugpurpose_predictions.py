# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.6
#   kernelspec:
#     display_name: ProtSol
#     language: python
#     name: protsol
# ---

import os
import pandas as pd
from DeepPurpose import utils
from DeepPurpose import DTI as models
import numpy as np
import torch
import time
import gc

# +
#Make the function to get dti virtual screening
common_path = '/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/Druggability_Score/'

def DTI_pred(drugs, targets, model_name, drug_names, target_names):
    model = models.model_pretrained(common_path+"Models/pretrained_models/"+model_name)
    model.config['batch_size']=8192
    print(model.config)

    #Get predictions using virtual screening
    y_pred = models.virtual_screening(drugs, targets, model, drug_names = drug_names, target_names = target_names, result_folder = common_path+"Results/"+model_name, convert_y = False, output_num_max = 10, verbose = True)
    return(y_pred)


# +
#Load the drug SMILES dataset and protein sequence dataset
drug_target_interactions_df = pd.read_csv("../Data/drug_protein_interactions.csv",header="infer",sep="|")

N_drugs = drug_target_interactions_df.shape[0]

#Get list of drug smiles and drug ids
drug_smiles = drug_target_interactions_df["smiles"].tolist()
drug_names = drug_target_interactions_df["drug_id"].tolist()

filtered_alphafold_df = pd.read_csv("../Data/Filtered_AlphaFold_UniProtKB_with_Sequence.csv",header="infer",sep=",")
filtered_alphafold_df.head()

duplicate_sequence_df = filtered_alphafold_df.loc[filtered_alphafold_df.groupby('Protein_Sequence').UniProtId.transform('nunique').gt(1)]
duplicate_sequence_df.to_csv("../Data/Filtered_AlphaFold_UniProtKB_with_Duplicate_Sequences.csv",index=False,sep=",")

rev_filtered_alphafold_df = filtered_alphafold_df.drop_duplicates(subset=['Protein_Sequence'])
print(rev_filtered_alphafold_df.shape)

N_proteins = rev_filtered_alphafold_df.shape[0]
protein_seqs = rev_filtered_alphafold_df["Protein_Sequence"].tolist()
protein_names = rev_filtered_alphafold_df["UniProtId"].tolist()

revised_protein_seqs = protein_seqs * N_drugs
revised_protein_names = protein_names * N_drugs
print(len(revised_protein_seqs))
# -

#Make the large list for smiles and drug names
revised_drug_smiles, revised_drug_names=[],[]
for i in range(N_drugs):
    for j in range(N_proteins):
        revised_drug_smiles.append(drug_smiles[i])
        revised_drug_names.append(drug_names[i])

# +
model_names = ["CNN_CNN_BindingDB_IC50","MPNN_CNN_DAVIS","Morgan_CNN_BindingDB_IC50","Morgan_AAC_BindingDB_IC50","Daylight_AAC_BindingDB_IC50","CNN_CNN_DAVIS",
               "CNN_CNN_BindingDB","Morgan_CNN_BindingDB","Morgan_CNN_DAVIS","MPNN_CNN_BindingDB","MPNN_CNN_KIBA","Transformer_CNN_BindingDB",
               "Daylight_AAC_DAVIS","Daylight_AAC_KIBA","Daylight_AAC_BindingDB","Morgan_AAC_BindingDB","Morgan_AAC_KIBA","Morgan_AAC_DAVIS"]

#Already ran code for CNN_CNN_BindingDB
for i in range(15,len(model_names)):
    model_name = model_names[i]
    start_time = time.time()
    predictions = DTI_pred(revised_drug_smiles, revised_protein_seqs, model_name, revised_drug_names, revised_protein_names)
    finish_time = time.time()-start_time

    #Move the results accordingly
    os.system("mv ../Results/"+model_name+"/virtual_screening.txt ../Results/"+model_name+"/"+model_name+"_virtual_screening.csv")
    os.system("mv ../Results/"+model_name+" ../../General_Data/virtual_screenings/")
    torch.cuda.empty_cache()
    gc.collect()

# -


