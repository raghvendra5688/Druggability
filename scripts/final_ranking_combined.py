# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.0
#   kernelspec:
#     display_name: BeatAML
#     language: python
#     name: beataml
# ---

# %%
#Import Modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gc
import re
# # %matplotlib inline
import scipy.stats as ss
import ranky as rk

# %%
model_names = ['Docking','ML_virtual_screening','Propagation']
model_files = ['../Results/Docking_Scores/Docking_scores.csv', '../Results/ML_virtual_screenings/ML_virtual_screening.csv', '../Results/Propagation_Scores/Propagation_Scores.csv']

#Process the dataframes
all_models_df = {}
common_drug_ids = None
common_targets = None
gene_names = None
for i in range(len(model_names)):
    model_name = model_names[i]
    if "Dock" in model_name:
        docking_df = pd.read_csv(model_files[i], header='infer', sep="\t")
        docking_df = docking_df.transpose()
        docking_df.columns = docking_df.iloc[0,:]
        docking_df = docking_df.iloc[1:,:]
        docking_df = docking_df.astype(float)

        #Convert matrix to ranks
        docking_ranked_df = docking_df.apply(lambda x: ss.rankdata(-x,method="min",nan_policy='omit'))
        #print(docking_ranked_df.head())
        print(f"Docking ranked dataset shape: {docking_ranked_df.shape}")

        all_models_df[model_name]=docking_ranked_df
        common_drug_ids = docking_ranked_df.columns.tolist()
        common_target = docking_ranked_df.index.tolist()
        del docking_df, docking_ranked_df
        gc.collect()
        
    if 'ML' in model_name:
        ml_vs_df = pd.read_csv(model_files[i], header='infer', sep=",")
        ml_ranked_df = ml_vs_df
        ml_ranked_df.index = ml_vs_df.iloc[:,0]
        ml_ranked_df = ml_vs_df.iloc[:,1:]
        
        all_models_df[model_name] = ml_ranked_df
        #print(ml_ranked_df.head())
        print(f"ML ranked dataset shape: {ml_ranked_df.shape}")
        common_drug_ids = list(set(ml_ranked_df.columns.tolist()) & set(common_drug_ids))
        common_target = list(set(common_target) & set(ml_ranked_df.index.tolist()))
        del ml_ranked_df, ml_vs_df
        gc.collect()

    if 'Propagation' in model_name:
        propagation_df = pd.read_csv(model_files[i], header='infer', sep='\t')
        propagation_df.index.name = 'Target'

        propagation_ranked_df = propagation_df.apply(lambda x: ss.rankdata(-x,method="min",nan_policy='omit'))
        all_models_df[model_name] = propagation_ranked_df
        print(f"Propagation ranked dataset shape: {propagation_ranked_df.shape}")

        common_drug_ids = list(set(common_drug_ids) & set(propagation_ranked_df.columns.tolist()))
        gene_names = propagation_ranked_df.index.tolist()

        del propagation_df, propagation_ranked_df
        gc.collect()
        
#Print the common drug ids with length
print(f"Length of common drug ids: {len(common_drug_ids)}")
print(f"Length of common targets: {len(common_target)}") 

# %%
#Load the mapping file between genes and AlphaFold ids
gene_alphafold_mapping_df = pd.read_csv("../Data/StringDB_genes_AlphaFold_Mapping.csv",header='infer',sep="\t")
gene_names = list(set(gene_alphafold_mapping_df["StringDB_Gene_Name"].tolist()))
alphafold_ids = list(set(gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"].tolist()))
print(f"Alphafold ids length: {len(alphafold_ids)},total genes: {len(gene_names)}")

#Get the common targets
common_targets = list(set(common_target) & set(alphafold_ids))
#Get the gene names mapped to AlphaFold ids
subset_gene_alphafold_mapping_df = gene_alphafold_mapping_df[gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"].isin(common_targets)]
print(f"Gene names mapping to common AlphaFold ids: {subset_gene_alphafold_mapping_df.shape}")

# %%
#Merge the AlphaFold Ids with gene names
final_models_df = {}
for model_name in all_models_df.keys():
    if "Docking" in model_name:
        temp_docking_df = all_models_df[model_name]
        temp_docking_df = temp_docking_df[common_drug_ids]
        temp_docking_df = temp_docking_df.loc[common_targets]
        print(f"Common alphafold ids and drug ids shape: {temp_docking_df.shape}")

        #Map each alphafold id to genes and if alphafold id maps to more than 1 gene, repeat the rows
        rev_gene_index = []
        all_rows = []
        for target in common_targets:
            gene_name_list = subset_gene_alphafold_mapping_df[subset_gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"]==target]["StringDB_Gene_Name"].tolist()
            temp_row = temp_docking_df.loc[target].tolist()
            for gene in gene_name_list:
                rev_gene_index.append(gene)
                all_rows.append(temp_row)

        rev_docking_df = pd.DataFrame(all_rows, columns = common_drug_ids)
        rev_docking_df.index  = rev_gene_index

        #For one gene there can be more than one alphafold id, we take the docking score rank which is best for each drug for a given gene
        final_docking_result = rev_docking_df.groupby(rev_docking_df.index).min()
        print(f"Final dataframe comprising gene names x drug ids for docking results: {final_docking_result.shape}")
        print(final_docking_result.head())

        #Add the docking result with gene name x drug ids
        final_models_df[model_name] = final_docking_result
    if "ML" in model_name:
        temp_ml_vs_df = all_models_df[model_name]
        temp_ml_vs_df = temp_ml_vs_df[common_drug_ids]
        temp_ml_vs_df = temp_ml_vs_df.loc[common_targets]
        print(f"Common alphafold id and drug ids shape: {temp_ml_vs_df.shape}")

        rev_gene_index = []
        all_rows = []
        for target in common_targets:
            gene_name_list = subset_gene_alphafold_mapping_df[subset_gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"]==target]["StringDB_Gene_Name"].tolist()
            temp_row = temp_ml_vs_df.loc[target].tolist()
            for gene in gene_name_list:
                rev_gene_index.append(gene)
                all_rows.append(temp_row)

        rev_ml_vs_df = pd.DataFrame(all_rows, columns = common_drug_ids)
        rev_ml_vs_df.index  = rev_gene_index

        #For one gene there can be more than one alphafold id, we take the ML score rank which is best for each drug for a given gene
        final_ml_vs_result = rev_ml_vs_df.groupby(rev_ml_vs_df.index).min()
        print(f"Final dataframe comprising gene names x drug ids for virtual screening: {final_ml_vs_result.shape}")
        print(final_ml_vs_result.head())

        #Add the virtual screening with gene name x drug ids
        final_models_df[model_name] = final_ml_vs_result
        
    if "Propagation" in model_name:
        final_models_df[model_name] = all_models_df[model_name]


# %%
common_drugs = final_models_df['Docking'].columns.tolist()
common_gene_names = final_models_df['Docking'].index.tolist()

prop_drugs = final_models_df['Propagation'].columns.tolist()
prop_genes = final_models_df['Propagation'].index.tolist()

common_gene_names = list(set(common_gene_names) & set(prop_genes))
common_drugs = list(set(common_drugs) & set(prop_drugs))

# %%
ultimate_models_df = {}
for model_name in model_names:
    temp_df = final_models_df[model_name]
    temp_df = temp_df[common_drugs]
    temp_df = temp_df.loc[common_gene_names]
    print(f"Common genes and drug ids shape: {temp_df.shape}")
    ultimate_models_df[model_name] = temp_df
del final_models_df
gc.collect()

# %%
#Create an empty drug ranking dataframe which will contain aggregated ranks
final_ranked_protein_drug_df = pd.DataFrame()
k=0
#Traverse overall the drugs
for drug in common_drugs:
    #Create an empty dataframe
    temp_df = pd.DataFrame()
    
    #Traverse over all models
    for model_name in model_names:
        ranked_df = ultimate_models_df[model_name]
        temp_df = pd.concat([temp_df,ranked_df[drug]],axis=1)
        
    #Set the columns as the different model names
    temp_df.columns = model_names
        
    #Get aggregated ranks for each drug
    temp_list = list(ss.rankdata(rk.borda(temp_df,method='median'),method='min',nan_policy='omit'))
    final_ranked_protein_drug_df = pd.concat([final_ranked_protein_drug_df,pd.DataFrame(temp_list)],axis=1)

    if (k%500==0):
        print("Done with "+str(k)+" drugs")
    k+=1

#Set the index of the final rank aggregated dataframe
final_ranked_protein_drug_df.columns = common_drugs
final_ranked_protein_drug_df.index  = ranked_df.index
del ultimate_models_df
gc.collect()

# %%
#Replace drug id by drug names in the final ranked matrix
drugbank_df = pd.read_csv("../Data/drug_protein_interactions.csv",header='infer',sep='|')
subset_drugbank_df = drugbank_df[["drug_id","name"]]
drug_names = []
for drug_id in common_drugs:
    drug_name = subset_drugbank_df[subset_drugbank_df["drug_id"]==drug_id]["name"].values.astype(str)
    drug_names.append(drug_name[0])

final_ranked_protein_drug_df.head()
#final_ranked_protein_drug_df.to_csv("../Results/Final_Ranking/drug_protein_ranking_with_drug_ids.csv",sep=",")
#final_ranked_protein_drug_df.columns = drug_names
#final_ranked_protein_drug_df.to_csv("../Results/Final_Ranking/drug_protein_ranking_with_drug_names.csv",sep="|")

# %%
