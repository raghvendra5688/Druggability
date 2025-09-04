# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.1
#   kernelspec:
#     display_name: BeatAML
#     language: python
#     name: python3
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
model_names = ['docking','ML_virtual_screening']
model_files = ['../Results/Docking_Scores/docking_scores.csv', '../Results/ML_virtual_screenings/ML_virtual_screening_by_protein.csv']

#Process the dataframes
all_models_df = {}
common_drug_ids = None
common_targets = None
gene_names = None
for i in range(len(model_names)):
    model_name = model_names[i]
    if "dock" in model_name:
        docking_df = pd.read_csv(model_files[i], header='infer', sep="\t")
        #Set first column as index for new docking_ranked_df
        rev_docking_df = pd.DataFrame(docking_df.iloc[:,1:])
        rev_docking_df.index = docking_df.iloc[:,0] # type: ignore
        print(rev_docking_df.head())

        #Convert matrix to ranks
        docking_ranked_df = rev_docking_df.apply(lambda x: ss.rankdata(-x,method="min",nan_policy='omit'))
        del rev_docking_df, docking_df
        gc.collect()
        print(f"Docking ranked dataset shape: {docking_ranked_df.shape}")

        all_models_df[model_name]=docking_ranked_df
        common_target = docking_ranked_df.columns.tolist()
        common_drug_ids = docking_ranked_df.index.tolist()
        del docking_ranked_df
        gc.collect()
        
    if 'ML' in model_name:
        ml_vs_df = pd.read_csv(model_files[i], header='infer', sep=",")

        #Set first column as index for new ml_ranked_df
        ml_ranked_df = pd.DataFrame(ml_vs_df.iloc[:,1:])
        ml_ranked_df.index = ml_vs_df.iloc[:,0] # type: ignore
        print(ml_ranked_df.head())
        
        all_models_df[model_name] = ml_ranked_df
        #print(ml_ranked_df.head())
        print(f"ML ranked dataset shape: {ml_ranked_df.shape}")
        common_drug_ids = list(set(ml_ranked_df.index.tolist()) & set(common_drug_ids))
        common_target = list(set(common_target) & set(ml_ranked_df.columns.tolist()))
        del ml_ranked_df, ml_vs_df
        gc.collect()


#Print the common drug ids with length
print(f"Length of common drug ids: {len(common_drug_ids)}") # type: ignore
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
    if "docking" in model_name:
        temp_docking_df = all_models_df[model_name]
        temp_docking_df = temp_docking_df.loc[common_drug_ids]
        temp_docking_df = temp_docking_df[common_targets]
        print(f"Common drug ids and alphafold ids shape: {temp_docking_df.shape}")

        #Map each alphafold id to genes and if alphafold id maps to more than 1 gene, repeat the rows
        rev_gene_index = []
        all_cols = []
        for target in common_targets:
            gene_name_list = subset_gene_alphafold_mapping_df[subset_gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"]==target]["StringDB_Gene_Name"].tolist()
            temp_col = temp_docking_df[target].tolist()
            for gene in gene_name_list:
                rev_gene_index.append(gene)
                all_cols.append(temp_col)

        rev_docking_df = pd.DataFrame(all_cols)
        rev_docking_df = rev_docking_df.T
        rev_docking_df.columns = rev_gene_index
        rev_docking_df.index  = common_drug_ids # type: ignore

        #For one gene there can be more than one alphafold id, we take the docking score rank which is best for any drug for a given gene
        temp_docking_df = rev_docking_df.T
        temp_docking_df = temp_docking_df.groupby(temp_docking_df.index).min()
        final_docking_result = temp_docking_df.T
        print(f"Final dataframe comprising drug ids x gene names for docking results: {final_docking_result.shape}")
        print(final_docking_result.head())

        #Add the docking result with gene name x drug ids
        final_models_df[model_name] = final_docking_result
    if "ML" in model_name:
        temp_ml_vs_df = all_models_df[model_name]
        temp_ml_vs_df = temp_ml_vs_df.loc[common_drug_ids]
        temp_ml_vs_df = temp_ml_vs_df[common_targets]
        print(f"Common drug ids and alphafold ids shape: {temp_ml_vs_df.shape}")

        rev_gene_index = []
        all_cols = []
        for target in common_targets:
            gene_name_list = subset_gene_alphafold_mapping_df[subset_gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"]==target]["StringDB_Gene_Name"].tolist()
            temp_col = temp_ml_vs_df[target].tolist()
            for gene in gene_name_list:
                rev_gene_index.append(gene)
                all_cols.append(temp_col)

        rev_ml_vs_df = pd.DataFrame(all_cols) #, columns = common_drug_ids)
        rev_ml_vs_df = rev_ml_vs_df.T
        rev_ml_vs_df.columns = rev_gene_index
        rev_ml_vs_df.index  = common_drug_ids

        #For one gene there can be more than one alphafold id, we take the ML score rank which is best for each drug for a given gene
        temp_ml_vs_df = rev_ml_vs_df.T
        temp_ml_vs_df = temp_ml_vs_df.groupby(temp_ml_vs_df.index).min()        
        final_ml_vs_result = temp_ml_vs_df.T
        print(f"Final dataframe comprising drug ids x gene names for virtual screening: {final_ml_vs_result.shape}")
        print(final_ml_vs_result.head())

        #Add the virtual screening with gene name x drug ids
        final_models_df[model_name] = final_ml_vs_result
        del rev_ml_vs_df, temp_ml_vs_df, final_ml_vs_result
        gc.collect()


# %%
#Find the alphafold id corresponding to a given gene name
def get_alphafold_id(gene_name):
    alphafold_id = subset_gene_alphafold_mapping_df[subset_gene_alphafold_mapping_df["StringDB_Gene_Name"]==gene_name]["AlphaFold_Uniprot_Id"].tolist()[0]
    return alphafold_id
#Find the gene name corresponding to a given alphafold id
def get_gene_name(alphafold_id):
    gene_name = subset_gene_alphafold_mapping_df[subset_gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"]==alphafold_id]["StringDB_Gene_Name"].tolist()[0]
    return gene_name
zbp1_alphafold_id = get_alphafold_id("ZBP1")
print(f"ZBP1 alphafold id: {zbp1_alphafold_id}")
zbp1_gene_name = get_gene_name(zbp1_alphafold_id)
print(f"{zbp1_alphafold_id} gene name: {zbp1_gene_name}")

#Find all the drug ids corresponding to a given gene name
def get_drug_ids_for_gene(gene_name):
    drug_ids = {}
    for model_name in model_names:
        model_df = final_models_df[model_name]
        drug_ids[model_name] = model_df[gene_name].sort_values(ascending=True).index.tolist()
        del model_df
        gc.collect()
    return drug_ids
zbp1_drug_ids = get_drug_ids_for_gene("ZBP1")
print(f"ZBP1 drug ids: {zbp1_drug_ids['docking'][:5]}, {zbp1_drug_ids['ML_virtual_screening'][:5] }")

# %%
common_drugs = final_models_df['docking'].index.tolist()
common_gene_names = final_models_df['docking'].columns.tolist()

ml_drugs = final_models_df['ML_virtual_screening'].index.tolist()
ml_gene_names = final_models_df['ML_virtual_screening'].columns.tolist()

common_gene_names = list(set(common_gene_names) & set(ml_gene_names))
common_drugs = list(set(common_drugs) & set(ml_drugs))
print(f"Final common drugs length: {len(common_drugs)}, common genes length: {len(common_gene_names)}")

# %%
ultimate_models_df = {}
for model_name in model_names:
    temp_df = final_models_df[model_name]
    temp_df = temp_df.loc[common_drugs]
    temp_df = temp_df[common_gene_names]
    print(f"Common genes and drug ids shape: {temp_df.shape}")
    ultimate_models_df[model_name] = temp_df
del final_models_df
gc.collect()

# %%
#Create an empty drug ranking dataframe which will contain aggregated ranks
final_ranked_protein_drug_df = pd.DataFrame()
k=0
#Traverse overall the drugs
for gene_name in common_gene_names:
    #Create an empty dataframe
    temp_df = pd.DataFrame()
    
    #Traverse over all models
    for model_name in model_names:
        ranked_df = ultimate_models_df[model_name]
        temp_df = pd.concat([temp_df,ranked_df[gene_name]],axis=1)
        
    #Set the columns as the different model names
    temp_df.columns = model_names
        
    #Get aggregated ranks for each drug
    temp_list = list(ss.rankdata(rk.borda(temp_df,method='mean'),method='min',nan_policy='omit'))
    final_ranked_protein_drug_df = pd.concat([final_ranked_protein_drug_df,pd.DataFrame(temp_list)],axis=1)

    if (k%500==0):
        print("Done with "+str(k)+" proteins")
    k+=1

#Set the index of the final rank aggregated dataframe
final_ranked_protein_drug_df.columns = ranked_df.columns
final_ranked_protein_drug_df.index  = common_drugs
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
#final_ranked_protein_drug_df.to_csv("../Results/Final_Ranking/protein_drug_ranking_with_drug_ids.csv",sep=",")
final_ranked_protein_drug_names_df = final_ranked_protein_drug_df.copy()
del final_ranked_protein_drug_df
gc.collect()

final_ranked_protein_drug_names_df.index = drug_names
#final_ranked_protein_drug_df.to_csv("../Results/Final_Ranking/protein_drug_ranking_with_drug_names.csv",sep="|")

# %%
#Make a function to get the top 10 drug candidates for a given protein
def get_top_drug_candidates_for_protein(gene_name, top_n=10):
    if gene_name not in final_ranked_protein_drug_names_df.columns:
        print(f"Gene name {gene_name} not found in the dataframe.")
        return None
    top_drugs = final_ranked_protein_drug_names_df[gene_name].nsmallest(top_n)
    return top_drugs

#Get the top 10 drug candidates for a given protein
top_drugs = get_top_drug_candidates_for_protein("NLRP3", top_n=10)
print(top_drugs)

top_drugs.to_csv("../Results/Raj_Analysis/top_10_drugs_for_NLRPS.csv",sep="\t")

# %%
final_ranked_protein_drug_names_df["NLRP3"].nsmallest(1000)

# %%
