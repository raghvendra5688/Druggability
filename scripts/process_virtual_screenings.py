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

#Import Modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gc
import re
# %matplotlib inline
import scipy.stats as ss
import ranky as rk

# +
#List of all models
model_names = ["CNN_CNN_BindingDB_IC50","MPNN_CNN_DAVIS","Morgan_CNN_BindingDB_IC50","Morgan_AAC_BindingDB_IC50","Daylight_AAC_BindingDB_IC50","CNN_CNN_DAVIS",
               "CNN_CNN_BindingDB","Morgan_CNN_BindingDB","Morgan_CNN_DAVIS","MPNN_CNN_BindingDB","MPNN_CNN_KIBA","Transformer_CNN_BindingDB",
               "Daylight_AAC_DAVIS","Daylight_AAC_KIBA","Daylight_AAC_BindingDB","Morgan_AAC_BindingDB","Morgan_AAC_KIBA","Morgan_AAC_DAVIS"]

all_model_df = {}
for model_name in model_names:
    #Read the big data
    big_df = pd.read_csv("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/General_Data/virtual_screenings/"+model_name+"/"+model_name+"_virtual_screening.csv", header="infer",sep="|")
    #Remove auxillary columns
    big_df = big_df.drop(columns = ["Unnamed: 0","Unnamed: 5"],axis=1)
    big_df.columns = ["Rank","Drug Name","Target Name","Binding Score"]

    #Remove the rankings
    big_df = big_df.drop(columns = ["Rank"],axis=1)

    #Remove white spaces from drug name
    big_df["Drug Name"] = big_df["Drug Name"].str.replace(r'\s+', '', regex=True)

    #Remove white spaces from target name
    big_df["Target Name"]=big_df["Target Name"].str.replace(r'\s+', '', regex=True)

    #Remove all NAs
    big_df = big_df.dropna(how="all")
    big_df["Binding Score"]=big_df["Binding Score"].astype(float)

    ##If after removing NAs still NAs remain as string, then replace with 0
    #vals = np.where(big_df["Binding Score"].isnull())
    #if (len(vals)>0):
    #    big_df.loc[vals[0],2]=0   #3rd column is Binding Score
    
    #Convert dataframe to matrix
    pivoted_big_df = big_df.pivot(index="Target Name",columns="Drug Name",values="Binding Score")
    
    #Convert matrix to ranks
    pivoted_big_ranked_df = pivoted_big_df.apply(lambda x: ss.rankdata(-x,method="min",nan_policy="omit"))
    print(pivoted_big_ranked_df.shape)

    #Remove the big data file
    del big_df, pivoted_big_df
    gc.collect()

    all_model_df[model_name]=pivoted_big_ranked_df

#print(all_model_df)

# +
#Print the keys which include all drug repurpose methods
print(all_model_df.keys())

#Get list of all drug ids
all_drug_ids = pivoted_big_ranked_df.columns.tolist()

#Create an empty drug ranking dataframe which will contain aggregated ranks
final_ranked_protein_drug_df = pd.DataFrame()

k=0
#Traverse overall the drugs
for drug in all_drug_ids:
    #Create an empty dataframe
    temp_df = pd.DataFrame()
    
    #Traverse over all models
    for model_name in model_names:
        pivoted_big_ranked_df = all_model_df[model_name]
        temp_df = pd.concat([temp_df,pivoted_big_ranked_df[drug]],axis=1)
    #Set the columns as the different model names
    temp_df.columns = model_names
        
    #Get aggregated ranks for each drug
    temp_list = list(ss.rankdata(rk.borda(temp_df,method='median'),method='min',nan_policy='omit'))
    final_ranked_protein_drug_df = pd.concat([final_ranked_protein_drug_df,pd.DataFrame(temp_list)],axis=1)

    if (k%500==0):
        print("Done with "+str(k)+" drugs")
    k+=1

#Set the index of the final rank aggregated dataframe
final_ranked_protein_drug_df.columns = all_drug_ids
final_ranked_protein_drug_df.index  = pivoted_big_ranked_df.index
final_ranked_protein_drug_df.to_csv("../Results/ML_virtual_screenings/ML_virtual_screening.csv",index=True,sep=",", index_label="Target", na_rep='nan')
# -

print(final_ranked_protein_drug_df.head(5))


