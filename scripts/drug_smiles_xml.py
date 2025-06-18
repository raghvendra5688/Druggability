# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: feature-pHLA
#     language: python
#     name: feature-phla
# ---

import untangle
import pandas as pd
import numpy as np
import os

#takes couple of minutes
filename="../../General_Data/drug_Full_DB.xml" # DrugBank Version 6
obj=untangle.parse(filename)

#Data Frame of Small Molecule Type Drugs
df_drug_sm=pd.DataFrame(columns=["drug_id","name","cas","smiles","logP ALOGPS", "logP ChemAxon", "solubility ALOGPS", "pKa (strongest acidic)", "pKa (strongest basic)"])
df_drug_sm

# Takes around 10 minutes to run.
i=-1
#iterate over drug entries to extract information
for drug in obj.drugbank.drug:
    drug_type= str(drug["type"])
    
    # select for small molecule drugs
    if drug_type in ["small molecule", "Small Molecule", "Small molecule"]:
        i=i+1    
        
        #Get drug_id
        for id in drug.drugbank_id:
            if str(id["primary"])=="true":
                df_drug_sm.loc[i, "drug_id"]=id.cdata
        
        #Drug name
        df_drug_sm.loc[i,"name"]=drug.name.cdata
        
        #Drug CAS
        df_drug_sm.loc[i, "cas"]=drug.cas_number.cdata
        
        #Get SMILES, logP, Solubility
        #Skip drugs with no structure. ("DB00386","DB00407","DB00702","DB00785","DB00840",
        #                                            "DB00893","DB00930","DB00965", "DB01109","DB01266",
        #                                           "DB01323", "DB01341"...)
        if len(drug.calculated_properties.cdata)==0: #If there is no calculated properties
            continue
        else:
            for property in drug.calculated_properties.property:
                if property.kind.cdata == "SMILES":
                    df_drug_sm.loc[i, "smiles"]=property.value.cdata
                    
                if property.kind.cdata == "logP":
                    if property.source.cdata == "ALOGPS":
                        df_drug_sm.loc[i, "logP ALOGPS"]=property.value.cdata
                    if property.source.cdata == "ChemAxon":
                        df_drug_sm.loc[i, "logP ChemAxon"]=property.value.cdata
                
                if property.kind.cdata == "Water Solubility":
                    df_drug_sm.loc[i, "solubility ALOGPS"]=property.value.cdata
                
                if property.kind.cdata == "pKa (strongest acidic)":
                    df_drug_sm.loc[i, "pKa (strongest acidic)"]=property.value.cdata
                
                if property.kind.cdata == "pKa (strongest basic)":
                    df_drug_sm.loc[i, "pKa (strongest basic)"]=property.value.cdata


df_drug_sm.head(10)


#Drop drugs without SMILES from the dataframe
df_drug_smiles = df_drug_sm.dropna()
df_drug_smiles= df_drug_smiles.reset_index(drop=True)
print(df_drug_smiles.shape)

#write to csv
df_drug_smiles.to_csv("../Data/smiles.csv", encoding='utf-8', sep="|", index=False)


