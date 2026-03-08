# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: BeatAML
#     language: python
#     name: beataml
# ---

# %%
import pandas as pd
import numpy as np

final_ranked_protein_drug_df = pd.read_csv("../Results/Final_Ranking/ranking_with_drug_names/out.csv",header='infer',sep="|")
gene_names = final_ranked_protein_drug_df.iloc[:,0]
final_ranked_protein_drug_df = final_ranked_protein_drug_df.iloc[:,1:]
final_ranked_protein_drug_df.index = gene_names
drug_names = final_ranked_protein_drug_df.columns.to_list()
print(len(drug_names))

# %%
#Load the drugs from Nigel's dataset
drugbank_df = pd.read_csv("../Data/drug_protein_interactions.csv",header='infer',sep='|')
nigel_drugs_df = pd.read_csv("../Data/Nigel_Compounds.csv",header="infer",sep=",")
nigel_compounds = nigel_drugs_df["Smiles"].tolist()
drug_smiles = []
for drug_name in drug_names:
    drug_smile = drugbank_df[drugbank_df["name"]==drug_name]["smiles"].values[0]
    drug_smiles.append(drug_smile)

nigel_drugs_df

# %%
#Calculate tanimoto coefficients between nigel_compounds and our drugs
from rdkit import Chem
from rdkit import DataStructs
def compute_tanimoto_coefficients(smiles_list1, smiles_list2):
    """
    Compute Tanimoto coefficients for a list of SMILES strings.
    
    Args:
    smiles_list (list): List of SMILES strings.
    
    Returns:
    pd.DataFrame: DataFrame containing Tanimoto coefficients.
    """
    mols_1 = [Chem.MolFromSmiles(smi) for smi in smiles_list1]
    fps_1 = [Chem.RDKFingerprint(mol) for mol in mols_1]

    mols_2 = [Chem.MolFromSmiles(smi) for smi in smiles_list2]
    fps_2 = [Chem.RDKFingerprint(mol) for mol in mols_2]
    
    tanimoto_matrix = np.zeros((len(fps_1), len(fps_2)))
    for i in range(len(fps_1)):
        tanimoto_matrix[i,:] = DataStructs.cDataStructs.BulkTanimotoSimilarity(fps_1[i], fps_2)
                
    return pd.DataFrame(tanimoto_matrix, index=smiles_list1, columns=smiles_list2)


# %%
def clean_and_validate_smiles(smiles):
    """Completely clean and validate SMILES, removing all problematic patterns"""
    if not isinstance(smiles, str) or len(smiles) == 0:
        return None
    
    # List of all problematic patterns we've seen
    bad_patterns = [
        '[R]', '[R1]', '[R2]', '[R3]', '[R4]', '[R5]', 
        "[R']", '[R"]', 'R1', 'R2', 'R3', 'R4', 'R5',
        # Additional patterns that cause issues
        '([R])', '([R1])', '([R2])', 
    ]
    
    # Check for any bad patterns
    for pattern in bad_patterns:
        if pattern in smiles:
            return None
    
    # Additional check: if it contains ] followed by [ without valid atoms, likely polymer notation
    if '][' in smiles and any(x in smiles for x in ['[R', 'R]']):
        return None
    
    # Try to parse with RDKit if available
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, canonical=True)
        else:
            return None
    except:
        return None
    
    # If RDKit not available, return cleaned SMILES
    return smiles


# %%
#Get the list of drug smiles and remove indices where molecule is None
rev_drug_smiles = [clean_and_validate_smiles(smi) for smi in drug_smiles]
indices = [i for i, val in enumerate(rev_drug_smiles) if val is None]

#Final list of drug smiles and corresponding 1-to-1 drug names
final_drug_smiles = [smiles for i, smiles in enumerate(rev_drug_smiles) if i not in indices]
final_drug_names = [name for i, name in enumerate(drug_names) if i not in indices]

#Now calculate Tanimoto coefficient
tc_df = compute_tanimoto_coefficients(nigel_compounds, final_drug_smiles)

#Save the similarity matrix
tc_df.columns = final_drug_names
tc_df.index = nigel_drugs_df["Compounds Name"].tolist()
#tc_df.to_csv("../Results/Nigel_Analysis/Tanimoto_Similarity_Matrix.csv", index=True)

# %%
def get_top_n_indices(row, threshold=0.5, n=50):
    """
    Returns the column indices corresponding to the top n maximum values in a row.
    """
    # # Sort the values in descending order and get their original indices
    # sorted_indices = row.argsort()[::-1]
    
    # # Return the top n indices
    # return sorted_indices[:n].tolist()
    filtered_row = row if threshold is None else row[row > threshold]
    if filtered_row.empty:
        return []
    n = min(n, filtered_row.size)
    top_n = filtered_row.nlargest(n)
    return top_n.index.tolist()

# Apply the function to each row
topk_tc = 10
df = pd.DataFrame()
df['top_indices']= tc_df.apply(get_top_n_indices, n=topk_tc, threshold=0.5, axis=1)
df['Group'] = [1,1,1,2,2,1,3,1,1,1,4,4,4]
df['Group2'] = [1,1,1,2,3,3,5,6,3,6,4,4,4]
print(df)
df.to_csv("../Results/Nigel_Analysis/Top_"+str(topk_tc)+"_drug_index_per_compound.csv",index=True)

#Get the top drugs for each group
all_drugs_of_interest_grp_1 = df[df['Group2']==1]['top_indices'].tolist()
all_drugs_of_interest_grp_3 = df[df['Group2']==3]['top_indices'].tolist()
all_drugs_of_interest_grp_4 = df[df['Group2']==4]['top_indices'].tolist()
all_drugs_of_interest_grp_6 = df[df['Group2']==6]['top_indices'].tolist()

#FInd common drugs for each group
intersection_result_group_1 = list(set.intersection(*map(set, all_drugs_of_interest_grp_1)))
intersection_result_group_3 = list(set.intersection(*map(set, all_drugs_of_interest_grp_3)))
intersection_result_group_4 = list(set.intersection(*map(set, all_drugs_of_interest_grp_4)))
intersection_result_group_6 = list(set.intersection(*map(set, all_drugs_of_interest_grp_6)))


#Find corresponding drug names
drug_names_of_interest_group_1 = [drug_name for drug_name in intersection_result_group_1]
drug_names_of_interest_group_3 = [drug_name for drug_name in intersection_result_group_3]
drug_names_of_interest_group_4 = [drug_name for drug_name in intersection_result_group_4]
drug_names_of_interest_group_6 = [drug_name for drug_name in intersection_result_group_6]


#Create the database 
drugs_group_1 = df[df["Group2"]==1].index.tolist()
drugs_group_3 = df[df["Group2"]==3].index.tolist()
drugs_group_4 = df[df["Group2"]==4].index.tolist()
drugs_group_6 = df[df["Group2"]==6].index.tolist()

final_nigel_drug_similarity_df = pd.DataFrame({'nigel_drugs':[drugs_group_1,drugs_group_3,drugs_group_4,drugs_group_6], 
                                               'drugbank_top_drugs':[drug_names_of_interest_group_1,
                                                                     drug_names_of_interest_group_3,
                                                                     drug_names_of_interest_group_4,
                                                                     drug_names_of_interest_group_6],
                                               'length_common_drugs':[len(drug_names_of_interest_group_1),
                                                                      len(drug_names_of_interest_group_3),
                                                                      len(drug_names_of_interest_group_4),
                                                                      len(drug_names_of_interest_group_6)],
                                               'Group': [1,3,4,6]})
final_nigel_drug_similarity_df.head()


# %%
def get_top_n_row_indices(col, n=50):
    """
    Returns the row indices corresponding to the top n minimum values in a col.
    """
    # Sort the values in descending order and get their original indices
    sorted_indices = col.argsort()
    # Return the top n indices
    return sorted_indices[:n].tolist()


# %%
#Get the list of protein ranking by drug id
def get_topk_in_dataframe(topk, final_ranked_protein_drug_df, final_nigel_drug_similarity_df):
    common_proteins = []
    for i in range(final_nigel_drug_similarity_df.shape[0]):
        df = pd.DataFrame()
        drugbank_drugs = final_nigel_drug_similarity_df.iloc[i,1]
        subset_final_ranked_protein_drug_df = final_ranked_protein_drug_df[drugbank_drugs]
        df = subset_final_ranked_protein_drug_df.apply(get_top_n_row_indices, n=topk, axis=0)
        all_proteins_list = []
        for i in range(df.shape[1]):
            all_proteins_list.append(df.iloc[:,i].tolist())
        intersection_result = list(set.intersection(*map(set, all_proteins_list)))
        common_proteins.append(final_ranked_protein_drug_df.index[intersection_result].tolist())
    final_nigel_drug_similarity_df['top_common_proteins'+str(topk)]=common_proteins
    return final_nigel_drug_similarity_df


# %%
final_nigel_drug_similarity_df = get_topk_in_dataframe(100, final_ranked_protein_drug_df, final_nigel_drug_similarity_df)
final_nigel_drug_similarity_df = get_topk_in_dataframe(200, final_ranked_protein_drug_df, final_nigel_drug_similarity_df)
final_nigel_drug_similarity_df = get_topk_in_dataframe(300, final_ranked_protein_drug_df, final_nigel_drug_similarity_df)

final_nigel_drug_similarity_df.to_csv("../Results/Nigel_Analysis/Top_Common_Drugs_and_Proteins_with_Similarity_from_Drugbank_v2.csv",index=None)
final_nigel_drug_similarity_df.head()

# %%
