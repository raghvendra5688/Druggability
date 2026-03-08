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

"""
Nigel Compound Analysis
=======================
Identifies common protein drug targets for Nigel's experimental compounds by:

1. Loading a pre-computed protein-drug druggability ranking matrix (proteins x drugs).
2. Loading Nigel's compounds and fetching their DrugBank SMILES strings.
3. Computing pairwise Tanimoto chemical similarity between Nigel's compounds and all DrugBank drugs.
4. Grouping Nigel's compounds by chemical scaffold and finding DrugBank analogs (Tanimoto > 0.5).
5. Identifying shared protein targets in the top-ranked proteins across all analogs within each group.

Results are exported to ../Results/Nigel_Analysis/.
"""

# %%
import pandas as pd
import numpy as np

# Load the protein-drug druggability ranking matrix.
# Format: rows = proteins (genes), columns = drugs, values = ranking scores.
# The first column contains gene names and is used as the DataFrame index.
final_ranked_protein_drug_df = pd.read_csv("../Results/Final_Ranking/ranking_with_drug_names/out.csv", header='infer', sep="|")
gene_names = final_ranked_protein_drug_df.iloc[:, 0]
final_ranked_protein_drug_df = final_ranked_protein_drug_df.iloc[:, 1:]
final_ranked_protein_drug_df.index = gene_names

# Extract the list of drug names from column headers for downstream SMILES lookup
drug_names = final_ranked_protein_drug_df.columns.to_list()
print(len(drug_names))

# %%
# Load DrugBank drug-protein interaction data (contains drug names and SMILES strings)
drugbank_df = pd.read_csv("../Data/drug_protein_interactions.csv", header='infer', sep='|')

# Load Nigel's experimental compounds (compound names, SMILES, and approval status)
nigel_drugs_df = pd.read_csv("../Data/Nigel_Compounds.csv", header="infer", sep=",")
nigel_compounds = nigel_drugs_df["Smiles"].tolist()

# Retrieve the SMILES string for each drug in the ranked matrix by looking it up in DrugBank
drug_smiles = []
for drug_name in drug_names:
    drug_smile = drugbank_df[drugbank_df["name"] == drug_name]["smiles"].values[0]
    drug_smiles.append(drug_smile)

nigel_drugs_df

# %%
# Calculate tanimoto coefficients between nigel_compounds and our drugs
from rdkit import Chem
from rdkit import DataStructs


def compute_tanimoto_coefficients(smiles_list1, smiles_list2):
    """
    Compute pairwise Tanimoto similarity coefficients between two sets of compounds.

    Uses RDKit topological (RDK) fingerprints for molecular representation.
    Tanimoto similarity ranges from 0 (no overlap) to 1 (identical fingerprints).

    Args:
        smiles_list1 (list of str): SMILES strings for query compounds (e.g., Nigel's compounds).
        smiles_list2 (list of str): SMILES strings for reference compounds (e.g., DrugBank drugs).

    Returns:
        pd.DataFrame: Matrix of shape (len(smiles_list1), len(smiles_list2)) containing
                      Tanimoto similarity scores, indexed by smiles_list1 and
                      columned by smiles_list2.
    """
    # Convert SMILES to RDKit molecule objects and generate topological fingerprints
    mols_1 = [Chem.MolFromSmiles(smi) for smi in smiles_list1]
    fps_1 = [Chem.RDKFingerprint(mol) for mol in mols_1]

    mols_2 = [Chem.MolFromSmiles(smi) for smi in smiles_list2]
    fps_2 = [Chem.RDKFingerprint(mol) for mol in mols_2]

    # Compute bulk Tanimoto similarity for each query compound against all references
    tanimoto_matrix = np.zeros((len(fps_1), len(fps_2)))
    for i in range(len(fps_1)):
        tanimoto_matrix[i, :] = DataStructs.cDataStructs.BulkTanimotoSimilarity(fps_1[i], fps_2)

    return pd.DataFrame(tanimoto_matrix, index=smiles_list1, columns=smiles_list2)


# %%
def clean_and_validate_smiles(smiles):
    """
    Validate and canonicalize a SMILES string, filtering out problematic entries.

    Rejects SMILES containing R-group placeholders (polymer/combinatorial chemistry
    notation such as [R], R1, etc.) that RDKit cannot parse as discrete molecules.
    Returns a canonical SMILES string for valid molecules, or None if invalid.

    Args:
        smiles (str): Input SMILES string to validate.

    Returns:
        str or None: Canonical SMILES string if valid, None if invalid or unparseable.
    """
    if not isinstance(smiles, str) or len(smiles) == 0:
        return None

    # Reject SMILES containing R-group placeholders (polymer/combinatorial chemistry notation)
    bad_patterns = [
        '[R]', '[R1]', '[R2]', '[R3]', '[R4]', '[R5]',
        "[R']", '[R"]', 'R1', 'R2', 'R3', 'R4', 'R5',
        # Bracketed forms also cause parsing errors
        '([R])', '([R1])', '([R2])',
    ]

    for pattern in bad_patterns:
        if pattern in smiles:
            return None

    # Additional check: bracket-pair notation combined with R-groups signals polymer SMILES
    if '][' in smiles and any(x in smiles for x in ['[R', 'R]']):
        return None

    # Attempt RDKit parsing; return canonical SMILES if successful
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, canonical=True)
        else:
            return None
    except:
        return None

    # Fallback (unreachable after try/except): return original SMILES
    return smiles


# %%
# Validate and canonicalize all DrugBank drug SMILES; track indices that fail
rev_drug_smiles = [clean_and_validate_smiles(smi) for smi in drug_smiles]
indices = [i for i, val in enumerate(rev_drug_smiles) if val is None]

# Build aligned lists of valid SMILES and drug names (invalid entries dropped)
final_drug_smiles = [smiles for i, smiles in enumerate(rev_drug_smiles) if i not in indices]
final_drug_names = [name for i, name in enumerate(drug_names) if i not in indices]

# Compute the Tanimoto similarity matrix: rows = Nigel's compounds, columns = DrugBank drugs
tc_df = compute_tanimoto_coefficients(nigel_compounds, final_drug_smiles)

# Replace SMILES indices with human-readable names for easier interpretation
tc_df.columns = final_drug_names
tc_df.index = nigel_drugs_df["Compounds Name"].tolist()
# tc_df.to_csv("../Results/Nigel_Analysis/Tanimoto_Similarity_Matrix.csv", index=True)

# %%
def get_top_n_indices(row, threshold=0.5, n=50):
    """
    Return the names of the top-n most similar DrugBank drugs for one Nigel compound.

    Filters entries above the given Tanimoto threshold, then returns up to n drug
    names with the highest similarity scores.

    Args:
        row (pd.Series): One row of the Tanimoto similarity DataFrame
                         (one Nigel compound vs all DrugBank drugs).
        threshold (float, optional): Minimum Tanimoto similarity to consider. Defaults to 0.5.
        n (int, optional): Maximum number of top similar drugs to return. Defaults to 50.

    Returns:
        list of str: Drug names (column labels) of the top-n most similar DrugBank drugs.
    """
    # Filter out drugs below the similarity threshold
    filtered_row = row if threshold is None else row[row > threshold]
    if filtered_row.empty:
        return []
    n = min(n, filtered_row.size)
    top_n = filtered_row.nlargest(n)
    return top_n.index.tolist()


# Find top-10 most similar DrugBank drugs per Nigel compound (Tanimoto > 0.5)
topk_tc = 10
df = pd.DataFrame()
df['top_indices'] = tc_df.apply(get_top_n_indices, n=topk_tc, threshold=0.5, axis=1)

# Manually assign chemical scaffold groups to each Nigel compound based on structural similarity.
# Group2 is the primary grouping used for downstream analysis.
df['Group'] = [1, 1, 1, 2, 2, 1, 3, 1, 1, 1, 4, 4, 4]
df['Group2'] = [1, 1, 1, 2, 3, 3, 5, 6, 3, 6, 4, 4, 4]
print(df)
df.to_csv("../Results/Nigel_Analysis/Top_" + str(topk_tc) + "_drug_index_per_compound.csv", index=True)

# Collect top-similar DrugBank drugs for each Group2 compound group
all_drugs_of_interest_grp_1 = df[df['Group2'] == 1]['top_indices'].tolist()
all_drugs_of_interest_grp_3 = df[df['Group2'] == 3]['top_indices'].tolist()
all_drugs_of_interest_grp_4 = df[df['Group2'] == 4]['top_indices'].tolist()
all_drugs_of_interest_grp_6 = df[df['Group2'] == 6]['top_indices'].tolist()

# Find DrugBank drugs common to ALL compounds within each group (set intersection)
intersection_result_group_1 = list(set.intersection(*map(set, all_drugs_of_interest_grp_1)))
intersection_result_group_3 = list(set.intersection(*map(set, all_drugs_of_interest_grp_3)))
intersection_result_group_4 = list(set.intersection(*map(set, all_drugs_of_interest_grp_4)))
intersection_result_group_6 = list(set.intersection(*map(set, all_drugs_of_interest_grp_6)))

# Collect drug names per group (already names, kept as list comprehension for clarity)
drug_names_of_interest_group_1 = [drug_name for drug_name in intersection_result_group_1]
drug_names_of_interest_group_3 = [drug_name for drug_name in intersection_result_group_3]
drug_names_of_interest_group_4 = [drug_name for drug_name in intersection_result_group_4]
drug_names_of_interest_group_6 = [drug_name for drug_name in intersection_result_group_6]

# Collect compound names per group for the summary DataFrame
drugs_group_1 = df[df["Group2"] == 1].index.tolist()
drugs_group_3 = df[df["Group2"] == 3].index.tolist()
drugs_group_4 = df[df["Group2"] == 4].index.tolist()
drugs_group_6 = df[df["Group2"] == 6].index.tolist()

# Build a summary DataFrame: one row per compound group with matched DrugBank drugs
final_nigel_drug_similarity_df = pd.DataFrame({
    'nigel_drugs': [drugs_group_1, drugs_group_3, drugs_group_4, drugs_group_6],
    'drugbank_top_drugs': [drug_names_of_interest_group_1, drug_names_of_interest_group_3,
                           drug_names_of_interest_group_4, drug_names_of_interest_group_6],
    'length_common_drugs': [len(drug_names_of_interest_group_1), len(drug_names_of_interest_group_3),
                            len(drug_names_of_interest_group_4), len(drug_names_of_interest_group_6)],
    'Group': [1, 3, 4, 6]
})
final_nigel_drug_similarity_df.head()


# %%
def get_top_n_row_indices(col, n=50):
    """
    Return the row positions of the top-n highest-ranked proteins for a given drug.

    Proteins are ranked by ascending score (lower score = higher rank = better druggability).
    Returns integer positional indices (not labels) into the ranking DataFrame.

    Args:
        col (pd.Series): One column of the protein-drug ranking DataFrame
                         (ranking scores for all proteins against a single drug).
        n (int, optional): Number of top-ranked proteins to return. Defaults to 50.

    Returns:
        list of int: Positional row indices of the top-n ranked proteins (ascending score order).
    """
    # argsort returns positions in ascending order: index 0 = protein with lowest (best) score
    sorted_indices = col.argsort()
    return sorted_indices[:n].tolist()


# %%
def get_topk_in_dataframe(topk, final_ranked_protein_drug_df, final_nigel_drug_similarity_df):
    """
    Find common top-k ranked proteins across all DrugBank analogs for each Nigel compound group.

    For each compound group:
    1. Retrieves the DrugBank drugs most similar to that group (pre-computed by Tanimoto).
    2. For each similar DrugBank drug, finds its top-k highest-ranked proteins.
    3. Takes the intersection — proteins that rank in the top-k for ALL similar drugs.

    Results are added as a new column 'top_common_proteins{topk}' in the summary DataFrame.

    Args:
        topk (int): Number of top-ranked proteins to consider per drug.
        final_ranked_protein_drug_df (pd.DataFrame): Protein-drug ranking matrix (proteins x drugs).
        final_nigel_drug_similarity_df (pd.DataFrame): Summary DataFrame with one row per compound
            group; must contain 'drugbank_top_drugs' column (list of DrugBank drug names).

    Returns:
        pd.DataFrame: Updated final_nigel_drug_similarity_df with a new column
                      'top_common_proteins{topk}' containing the intersected protein lists.
    """
    common_proteins = []
    for i in range(final_nigel_drug_similarity_df.shape[0]):
        # Get the DrugBank drugs most similar to this Nigel compound group
        drugbank_drugs = final_nigel_drug_similarity_df.iloc[i, 1]

        # Subset the ranking matrix to only the relevant DrugBank drugs
        subset_final_ranked_protein_drug_df = final_ranked_protein_drug_df[drugbank_drugs]

        # Get top-k protein row indices for each DrugBank drug (one column per drug)
        df = subset_final_ranked_protein_drug_df.apply(get_top_n_row_indices, n=topk, axis=0)

        # Collect each drug's top-k protein index list
        all_proteins_list = []
        for j in range(df.shape[1]):
            all_proteins_list.append(df.iloc[:, j].tolist())

        # Intersect: keep only proteins in the top-k for ALL drugs in this group
        intersection_result = list(set.intersection(*map(set, all_proteins_list)))
        common_proteins.append(final_ranked_protein_drug_df.index[intersection_result].tolist())

    final_nigel_drug_similarity_df['top_common_proteins' + str(topk)] = common_proteins
    return final_nigel_drug_similarity_df


# %%
# Run the protein target analysis at three increasing top-k thresholds.
# Larger k = more proteins considered per drug = fewer proteins in common (stricter overlap).
final_nigel_drug_similarity_df = get_topk_in_dataframe(100, final_ranked_protein_drug_df, final_nigel_drug_similarity_df)
final_nigel_drug_similarity_df = get_topk_in_dataframe(200, final_ranked_protein_drug_df, final_nigel_drug_similarity_df)
final_nigel_drug_similarity_df = get_topk_in_dataframe(300, final_ranked_protein_drug_df, final_nigel_drug_similarity_df)

# Save the final results: compound groups, their DrugBank analogs, and common protein targets
final_nigel_drug_similarity_df.to_csv(
    "../Results/Nigel_Analysis/Top_Common_Drugs_and_Proteins_with_Similarity_from_Drugbank_v2.csv",
    index=None
)
final_nigel_drug_similarity_df.head()

# %%
