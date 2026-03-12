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
Pathik Compound Analysis
========================
Identifies common protein drug targets for Pathik's experimental compounds by:

1. Loading Pathik's compounds (with group assignments) and DrugBank drug-protein interactions.
2. Computing pairwise Tanimoto chemical similarity between Pathik's compounds and all DrugBank drugs.
3. For each compound group, finding DrugBank analogs that consistently have high Tanimoto
   similarity (above a configurable threshold) across all compounds in the group.
4. Identifying shared protein targets in the top-ranked proteins across all analogs within each group.
5. Mapping top-hit proteins against a priority target list to produce a condensed result per group.

Results are exported to ../Results/Pathik_Analysis/.
"""

# %%
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import DataStructs

# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------
TANIMOTO_THRESHOLD = 0.3   # minimum Tanimoto similarity to consider a DrugBank drug an analog
TOPK_TC = 20               # top-k most similar DrugBank drugs per compound
TOPK_PROTEINS = [200, 500, 1000]  # protein rank thresholds for intersection analysis

OUTPUT_DIR = "../Results/Pathik_Analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# %%
# Load the protein-drug druggability ranking matrix.
# Format: rows = proteins (genes), columns = drugs, values = ranking scores (lower = better).
# The first column contains gene names and is used as the DataFrame index.
final_ranked_protein_drug_df = pd.read_csv(
    "../Results/Final_Ranking/ranking_with_drug_names/out.csv",
    header="infer",
    sep="|"
)
gene_names = final_ranked_protein_drug_df.iloc[:, 0]
final_ranked_protein_drug_df = final_ranked_protein_drug_df.iloc[:, 1:]
final_ranked_protein_drug_df.index = gene_names

# Extract the list of drug names from column headers for downstream SMILES lookup
drug_names = final_ranked_protein_drug_df.columns.to_list()
print(f"Ranked matrix loaded: {final_ranked_protein_drug_df.shape[0]} proteins x {len(drug_names)} drugs")

# %%
# Load DrugBank drug-protein interaction data (contains drug names and SMILES strings)
drugbank_df = pd.read_csv("../Data/drug_protein_interactions.csv", header="infer", sep="|")

# Load Pathik's experimental compounds (compound names, SMILES, and group assignments)
pathik_drugs_df = pd.read_csv(
    "../Results/Pathik_Analysis/Pathik_Drugs_To_Test.csv",
    header="infer",
    sep=","
)
# Strip any accidental whitespace from column names
pathik_drugs_df.columns = pathik_drugs_df.columns.str.strip()
pathik_drugs_df["Compounds"] = pathik_drugs_df["Compounds"].str.strip()

print(f"Pathik compounds loaded: {len(pathik_drugs_df)} compounds across groups "
      f"{sorted(pathik_drugs_df['Groups'].unique())}")
pathik_drugs_df.head()

# %%
# Retrieve the SMILES string for each drug in the ranked matrix by looking it up in DrugBank
drug_smiles = []
for drug_name in drug_names:
    match = drugbank_df[drugbank_df["name"] == drug_name]["smiles"].values
    drug_smiles.append(match[0] if len(match) > 0 else None)


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
        '([R])', '([R1])', '([R2])',
    ]
    for pattern in bad_patterns:
        if pattern in smiles:
            return None

    if '][' in smiles and any(x in smiles for x in ['[R', 'R]']):
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, canonical=True)
        else:
            return None
    except Exception:
        return None


# %%
def compute_tanimoto_coefficients(smiles_list1, smiles_list2):
    """
    Compute pairwise Tanimoto similarity coefficients between two sets of compounds.

    Uses RDKit topological (RDK) fingerprints for molecular representation.
    Tanimoto similarity ranges from 0 (no overlap) to 1 (identical fingerprints).

    Args:
        smiles_list1 (list of str): SMILES strings for query compounds (Pathik's compounds).
        smiles_list2 (list of str): SMILES strings for reference compounds (DrugBank drugs).

    Returns:
        pd.DataFrame: Matrix of shape (len(smiles_list1), len(smiles_list2)) containing
                      Tanimoto similarity scores, indexed by smiles_list1 and
                      columned by smiles_list2.
    """
    mols_1 = [Chem.MolFromSmiles(smi) for smi in smiles_list1]
    fps_1 = [Chem.RDKFingerprint(mol) for mol in mols_1]

    mols_2 = [Chem.MolFromSmiles(smi) for smi in smiles_list2]
    fps_2 = [Chem.RDKFingerprint(mol) for mol in mols_2]

    tanimoto_matrix = np.zeros((len(fps_1), len(fps_2)))
    for i in range(len(fps_1)):
        tanimoto_matrix[i, :] = DataStructs.cDataStructs.BulkTanimotoSimilarity(fps_1[i], fps_2)

    return pd.DataFrame(tanimoto_matrix, index=smiles_list1, columns=smiles_list2)


# %%
# Validate and canonicalize all DrugBank drug SMILES; drop invalid entries
rev_drug_smiles = [clean_and_validate_smiles(smi) for smi in drug_smiles]
invalid_indices = {i for i, val in enumerate(rev_drug_smiles) if val is None}

final_drug_smiles = [smi for i, smi in enumerate(rev_drug_smiles) if i not in invalid_indices]
final_drug_names = [name for i, name in enumerate(drug_names) if i not in invalid_indices]

print(f"Valid DrugBank SMILES: {len(final_drug_smiles)} / {len(drug_smiles)}")

# Compute the Tanimoto similarity matrix: rows = Pathik's compounds, columns = DrugBank drugs
tc_df = compute_tanimoto_coefficients(pathik_drugs_df["SMILES"].tolist(), final_drug_smiles)
tc_df.columns = final_drug_names
tc_df.index = pathik_drugs_df["Compounds"].tolist()
tc_df.to_csv(os.path.join(OUTPUT_DIR, "Tanimoto_Similarity_Matrix.csv"), index=True)
print("Tanimoto similarity matrix saved.")


# %%
def get_top_n_similar_drugs(row, threshold=TANIMOTO_THRESHOLD, n=TOPK_TC):
    """
    Return the names of the top-n most similar DrugBank drugs for one Pathik compound.

    Filters entries above the given Tanimoto threshold, then returns up to n drug
    names with the highest similarity scores.

    Args:
        row (pd.Series): One row of the Tanimoto similarity DataFrame
                         (one Pathik compound vs all DrugBank drugs).
        threshold (float): Minimum Tanimoto similarity to consider. Default: TANIMOTO_THRESHOLD.
        n (int): Maximum number of top similar drugs to return. Default: TOPK_TC.

    Returns:
        list of str: Drug names of the top-n most similar DrugBank drugs.
    """
    filtered = row if threshold is None else row[row > threshold]
    if filtered.empty:
        return []
    top_n = filtered.nlargest(min(n, filtered.size))
    return top_n.index.tolist()


# %%
# Build per-compound top-similar DrugBank drug lists
tc_result_df = pd.DataFrame()
tc_result_df["top_similar_drugs"] = tc_df.apply(
    get_top_n_similar_drugs, threshold=TANIMOTO_THRESHOLD, n=TOPK_TC, axis=1
)
tc_result_df["Group"] = pathik_drugs_df["Groups"].values
tc_result_df.index = pathik_drugs_df["Compounds"].tolist()
tc_result_df.to_csv(
    os.path.join(OUTPUT_DIR, f"Top_{TOPK_TC}_similar_DrugBank_drugs_per_compound.csv"),
    index=True
)
print(tc_result_df)


# %%
# For each group, find DrugBank drugs consistently similar to ALL compounds in that group
# (set intersection of top-similar drugs across all compounds within the group)
groups = sorted(pathik_drugs_df["Groups"].unique())
group_summary_rows = []

for grp in groups:
    grp_mask = tc_result_df["Group"] == grp
    compounds_in_group = tc_result_df[grp_mask].index.tolist()
    drug_lists = tc_result_df[grp_mask]["top_similar_drugs"].tolist()

    # Only intersect if every compound has at least one similar drug
    non_empty = [d for d in drug_lists if len(d) > 0]
    if len(non_empty) == 0:
        common_drugs = []
    elif len(non_empty) == 1:
        common_drugs = non_empty[0]
    else:
        common_drugs = list(set.intersection(*map(set, non_empty)))

    group_summary_rows.append({
        "Group": grp,
        "pathik_compounds": compounds_in_group,
        "n_compounds": len(compounds_in_group),
        "drugbank_common_drugs": common_drugs,
        "n_common_drugs": len(common_drugs),
    })
    print(f"Group {grp}: {len(compounds_in_group)} compounds → "
          f"{len(common_drugs)} common DrugBank analogs (threshold={TANIMOTO_THRESHOLD})")

group_summary_df = pd.DataFrame(group_summary_rows)
group_summary_df.to_csv(
    os.path.join(OUTPUT_DIR, "Group_Common_DrugBank_Drugs.csv"),
    index=False
)
group_summary_df.head()


# %%
def get_top_n_row_indices(col, n=50):
    """
    Return the row positions of the top-n highest-ranked proteins for a given drug.

    Proteins are ranked by ascending score (lower score = higher rank = better druggability).

    Args:
        col (pd.Series): One column of the protein-drug ranking DataFrame.
        n (int): Number of top-ranked proteins to return. Defaults to 50.

    Returns:
        list of int: Positional row indices of the top-n ranked proteins.
    """
    sorted_indices = col.argsort()
    return sorted_indices[:n].tolist()


def get_topk_common_proteins(topk, final_ranked_protein_drug_df, group_summary_df):
    """
    Find common top-k ranked proteins across all DrugBank analogs for each Pathik compound group.

    For each compound group:
    1. Retrieves the DrugBank drugs identified as common analogs.
    2. For each analog, finds its top-k highest-ranked proteins (lowest scores).
    3. Takes the intersection — proteins ranking in the top-k for ALL analogs.

    Results are added as a new column 'top_common_proteins{topk}' in group_summary_df.

    Args:
        topk (int): Number of top-ranked proteins to consider per drug.
        final_ranked_protein_drug_df (pd.DataFrame): Protein-drug ranking matrix (proteins x drugs).
        group_summary_df (pd.DataFrame): Summary DataFrame with 'drugbank_common_drugs' column.

    Returns:
        pd.DataFrame: Updated group_summary_df with 'top_common_proteins{topk}' column.
    """
    common_proteins_per_group = []

    for _, row in group_summary_df.iterrows():
        drugbank_drugs = row["drugbank_common_drugs"]

        # Filter to drugs that exist in the ranking matrix
        available_drugs = [d for d in drugbank_drugs if d in final_ranked_protein_drug_df.columns]

        if len(available_drugs) == 0:
            common_proteins_per_group.append([])
            continue

        subset_df = final_ranked_protein_drug_df[available_drugs]
        top_indices_df = subset_df.apply(get_top_n_row_indices, n=topk, axis=0)

        all_protein_index_lists = [top_indices_df.iloc[:, j].tolist()
                                   for j in range(top_indices_df.shape[1])]

        if len(all_protein_index_lists) == 1:
            intersection = all_protein_index_lists[0]
        else:
            intersection = list(set.intersection(*map(set, all_protein_index_lists)))

        common_proteins_per_group.append(
            final_ranked_protein_drug_df.index[intersection].tolist()
        )

    col_name = f"top_common_proteins{topk}"
    group_summary_df[col_name] = common_proteins_per_group
    return group_summary_df


# %%
# Run protein target analysis at increasing top-k thresholds
for topk in TOPK_PROTEINS:
    group_summary_df = get_topk_common_proteins(topk, final_ranked_protein_drug_df, group_summary_df)
    for _, row in group_summary_df.iterrows():
        n_proteins = len(row[f"top_common_proteins{topk}"])
        print(f"  Group {row['Group']} | top-{topk}: {n_proteins} common proteins")

group_summary_df.to_csv(
    os.path.join(OUTPUT_DIR, "Group_Common_Proteins_All_Thresholds.csv"),
    index=False
)
print("\nGroup summary with protein targets saved.")
group_summary_df.head()


# %%
# Load priority target list
priority_targets_df = pd.read_csv(
    "../Results/Pathik_Analysis/Pathik_Targets_To_Priotize.csv",
    header="infer",
    sep="\t"
)
priority_targets_df.columns = priority_targets_df.columns.str.strip()
priority_protein_set = set(priority_targets_df["Protein Name"].str.strip().tolist())

print(f"Priority proteins loaded: {len(priority_protein_set)}")
print(priority_targets_df.head())


# %%
# For each group and each top-k threshold, intersect common proteins with priority targets
# and build a condensed mapping table
condensed_rows = []

for _, row in group_summary_df.iterrows():
    grp = row["Group"]
    compounds = row["pathik_compounds"]

    for topk in TOPK_PROTEINS:
        col_name = f"top_common_proteins{topk}"
        common_proteins = row[col_name]

        # Intersect with priority targets
        priority_hits = [p for p in common_proteins if p in priority_protein_set]

        # Enrich with category information from priority targets table
        hit_details = priority_targets_df[
            priority_targets_df["Protein Name"].str.strip().isin(priority_hits)
        ][["Target Type", "Category", "Protein Name"]].drop_duplicates()

        if hit_details.empty:
            condensed_rows.append({
                "Group": grp,
                "Pathik_Compounds": "; ".join(compounds),
                "Top_K_Threshold": topk,
                "Target_Type": None,
                "Category": None,
                "Priority_Protein": None,
                "N_Common_DrugBank_Drugs": row["n_common_drugs"],
            })
        else:
            for _, hit_row in hit_details.iterrows():
                condensed_rows.append({
                    "Group": grp,
                    "Pathik_Compounds": "; ".join(compounds),
                    "Top_K_Threshold": topk,
                    "Target_Type": hit_row["Target Type"],
                    "Category": hit_row["Category"],
                    "Priority_Protein": hit_row["Protein Name"],
                    "N_Common_DrugBank_Drugs": row["n_common_drugs"],
                })

condensed_df = pd.DataFrame(condensed_rows)
condensed_df.to_csv(
    os.path.join(OUTPUT_DIR, "Priority_Target_Hits_Per_Group.csv"),
    index=False
)
print(f"\nCondensed priority target mapping saved: {len(condensed_df)} rows")
condensed_df.head(20)


# %%
# Summary printout: hits per group at each top-k threshold
for topk in TOPK_PROTEINS:
    print(f"\n=== Priority Target Hits (top-k={topk}) ===")
    for grp in groups:
        subset = condensed_df[(condensed_df["Group"] == grp) & (condensed_df["Top_K_Threshold"] == topk)]
        n_compounds = (pathik_drugs_df["Groups"] == grp).sum()
        print(f"\nGroup {grp} ({n_compounds} compounds):")
        hits = subset.dropna(subset=["Priority_Protein"])
        if hits.empty:
            print("  No priority targets found at this threshold.")
        else:
            for _, r in hits.iterrows():
                print(f"  [{r['Target_Type']} | {r['Category']}] {r['Priority_Protein']}")
