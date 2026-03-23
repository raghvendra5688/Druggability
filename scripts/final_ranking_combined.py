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
"""
Cell 1: Import Required Libraries
----------------------------------
pandas, numpy  : data loading and matrix manipulation
matplotlib     : plotting (available but not used in core pipeline)
gc             : manual garbage collection to manage memory on large dataframes
re             : regular expressions (available for string parsing if needed)
scipy.stats    : rankdata() used to convert raw scores to per-gene drug ranks
ranky          : imported for compatibility; rank aggregation uses mean rank method
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gc
import re
# # %matplotlib inline
import scipy.stats as ss
import ranky as rk

# %%
"""
Cell 2: Load and Preprocess Docking, ML Virtual Screening, and Propagation Data
---------------------------------------------------------------------------------
All three models are converted to a common layout and rank convention:
  layout : targets x drugs  (rows = protein targets/genes, columns = drug IDs)
  rank 1 = best target for a given drug (column-wise ranking per drug)

Docking Scores (tab-separated, drugs x targets layout):
  - First column = drug IDs, remaining columns = AlphaFold target IDs.
  - Transposed to targets x drugs after loading.
  - Raw docking scores in kcal/mol; more negative = stronger predicted binding.
  - Ranked column-wise (per drug) with rankdata(x) ascending:
      rank 1 = minimum (most negative) score = best binding target for that drug.
  - Ties: method="min"; NaN: nan_policy="omit".

ML Virtual Screening (comma-separated, targets x drugs layout):
  - First column = AlphaFold target IDs (used as index), remaining = drug IDs.
  - Already in targets x drugs layout — no transpose needed.
  - Higher ML value = stronger predicted interaction.
  - Ranked column-wise (per drug) with rankdata(-x):
      rank 1 = target with highest ML score = best predicted target for that drug.

Propagation Scores (tab-separated, genes x drugs layout):
  - Rows = StringDB gene names (auto-detected as index), Columns = drug IDs.
  - Already in genes x drugs layout — no transpose needed.
  - Higher propagation score = stronger drug-gene association.
  - Ranked column-wise (per drug) with rankdata(-x):
      rank 1 = gene with highest propagation score = best target for that drug.

At the end, common drug IDs (columns) and common AlphaFold targets (rows) are
identified across Docking and ML. Propagation gene names are aligned in Cell 6.
"""

model_names = ['Docking', 'ML_virtual_screening', 'Propagation']
model_files = [
    '../Results/Docking_Scores/Docking_scores.csv',
    '../Results/ML_virtual_screenings/ML_virtual_screening.csv',
    '../Results/Propagation_Scores/Propagation_Scores.csv'
]

all_models_df = {}    # stores the processed (ranked) dataframe for each model
common_drug_ids = None
common_target = None

for i, model_name in enumerate(model_names):

    if "Dock" in model_name:
        # Load tab-separated docking score matrix (drugs x targets).
        # First column = drug IDs, remaining columns = AlphaFold target IDs.
        docking_df = pd.read_csv(model_files[i], header='infer', sep="\t")

        # Set drug IDs as index, then transpose to targets x drugs
        rev_docking_df = pd.DataFrame(docking_df.iloc[:, 1:].astype(float))
        rev_docking_df.index = docking_df.iloc[:, 0]  # type: ignore
        docking_df_T = rev_docking_df.T                # targets x drugs
        del rev_docking_df, docking_df
        gc.collect()
        print(docking_df_T.head())

        # Rank column-wise (per drug): for each drug, rank all targets ascending.
        # rank 1 = minimum (most negative) docking score = best binding target for that drug.
        docking_ranked_df = docking_df_T.apply(
            lambda x: ss.rankdata(x, method="min", nan_policy='omit')
        )
        del docking_df_T
        gc.collect()
        print(f"Docking ranked dataset shape: {docking_ranked_df.shape}")

        all_models_df[model_name] = docking_ranked_df
        common_target = docking_ranked_df.index.tolist()    # AlphaFold target IDs (rows)
        common_drug_ids = docking_ranked_df.columns.tolist()  # drug IDs (columns)
        del docking_ranked_df
        gc.collect()

    if 'ML' in model_name:
        # Load comma-separated ML virtual screening matrix (targets x drugs).
        # First column = AlphaFold target IDs (used as index), remaining = drug IDs.
        ml_vs_df = pd.read_csv(model_files[i], header='infer', sep=",")

        # Set target IDs as index — already in targets x drugs layout
        ml_vs_df.index = ml_vs_df.iloc[:, 0]
        ml_vs_df = ml_vs_df.iloc[:, 1:].astype(float)
        print(ml_vs_df.head())

        # Rank column-wise (per drug): for each drug, rank all targets.
        # Higher ML value = better predicted interaction for that target.
        # rankdata(-x): rank 1 = target with highest ML score = best target for that drug.
        ml_ranked_df = ml_vs_df.apply(
            lambda x: ss.rankdata(-x, method="min", nan_policy='omit')
        )
        del ml_vs_df
        gc.collect()
        print(f"ML ranked dataset shape: {ml_ranked_df.shape}")

        all_models_df[model_name] = ml_ranked_df

        # Retain only drugs and targets present in both docking and ML datasets
        common_drug_ids = list(set(ml_ranked_df.columns.tolist()) & set(common_drug_ids))
        common_target = list(set(common_target) & set(ml_ranked_df.index.tolist()))
        del ml_ranked_df
        gc.collect()

    if 'Propagation' in model_name:
        # Load tab-separated propagation score matrix (genes x drugs).
        # Gene names are auto-detected as index; drug IDs are column headers.
        propagation_df = pd.read_csv(model_files[i], header='infer', sep='\t')
        print(propagation_df.head())

        # Rank column-wise (per drug): for each drug, rank all genes.
        # Higher propagation score = stronger drug-gene association = better target.
        # rankdata(-x): rank 1 = gene with highest propagation score = best target for that drug.
        propagation_ranked_df = propagation_df.apply(
            lambda x: ss.rankdata(-x.astype(float), method='min', nan_policy='omit')
        )
        del propagation_df
        gc.collect()
        print(f"Propagation ranked dataset shape: {propagation_ranked_df.shape}")

        all_models_df[model_name] = propagation_ranked_df
        del propagation_ranked_df
        gc.collect()

print(f"Length of common drug ids (Docking ∩ ML): {len(common_drug_ids)}")   # type: ignore
print(f"Length of common AlphaFold targets (Docking ∩ ML): {len(common_target)}")


# %%
"""
Cell 3: Load Gene Name to AlphaFold ID Mapping
------------------------------------------------
Docking and ML screening results use UniProt/AlphaFold IDs as column headers.
This cell loads a mapping file that links AlphaFold IDs to human-readable
StringDB gene names (e.g., P10415 -> BCL2).

The intersection of common targets (from Cell 2) and known AlphaFold IDs
defines the final set of proteins carried forward for Docking and ML.
Propagation uses StringDB gene names directly and is aligned separately in Cell 6.
"""

# Load the mapping between StringDB gene names and AlphaFold UniProt IDs
gene_alphafold_mapping_df = pd.read_csv(
    "../Data/StringDB_genes_AlphaFold_Mapping.csv", header='infer', sep="\t"
)
gene_names = list(set(gene_alphafold_mapping_df["StringDB_Gene_Name"].tolist()))
alphafold_ids = list(set(gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"].tolist()))
print(f"AlphaFold ids length: {len(alphafold_ids)}, total genes: {len(gene_names)}")

# Intersect with AlphaFold targets present in both Docking and ML datasets
common_targets = list(set(common_target) & set(alphafold_ids))

# Subset the mapping to only those AlphaFold IDs present in common targets
subset_gene_alphafold_mapping_df = gene_alphafold_mapping_df[
    gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"].isin(common_targets)
]
print(f"Gene names mapping to common AlphaFold ids: {subset_gene_alphafold_mapping_df.shape}")


# %%
"""
Cell 4: Replace AlphaFold IDs with Gene Names and Resolve Multi-Mapping
------------------------------------------------------------------------
AlphaFold IDs (row labels in Docking and ML) are replaced with human-readable
gene names. Two sources of redundancy are handled:
  1. One AlphaFold ID -> multiple gene names: the rank row is assigned to each gene.
  2. One gene name <- multiple AlphaFold IDs: for each drug, only the minimum rank
     across all AlphaFold structures for that gene is kept.
     (minimum rank = best-scoring structure for that drug, rank 1 = best)

Propagation already uses StringDB gene names so it is not passed through this step.

Output: final_models_df — dict keyed by model name. All values are DataFrames
of shape (gene_names x drugs).
"""

# Pre-build lookup dict once: AlphaFold ID -> [gene_name, ...]
af_to_genes = (
    subset_gene_alphafold_mapping_df
    .groupby("AlphaFold_Uniprot_Id")["StringDB_Gene_Name"]
    .apply(list)
    .to_dict()
)

def expand_and_resolve(model_df, common_drug_ids, common_targets, af_to_genes):
    """
    Subset to common targets/drugs, expand AlphaFold IDs to gene names,
    and resolve multi-mapping by taking the minimum rank across all AlphaFold
    structures for each gene:
      - One AlphaFold ID -> multiple genes: each gene receives the same rank row.
      - Multiple AlphaFold IDs -> one gene: element-wise minimum rank per drug
        (best-scoring structure wins for each drug independently).

    Parameters
    ----------
    model_df        : DataFrame (AlphaFold IDs x drugs) with rank values.
    common_drug_ids : list of drug IDs to retain (column subset).
    common_targets  : list of AlphaFold IDs to expand (row subset).
    af_to_genes     : dict {alphafold_id: [gene_name, ...]} pre-built from mapping.

    Returns
    -------
    DataFrame of shape (unique_gene_names x drugs).
    """
    df = model_df.loc[common_targets, common_drug_ids]
    print(f"  Subset shape (AlphaFold IDs x drugs): {df.shape}")

    gene_arrays = {}  # gene -> minimum rank array across all its AlphaFold structures

    for target in common_targets:
        row = df.loc[target].to_numpy()   # drug ranks for this target, shape (n_drugs,)

        for gene in af_to_genes.get(target, []):
            if gene in gene_arrays:
                # Keep the minimum rank per drug across all structures for this gene
                np.minimum(gene_arrays[gene], row, out=gene_arrays[gene])
            else:
                gene_arrays[gene] = row.copy()

    # gene_arrays[gene] = array of length n_drugs → build genes x drugs DataFrame
    result = pd.DataFrame(gene_arrays, index=common_drug_ids).T
    print(f"  Expanded shape (gene names x drugs): {result.shape}")
    return result


final_models_df = {}

for model_name in model_names:
    if model_name == 'Propagation':
        # Propagation already uses StringDB gene names — no AlphaFold mapping required.
        # Store as-is; drug/gene alignment with the other models is done in Cell 6.
        final_models_df[model_name] = all_models_df[model_name]
        print(f"Propagation stored directly (shape: {all_models_df[model_name].shape})")
    else:
        print(f"Processing model: {model_name}")
        result = expand_and_resolve(
            all_models_df[model_name], common_drug_ids, common_targets, af_to_genes
        )
        print(result.head())
        final_models_df[model_name] = result
        del result
        gc.collect()

del all_models_df
gc.collect()


# %%
"""
Cell 5: Utility Functions for Gene/AlphaFold ID Lookup
-------------------------------------------------------
Defines helper functions to convert between AlphaFold UniProt IDs and
StringDB gene names, and to retrieve per-model drug rankings for a gene.

NOTE: get_drug_ids_for_gene() references final_models_df, which exists only
between cells 4 and 7. It will raise NameError if called after cell 7 runs.

BCL2 (UniProt: P10415) is used as a sanity-check example.
The mapping is many-to-many (one gene -> many IDs, one ID -> many genes),
so both lookup functions return the full list of matches.
"""

def get_alphafold_ids(gene_name):
    """Return all AlphaFold UniProt IDs mapped to a given gene name."""
    ids = subset_gene_alphafold_mapping_df[
        subset_gene_alphafold_mapping_df["StringDB_Gene_Name"] == gene_name
    ]["AlphaFold_Uniprot_Id"].tolist()
    if not ids:
        print(f"WARNING: gene '{gene_name}' not found in mapping.")
    return ids

def get_gene_names(alphafold_id):
    """Return all gene names mapped to a given AlphaFold UniProt ID."""
    names = subset_gene_alphafold_mapping_df[
        subset_gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"] == alphafold_id
    ]["StringDB_Gene_Name"].tolist()
    if not names:
        print(f"WARNING: AlphaFold ID '{alphafold_id}' not found in mapping.")
    return names

def get_gene_names_for_drug(drug_id):
    """
    Return genes sorted by rank (ascending) for a given drug, for each model.
    Rank 1 = best predicted target for that drug.

    Returns
    -------
    dict of {model_name: pd.Series}
        Each Series has gene_name as index and rank value as values, sorted
        ascending by rank. Use series.index to get ordered gene names;
        use series.loc[gene_name] to look up a specific gene's rank.

    NOTE: Uses final_models_df — only valid between cells 4 and 7.
    After cell 7 runs (which deletes final_models_df), use
    final_ranked_protein_drug_df or final_ranked_protein_drug_names_df instead.
    """
    result = {}
    for model_name in model_names:
        if drug_id in final_models_df[model_name].columns:
            result[model_name] = (
                final_models_df[model_name][drug_id]
                .sort_values(ascending=True)
                .rename("rank")
            )
        else:
            print(f"WARNING: drug '{drug_id}' not found in {model_name} model.")
    return result

# Sanity check: for Venetoclax (DB11581), BCL2 should rank near the top in all models.
# Venetoclax is a well-known selective BCL2 inhibitor.
venetoclax_genes = get_gene_names_for_drug("DB11581")
for model_name, series in venetoclax_genes.items():
    print(f"\n{model_name} — top 5 genes for DB11581 (Venetoclax):")
    print(series.head())
    if 'BCL2' in series.index:
        print(f"  BCL2 rank: {series.loc['BCL2']}, position: {series.index.get_loc('BCL2') + 1}")


# %%
"""
Cell 6: Identify Final Common Drugs and Gene Names Across All Three Models
---------------------------------------------------------------------------
After gene-name mapping (Cell 4), all three models are in gene_names x drugs format
(rows = gene names, columns = drug IDs).

This cell computes the three-way intersection of drug IDs and gene names to ensure
all model matrices are perfectly aligned before rank aggregation.

Note: Propagation may cover a different set of drugs and/or genes than Docking/ML,
so the final intersection can be smaller than the pairwise Docking ∩ ML set.
"""

# Collect drugs (columns) and genes (index/rows) from all three models
docking_drugs      = final_models_df['Docking'].columns.tolist()
docking_genes      = final_models_df['Docking'].index.tolist()
ml_drugs           = final_models_df['ML_virtual_screening'].columns.tolist()
ml_genes           = final_models_df['ML_virtual_screening'].index.tolist()
prop_drugs         = final_models_df['Propagation'].columns.tolist()
prop_genes         = final_models_df['Propagation'].index.tolist()

# Three-way intersection of drugs and gene names
common_drugs      = list(set(docking_drugs) & set(ml_drugs) & set(prop_drugs))
common_gene_names = list(set(docking_genes) & set(ml_genes) & set(prop_genes))

print(f"Docking:     {len(docking_drugs)} drugs, {len(docking_genes)} genes")
print(f"ML:          {len(ml_drugs)} drugs, {len(ml_genes)} genes")
print(f"Propagation: {len(prop_drugs)} drugs, {len(prop_genes)} genes")
print(f"Final common drugs: {len(common_drugs)}, common genes: {len(common_gene_names)}")


# %%
"""
Cell 7: Align All Three Model Matrices to Common Genes and Drugs
-----------------------------------------------------------------
Subsets each model dataframe (gene_names x drugs) to the shared set of gene
names and drug IDs determined in Cell 6, ready for rank aggregation.

Rank validity guard:
  - If genes are dropped: rank values are stale (rank 5 out of 15 000 genes is no
    longer meaningful if only 12 000 genes remain). Ranks are recomputed column-wise
    (per drug) within the subset so rank 1 = best among the remaining genes.
  - If drugs are dropped: no recomputation needed — ranks are per-column (per-drug)
    and dropping a drug column does not affect any other column's ranks.
"""

ultimate_models_df = {}

for model_name in model_names:
    temp_df = final_models_df[model_name]

    genes_dropped = len(temp_df) - len(common_gene_names)        # rows = genes
    drugs_dropped = len(temp_df.columns) - len(common_drugs)     # cols = drugs

    # Subset to common gene names (rows) and common drugs (columns)
    temp_df = temp_df.loc[common_gene_names]
    temp_df = temp_df[common_drugs]

    if genes_dropped > 0:
        # Ranks were computed over a larger gene pool — recompute within the subset.
        # Re-applying rankdata on existing rank values is valid because relative ordering
        # is preserved: rank 5 < rank 10 before subsetting implies the same after.
        print(f"{model_name}: {genes_dropped} genes dropped — recomputing ranks over {len(common_gene_names)} remaining genes...")
        temp_df = pd.DataFrame(
            temp_df.apply(lambda x: ss.rankdata(x, method="min", nan_policy='omit')),
            index=common_gene_names,
            columns=common_drugs
        )
        print(f"{model_name}: ranks recomputed ✓")
    else:
        print(f"{model_name}: no genes dropped — rank values are valid ✓")

    if drugs_dropped > 0:
        print(f"{model_name}: {drugs_dropped} drugs dropped — per-drug ranks unaffected ✓")
    else:
        print(f"{model_name}: no drugs dropped ✓")

    print(f"{model_name} aligned shape: {temp_df.shape}")
    ultimate_models_df[model_name] = temp_df

del final_models_df
gc.collect()


# %%
"""
Cell 8: Aggregate Docking, ML, and Propagation Ranks Using Mean Rank Method
----------------------------------------------------------------------------
For each drug, combines the three model ranks across all genes:

  1. Align all model DataFrames to an identical sorted drug column order once
     before the drug loop — avoids pandas alignment overhead per iteration.
  2. For each drug, stack Docking, ML, and Propagation rank columns into a
     numpy array of shape (n_genes x 3) — rows = gene targets, columns = models.
  3. Compute the mean rank per gene across all three models (axis=1).
       mean rank = (docking_rank + ML_rank + propagation_rank) / 3
       Lower mean = gene ranks well in ALL models = better target for that drug.
  4. Re-rank the mean scores with scipy.stats.rankdata (ascending, method='min'):
       rank 1 = lowest mean rank = best combined target for that drug.

NOTE on rk.borda: ranky.borda(method='median') computes traditional Borda scores
  (N + 1 - rank) per model before aggregating, so higher output = better target.
  Feeding that into rankdata ascending inverts the result (rank 1 = worst target).
  Using plain mean(axis=1) stays in the convention (lower = better) throughout,
  so rankdata ascending gives rank 1 = best with no inversion needed.

Correctness guarantee:
  Gene with rank 1 in all three models  → mean 1.0  → final rank 1  (best)
  Gene with rank N in all three models  → mean N    → final rank N  (worst)

Output: final_ranked_protein_drug_df — shape (genes x drugs),
  index   = gene_order (StringDB gene names)
  columns = drug_order (sorted DrugBank IDs)
  values  = final aggregated integer rank (1 = best target for that drug)
"""

# Step 1: Align all models to an identical sorted drug column order.
# common_drugs comes from list(set(...)) so its order is arbitrary.
# Sorting once here means the loop can use .to_numpy() directly.
drug_order = sorted(common_drugs)
for model_name in model_names:
    ultimate_models_df[model_name] = ultimate_models_df[model_name][drug_order]

# Verify column alignment before entering the expensive drug loop
assert all(
    ultimate_models_df[model_names[0]].columns.equals(ultimate_models_df[m].columns)
    for m in model_names[1:]
), "Drug column mismatch between models after alignment — cannot safely stack arrays"
print(f"Drug columns aligned across all models ✓  ({len(drug_order)} drugs)")

gene_order = ultimate_models_df[model_names[0]].index.tolist()  # consistent gene row order

# Step 2-4: Mean rank aggregation over all drugs
all_drug_ranks = []  # collect per-drug rank arrays; single DataFrame build at the end

for k, drug in enumerate(drug_order):
    # Stack Docking, ML, and Propagation rank rows for this drug: shape (n_genes x 3).
    # .to_numpy() is safe because all DataFrames share the same gene_order index.
    gene_matrix = np.column_stack(
        [ultimate_models_df[m][drug].to_numpy() for m in model_names]
    )

    # Mean rank across models: lower mean = gene scores well across all models.
    mean_ranks = gene_matrix.mean(axis=1)

    # Re-rank ascending: rank 1 = smallest mean rank = best target for this drug.
    drug_gene_ranks = ss.rankdata(mean_ranks, method='min', nan_policy='omit')
    all_drug_ranks.append(drug_gene_ranks)

    if k % 500 == 0:
        print(f"Processed {k} / {len(drug_order)} drugs...")

# Build the final matrix in one operation.
# Rows = genes (gene_order), columns = drugs (drug_order).
final_ranked_protein_drug_df = pd.DataFrame(
    np.column_stack(all_drug_ranks),
    index=gene_order,
    columns=drug_order
)

del all_drug_ranks, ultimate_models_df
gc.collect()


# %%
# Sanity check: rank of BCL2 for Venetoclax (DB11581) in the final combined ranking.
# Venetoclax is a known BCL2 inhibitor — BCL2 should appear near the top for DB11581.
final_ranked_protein_drug_df.loc["BCL2", "DB11581"]

# %%
"""
Cell 9: Replace DrugBank IDs with Human-Readable Drug Names
------------------------------------------------------------
The rank matrix has DrugBank IDs as column headers (genes x drugs layout).
This cell maps each column ID to its drug name, making results easier to interpret.

Alignment guarantee: drug_names is built by iterating final_ranked_protein_drug_df.columns
(which is drug_order = sorted(common_drugs) set in cell 8), NOT common_drugs whose
order is arbitrary from list(set(...)). Using common_drugs would misalign names
with columns — e.g. Venetoclax (DB11581) would be pasted onto the wrong column.

Fallback: if a DrugBank ID has no name entry, the ID itself is kept as the label.
Duplicate names are reported so downstream .loc[] lookups are not silently ambiguous.
"""

# Load DrugBank ID -> drug name mapping
drugbank_df = pd.read_csv("../Data/drug_protein_interactions.csv", header='infer', sep='|')

# Build O(1) lookup: drug_id -> name (keep first occurrence per drug_id)
drug_id_to_name = (
    drugbank_df[["drug_id", "name"]]
    .drop_duplicates(subset="drug_id")
    .set_index("drug_id")["name"]
)

# Map each drug ID to its name following the exact column order of the rank matrix.
# final_ranked_protein_drug_df.columns = drug_order (sorted), NOT common_drugs (set order).
drug_names = [
    str(drug_id_to_name.get(drug_id, drug_id))
    for drug_id in final_ranked_protein_drug_df.columns
]

# Warn if any IDs were not found in the mapping
missing = [d for d in final_ranked_protein_drug_df.columns if d not in drug_id_to_name]
if missing:
    print(f"WARNING: {len(missing)} drug IDs have no name entry and will keep their ID as label: {missing[:5]}{'...' if len(missing) > 5 else ''}")

# Warn if duplicate names exist (would make .loc[] return multiple columns downstream)
name_counts = pd.Series(drug_names).value_counts()
duplicates = name_counts[name_counts > 1]
if not duplicates.empty:
    print(f"WARNING: {len(duplicates)} drug names appear more than once (different DB IDs, same name).")
    print(duplicates.head(10))

print(f"Total drug labels: {len(drug_names)}, unique: {len(set(drug_names))}")
print(final_ranked_protein_drug_df.head())

# Optional: save with DrugBank IDs (genes as rows, drug IDs as columns)
final_ranked_protein_drug_df.to_csv(
    "../Results/Final_Ranking/drug_protein_ranking_combined_with_drug_ids.csv", sep=",")

# Replace column headers with drug names (copy first to preserve the DrugBank-ID version until del)
final_ranked_protein_drug_names_df = final_ranked_protein_drug_df.copy()
del final_ranked_protein_drug_df
gc.collect()

final_ranked_protein_drug_names_df.columns = drug_names

# Sanity check: verify Venetoclax (DB11581) name lands on the correct column
venetoclax_name = str(drug_id_to_name.get("DB11581", "DB11581"))
if venetoclax_name in final_ranked_protein_drug_names_df.columns:
    print(f"\nSanity check — BCL2 rank for {venetoclax_name}: {final_ranked_protein_drug_names_df.loc['BCL2', venetoclax_name]}")

# Optional: save with drug names (genes as rows, drug names as columns)
final_ranked_protein_drug_names_df.to_csv(
    "../Results/Final_Ranking/drug_protein_ranking_combined_with_drug_names.csv", sep="|")


# %%
