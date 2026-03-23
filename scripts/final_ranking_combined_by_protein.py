# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: BeatAML2.0
#     language: python
#     name: python3
# ---

# %%
"""
Cell 1: Import Required Libraries
----------------------------------
pandas, numpy  : data loading and matrix manipulation
matplotlib     : plotting (available but not used in core pipeline)
gc             : manual garbage collection to manage memory on large dataframes
re             : regular expressions (available for string parsing if needed)
scipy.stats    : rankdata() used to convert raw docking scores to ranks
ranky          : Borda rank aggregation across multiple scoring models
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
Cell 2: Load and Preprocess Docking Scores and ML Virtual Screening Data
-------------------------------------------------------------------------
Docking Scores (tab-separated):
  - Rows = drug IDs, Columns = protein targets (UniProt/AlphaFold IDs)
  - Raw docking scores are in kcal/mol; more negative = stronger predicted binding.
  - Scores are converted to integer ranks per protein using scipy.stats.rankdata
    with ascending order: rank 1 = minimum (most negative) docking score = best binder.
  - Ties are broken by assigning the minimum rank (method="min").
  - NaN values are ignored during ranking (nan_policy="omit").

ML Virtual Screening (comma-separated):
  - Rows = drug IDs, Columns = protein targets (UniProt/AlphaFold IDs)
  - Values are a pre-computed ranked list where rank 1 = highest ML score = best binder.
  - The ranking is inverted (N + 1 - rank) so that rank 1 = best binder,
    consistent with the docking convention used downstream.

At the end, common drug IDs and common target columns are identified
across both datasets to ensure alignment in downstream steps.
"""

model_names = ['docking', 'ML_virtual_screening']
model_files = [
    '../Results/Docking_Scores/docking_scores.csv',
    '../Results/ML_virtual_screenings/ML_virtual_screening_by_protein.csv'
]

all_models_df = {}    # stores the processed dataframe for each model
common_drug_ids = None
common_targets = None
gene_names = None

for i in range(len(model_names)):
    model_name = model_names[i]

    if "dock" in model_name:
        # Load tab-separated docking score matrix (drugs x targets)
        docking_df = pd.read_csv(model_files[i], header='infer', sep="\t")

        # Set the first column (drug IDs) as the index
        rev_docking_df = pd.DataFrame(docking_df.iloc[:, 1:])
        rev_docking_df.index = docking_df.iloc[:, 0]  # type: ignore
        print(rev_docking_df.head())

        # Convert raw docking scores to ranks per protein (column-wise).
        # rankdata uses ascending order: rank 1 = minimum (most negative) score = best binder.
        # method="min" assigns the lowest rank to ties; nan_policy="omit" skips missing values.
        docking_ranked_df = rev_docking_df.apply(
            lambda x: ss.rankdata(x, method="min", nan_policy='omit')
        )
        del rev_docking_df, docking_df
        gc.collect()
        print(f"Docking ranked dataset shape: {docking_ranked_df.shape}")

        all_models_df[model_name] = docking_ranked_df
        common_target = docking_ranked_df.columns.tolist()
        common_drug_ids = docking_ranked_df.index.tolist()
        del docking_ranked_df
        gc.collect()

    if "ML" in model_name:
        # Load comma-separated ML virtual screening matrix (drugs x targets).
        ml_vs_df = pd.read_csv(model_files[i], header='infer', sep=",")

        # Set the first column (drug IDs) as the index
        ml_ranked_df = pd.DataFrame(ml_vs_df.iloc[:, 1:])
        ml_ranked_df.index = ml_vs_df.iloc[:, 0]  # type: ignore

        # Invert ranks so rank 1 = best predicted binder, consistent with docking.
        # Original ML ranks: rank 1 = lowest ML score = worst binder.
        # After inversion: rank 1 = best binder (rank N = worst).
        n_drugs = ml_ranked_df.shape[0]
        ml_ranked_df = (n_drugs + 1) - ml_ranked_df
        print(ml_ranked_df.head())

        all_models_df[model_name] = ml_ranked_df
        print(f"ML ranked dataset shape: {ml_ranked_df.shape}")

        # Retain only drugs and targets present in both docking and ML datasets
        common_drug_ids = list(set(ml_ranked_df.index.tolist()) & set(common_drug_ids))
        common_target = list(set(common_target) & set(ml_ranked_df.columns.tolist()))
        del ml_ranked_df, ml_vs_df
        gc.collect()

print(f"Length of common drug ids: {len(common_drug_ids)}")  # type: ignore
print(f"Length of common targets: {len(common_target)}")


# %%
"""
Cell 3: Load Gene Name to AlphaFold ID Mapping
------------------------------------------------
Docking and ML screening results use UniProt/AlphaFold IDs as column headers.
This cell loads a mapping file that links AlphaFold IDs to human-readable
StringDB gene names (e.g., P10415 -> BCL2).

The intersection of common targets (from Cell 2) and known AlphaFold IDs
defines the final set of proteins carried forward in the pipeline.
"""

# Load the mapping between StringDB gene names and AlphaFold UniProt IDs
gene_alphafold_mapping_df = pd.read_csv(
    "../Data/StringDB_genes_AlphaFold_Mapping.csv", header='infer', sep="\t"
)
gene_names = list(set(gene_alphafold_mapping_df["StringDB_Gene_Name"].tolist()))
alphafold_ids = list(set(gene_alphafold_mapping_df["AlphaFold_Uniprot_Id"].tolist()))
print(f"Alphafold ids length: {len(alphafold_ids)}, total genes: {len(gene_names)}")

# Intersect with targets present in both docking and ML datasets
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
AlphaFold IDs are replaced with human-readable gene names for both models.
Two sources of redundancy are handled:
  1. One AlphaFold ID -> multiple gene names: the column is used for each gene.
  2. One gene name <- multiple AlphaFold IDs: resolved by canonical preference:
       - Canonical UniProt/Swiss-Prot IDs are exactly 6 characters (e.g. P10415).
       - TrEMBL/isoform entries are longer (e.g. A0A024RBG1 = 10 chars).
       - If any canonical ID exists for a gene, use min rank across canonical
         structures only (ignore non-canonical structures for that gene).
       - If no canonical ID exists, fall back to min rank across all structures.

This ensures BCL2 uses P10415 rather than an arbitrary TrEMBL entry, giving
biologically meaningful results for well-characterised human proteins.

Output: final_models_df — dict keyed by model name, each value is a
DataFrame of shape (drugs x gene_names) with resolved rank values.
"""

# Pre-build lookup dict once: AlphaFold ID -> [gene_name, ...]
af_to_genes = (
    subset_gene_alphafold_mapping_df
    .groupby("AlphaFold_Uniprot_Id")["StringDB_Gene_Name"]
    .apply(list)
    .to_dict()
)

def is_canonical(af_id):
    """
    Return True if af_id is a canonical Swiss-Prot UniProt accession.
    Canonical IDs are exactly 6 characters (e.g. P10415, Q9BXK5).
    TrEMBL/isoform entries are longer (e.g. A0A024RBG1).
    """
    return len(af_id) == 6

def expand_and_resolve(model_df, common_drug_ids, common_targets, af_to_genes):
    """
    Subset to common drugs/targets, expand AlphaFold IDs to gene names,
    and resolve multi-mapping with canonical preference:
      - canonical structures (Swiss-Prot, 6-char ID): min rank across these
      - non-canonical structures: used only when no canonical ID exists for the gene

    Parameters
    ----------
    model_df        : DataFrame (drugs x AlphaFold IDs) with rank values.
    common_drug_ids : list of drug IDs to retain (row subset).
    common_targets  : list of AlphaFold IDs to expand (column subset).
    af_to_genes     : dict {alphafold_id: [gene_name, ...]} pre-built from mapping.

    Returns
    -------
    DataFrame of shape (drugs x unique_gene_names).
    """
    df = model_df.loc[common_drug_ids, common_targets]
    print(f"  Subset shape (drugs x AlphaFold IDs): {df.shape}")

    canonical_arrays = {}     # gene -> min rank across canonical structures only
    noncanonical_arrays = {}  # gene -> min rank across all structures (fallback)

    for target in common_targets:
        col = df[target].to_numpy()
        canonical = is_canonical(target)

        for gene in af_to_genes.get(target, []):
            # Always update the all-structures fallback
            if gene in noncanonical_arrays:
                np.minimum(noncanonical_arrays[gene], col, out=noncanonical_arrays[gene])
            else:
                noncanonical_arrays[gene] = col.copy()

            # Update canonical-only pool if this structure is canonical
            if canonical:
                if gene in canonical_arrays:
                    np.minimum(canonical_arrays[gene], col, out=canonical_arrays[gene])
                else:
                    canonical_arrays[gene] = col.copy()

    # Prefer canonical ranks; fall back to all-structures min if no canonical exists
    gene_arrays = {
        gene: canonical_arrays.get(gene, noncanonical_arrays[gene])
        for gene in noncanonical_arrays
    }

    n_canonical = sum(1 for g in noncanonical_arrays if g in canonical_arrays)
    n_fallback  = len(noncanonical_arrays) - n_canonical
    print(f"  Genes resolved via canonical ID: {n_canonical}, via fallback (no canonical): {n_fallback}")

    result = pd.DataFrame(gene_arrays, index=common_drug_ids)
    print(f"  Expanded shape (drugs x gene names): {result.shape}")
    return result


final_models_df = {}

for model_name in all_models_df.keys():
    print(f"Processing model: {model_name}")
    result = expand_and_resolve(
        all_models_df[model_name], common_drug_ids, common_targets, af_to_genes
    )
    print(result.head())
    final_models_df[model_name] = result
    del result
    gc.collect()


# %%
drug_id = "DB11581"
gene   = "BCL2"

for model_name in model_names:
    series    = final_models_df[model_name][gene]
    rank_val  = series.loc[drug_id]                          # rank assigned by rankdata
    sorted_ids = series.sort_values(ascending=True).index.tolist()
    position  = sorted_ids.index(drug_id)                   # 0-indexed position in sorted list

    # How many drugs share the same rank (ties)?
    n_tied = (series == rank_val).sum()

    print(f"--- {model_name} ---")
    print(f"  rank value (from rankdata) : {rank_val}")
    print(f"  position in sorted list    : {position} (0-indexed)  →  {position+1} (1-indexed)")
    print(f"  drugs sharing rank {rank_val}      : {n_tied}")
    if n_tied > 1:
        print(f"  → tie causes position ({position+1}) to differ from rank ({rank_val}) by up to {n_tied-1}")
    else:
        print(f"  → no tie: position+1 == rank  ✓")
    print()


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

def get_drug_ids_for_gene(gene_name):
    """
    Return drugs sorted by rank (ascending) for a given gene, for each model.
    Rank 1 = best predicted binder.

    Returns
    -------
    dict of {model_name: pd.Series}
        Each Series has drug_id as index and rank value as values, sorted
        ascending by rank. Use series.index to get ordered drug IDs;
        use series.loc[drug_id] to look up a specific drug's rank.

    NOTE: Uses final_models_df — only valid between cells 4 and 7.
    After cell 7 runs (which deletes final_models_df), use
    final_ranked_protein_drug_df or final_ranked_protein_drug_names_df instead.
    """
    result = {}
    for model_name in model_names:
        result[model_name] = (
            final_models_df[model_name][gene_name]
            .sort_values(ascending=True)
            .rename("rank")
        )
    return result

# Sanity check: BCL2 maps to multiple AlphaFold IDs (many-to-many mapping).
# P10415 is the canonical human BCL2 UniProt ID.
bcl2_alphafold_ids = get_alphafold_ids("BCL2")
print(f"BCL2 AlphaFold IDs ({len(bcl2_alphafold_ids)}): {bcl2_alphafold_ids}")

# Each AlphaFold ID may itself map to multiple gene names — show all genes per ID
for af_id in bcl2_alphafold_ids:
    gene_names_for_id = get_gene_names(af_id)
    print(f"  {af_id} -> {gene_names_for_id}")

# Confirm P10415 (canonical BCL2) is present
canonical_id = "P10415"
print(f"\nCanonical BCL2 ID {canonical_id} in mapping: {canonical_id in bcl2_alphafold_ids}")

# Show top 5 drugs for BCL2 from each model, with their rank values
bcl2_drug_ids = get_drug_ids_for_gene("BCL2")
for model_name, series in bcl2_drug_ids.items():
    print(f"\n{model_name} — top 5 drugs for BCL2:")
    print(series.head())
    print(f"  DB11581 rank: {series.loc['DB11581']}, position: {series.index.get_loc('DB11581') + 1}")


# %%
"""
Cell 6: Identify Final Common Drugs and Gene Names
---------------------------------------------------
After gene-name mapping (Cell 4), the docking and ML dataframes may differ
slightly in covered drugs and gene names. This cell computes the final
intersection to ensure both matrices are perfectly aligned before aggregation.
"""

# Collect drugs and genes from both models
common_drugs = final_models_df['docking'].index.tolist()
common_gene_names = final_models_df['docking'].columns.tolist()
ml_drugs = final_models_df['ML_virtual_screening'].index.tolist()
ml_gene_names = final_models_df['ML_virtual_screening'].columns.tolist()

# Keep only drugs and genes present in both models
common_gene_names = list(set(common_gene_names) & set(ml_gene_names))
common_drugs = list(set(common_drugs) & set(ml_drugs))
print(f"Final common drugs: {len(common_drugs)}, common genes: {len(common_gene_names)}")

# %%
"""
Cell 7: Align Both Model Matrices to Common Drugs and Genes
------------------------------------------------------------
Subsets each model dataframe to the shared set of drug IDs and gene names
determined in Cell 6, producing aligned matrices ready for rank aggregation.

Rank validity guard:
  - If drugs are dropped: rank values are stale (rank 194 out of 4275 is no longer
    meaningful if only 4000 drugs remain). Ranks are recomputed column-wise within
    the subset so rank 1 = best among the remaining drugs.
  - If genes are dropped: no recomputation needed — ranks are per-column (per-gene)
    and dropping a gene column does not affect any other column's ranks.
"""

ultimate_models_df = {}

for model_name in model_names:
    temp_df = final_models_df[model_name]

    drugs_dropped = len(temp_df) - len(common_drugs)
    genes_dropped = len(temp_df.columns) - len(common_gene_names)

    # Subset to common drugs (rows) and common gene names (columns)
    temp_df = temp_df.loc[common_drugs]
    temp_df = temp_df[common_gene_names]

    if drugs_dropped > 0:
        # Ranks were computed over a larger drug pool — recompute within the subset.
        # Re-applying rankdata on existing rank values is valid because relative ordering
        # is preserved: rank 5 < rank 10 before subsetting implies the same after.
        print(f"{model_name}: {drugs_dropped} drugs dropped — recomputing ranks over {len(common_drugs)} remaining drugs...")
        temp_df = pd.DataFrame(
            temp_df.apply(lambda x: ss.rankdata(x, method="min", nan_policy='omit')),
            index=common_drugs,
            columns=common_gene_names
        )
        print(f"{model_name}: ranks recomputed ✓")
    else:
        print(f"{model_name}: no drugs dropped — rank values are valid ✓")

    if genes_dropped > 0:
        print(f"{model_name}: {genes_dropped} genes dropped — per-gene ranks unaffected ✓")
    else:
        print(f"{model_name}: no genes dropped ✓")

    print(f"{model_name} aligned shape: {temp_df.shape}")
    ultimate_models_df[model_name] = temp_df

del final_models_df
gc.collect()


# %%
#Find rank of DB11581 for BCL2 in both models
for model_name in model_names:
    rank = ultimate_models_df[model_name]["BCL2"].loc["DB11581"]
    total = len(ultimate_models_df[model_name])
    print(f"{model_name}: {drug_id} rank for BCL2 = {rank} / {total} (top {100*rank/total:.1f}%)")

# %%
"""
Cell 8: Aggregate Docking and ML Ranks Using Mean Rank Method
--------------------------------------------------------------
For each protein (gene), combines the docking and ML ranks across all drugs:

  1. Align both model DataFrames to an identical sorted drug order once
     before the gene loop — avoids pandas index-alignment overhead on
     every one of the ~15K gene iterations.
  2. For each gene, stack docking and ML rank columns into a numpy array
     of shape (n_drugs x 2) — rows = drug candidates, columns = models.
  3. Compute the mean rank per drug across both models (axis=1).
       mean rank = (docking_rank + ML_rank) / 2
       Lower mean = drug ranks well in BOTH models = better candidate.
  4. Re-rank the mean scores with scipy.stats.rankdata (ascending, method='min'):
       rank 1 = lowest mean rank = best combined candidate for that protein.

NOTE on rk.borda: ranky.borda(method='mean') computes traditional Borda scores
  (N + 1 - rank) per model before averaging, so higher output = better binder.
  Feeding that into rankdata ascending inverts the result (rank 1 = worst binder).
  Example — DB11581/BCL2 (docking 194, ML 84, N=4275):
    rk.borda mean = ((4275+1-194) + (4275+1-84)) / 2 = 4137  (high = good)
    rankdata ascending of 4137 among 4275 drugs → rank 4052  (looks bad, but is good)
  Using plain mean(axis=1) stays in our convention (lower = better) throughout,
  so rankdata ascending gives rank 1 = best with no inversion.

Correctness guarantee:
  Drug with docking rank 1 and ML rank 1   → mean 1.0   → final rank 1   (best)
  Drug with docking rank N and ML rank N   → mean N     → final rank N   (worst)
  Drug with docking rank 194 and ML rank 84 → mean 139  → final rank ~139 ✓

Output: final_ranked_protein_drug_df — shape (drugs x genes),
  index   = drug_order (sorted DrugBank IDs)
  columns = common_gene_names (StringDB gene names)
  values  = final aggregated integer rank (1 = best)
"""

# Step 1: Align both models to an identical sorted drug order.
# common_drugs comes from list(set(...)) so its order is arbitrary and may differ
# between the two models. Sorting once here means the loop can use .to_numpy()
# directly (no per-iteration index alignment).
drug_order = sorted(common_drugs)
for model_name in model_names:
    ultimate_models_df[model_name] = ultimate_models_df[model_name].loc[drug_order]

# Verify alignment before entering the expensive gene loop
assert all(
    ultimate_models_df[model_names[0]].index.equals(ultimate_models_df[m].index)
    for m in model_names[1:]
), "Drug index mismatch between models after alignment — cannot safely stack arrays"
print(f"Drug indices aligned across all models ✓  ({len(drug_order)} drugs)")

# Step 2-4: Mean rank aggregation over all genes
all_gene_ranks = []  # collect per-gene rank arrays; single DataFrame build at the end

for k, gene_name in enumerate(common_gene_names):
    # Stack docking and ML rank columns: shape (n_drugs x 2).
    # .to_numpy() is safe because both DataFrames share the same drug_order index.
    gene_matrix = np.column_stack(
        [ultimate_models_df[m][gene_name].to_numpy() for m in model_names]
    )

    # Mean rank across models: lower mean = drug scores well in both models.
    mean_ranks = gene_matrix.mean(axis=1)

    # Re-rank ascending: rank 1 = smallest mean rank = best drug for this protein.
    gene_ranks = ss.rankdata(mean_ranks, method='min', nan_policy='omit')
    all_gene_ranks.append(gene_ranks)

    if k % 500 == 0:
        print(f"Processed {k} / {len(common_gene_names)} proteins...")

# Build the final matrix in one operation.
# Rows = drugs (drug_order), columns = genes (common_gene_names iteration order).
final_ranked_protein_drug_df = pd.DataFrame(
    np.column_stack(all_gene_ranks),
    index=drug_order,
    columns=common_gene_names
)

del all_gene_ranks, ultimate_models_df
gc.collect()


# %%
final_ranked_protein_drug_df["BCL2"].loc["DB11581"]  # rank of DB11581 for BCL2 in the final combined ranking

# %%
"""
Cell 9: Replace DrugBank IDs with Human-Readable Drug Names
------------------------------------------------------------
The rank matrix currently uses DrugBank IDs (e.g., DB00404) as the row index.
This cell maps each ID to its drug name using the drug-protein interaction
reference file, making results easier to interpret.

Alignment guarantee: drug_names is built by iterating final_ranked_protein_drug_df.index
(which is drug_order = sorted(common_drugs) set in cell 8), NOT common_drugs whose
order is arbitrary from list(set(...)). Using common_drugs would misalign names
with rows — e.g. Venetoclax (DB11581) would be pasted onto the wrong row.

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

# Map each drug ID to its name following the exact row order of the rank matrix.
# final_ranked_protein_drug_df.index = drug_order (sorted), NOT common_drugs (set order).
drug_names = [
    str(drug_id_to_name.get(drug_id, drug_id))
    for drug_id in final_ranked_protein_drug_df.index
]

# Warn if any IDs were not found in the mapping
missing = [d for d in final_ranked_protein_drug_df.index if d not in drug_id_to_name]
if missing:
    print(f"WARNING: {len(missing)} drug IDs have no name entry and will keep their ID as label: {missing[:5]}{'...' if len(missing) > 5 else ''}")

# Warn if duplicate names exist (would make .loc[] return multiple rows downstream)
name_counts = pd.Series(drug_names).value_counts()
duplicates = name_counts[name_counts > 1]
if not duplicates.empty:
    print(f"WARNING: {len(duplicates)} drug names appear more than once (different DB IDs, same name).")
    print(duplicates.head(10))

print(f"Total drug labels: {len(drug_names)}, unique: {len(set(drug_names))}")
print(final_ranked_protein_drug_df.head())

# Optional: save with DrugBank IDs
final_ranked_protein_drug_df.to_csv(
    "../Results/Final_Ranking/protein_drug_ranking_with_drug_ids.csv", sep=",")

# Replace index with drug names (copy first so the DrugBank-ID version is preserved until del)
final_ranked_protein_drug_names_df = final_ranked_protein_drug_df.copy()
del final_ranked_protein_drug_df
gc.collect()

final_ranked_protein_drug_names_df.index = drug_names

# Sanity check: verify Venetoclax (DB11581) name lands on the correct row
venetoclax_name = str(drug_id_to_name.get("DB11581", "DB11581"))
if venetoclax_name in final_ranked_protein_drug_names_df.index:
    print(f"\nSanity check — {venetoclax_name} BCL2 rank: {final_ranked_protein_drug_names_df.loc[venetoclax_name, 'BCL2']}")

# Optional: save with drug names
final_ranked_protein_drug_names_df.to_csv(
    "../Results/Final_Ranking/protein_drug_ranking_with_drug_names.csv", sep="|")


# %%
"""
Cell 10: Query Top Drug Candidates for a Given Protein
-------------------------------------------------------
Defines get_top_drug_candidates_for_protein() which retrieves the top N
drugs for a specified gene from the aggregated rank matrix.
Rank 1 = best combined candidate (lowest aggregated Borda rank).

Example: top 10 drugs for NLRP3 are retrieved and saved to CSV.
"""

def get_top_drug_candidates_for_protein(gene_name, top_n=10):
    """
    Return the top N drug candidates for a given protein target.

    Parameters
    ----------
    gene_name : str
        StringDB gene name (e.g., "BCL2", "NLRP3").
    top_n : int
        Number of top candidates to return (default 10).

    Returns
    -------
    pd.Series or None
        Drug names with their aggregated rank (ascending); None if gene not found.
    """
    if gene_name not in final_ranked_protein_drug_names_df.columns:
        print(f"Gene name '{gene_name}' not found in the ranked matrix.")
        return None
    # nsmallest: rank 1 = best combined candidate
    return final_ranked_protein_drug_names_df[gene_name].nsmallest(top_n)

# Example: retrieve and save top 10 drugs for NLRP3
top_drugs = get_top_drug_candidates_for_protein("ZBP1", top_n=20)
print(top_drugs)


# %%
"""
Cell 11: Retrieve Top 1000 Drug Candidates for NLRP3
-----------------------------------------------------
Extends the query from Cell 10 to return the top 1000 ranked drugs
for NLRP3. Useful for downstream enrichment or further filtering.
"""

final_ranked_protein_drug_names_df["NLRP3"].nsmallest(1000)

# %%
