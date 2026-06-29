#!/export/home/rmall/micromamba/envs/effichem/bin/python
# Run via: scripts/run_panviral_predictions.sh

import os
import pandas as pd
import numpy as np
import torch
import gc
import time
import scipy.stats as ss
from DeepPurpose import DTI as models

def borda_median(df):
    """Borda aggregation via column-wise median rank (equivalent to rk.borda(df, method='median'))."""
    return df.median(axis=1)

BASE_DIR    = "/export/qcai-omics/Raghvendra/Druggability"
MODELS_DIR  = os.path.join(BASE_DIR, "Models/pretrained_models")
RESULTS_DIR = os.path.join(BASE_DIR, "Results/Chirag_Analysis/PanViral_screenings")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Load compounds and proteins
compound_df   = pd.read_csv(os.path.join(BASE_DIR, "Data/PanViral/compound_info.csv"))
protein_df    = pd.read_csv(os.path.join(BASE_DIR, "Data/PanViral/protein_info.csv"))

drug_smiles   = compound_df["Canonical_SMILES"].tolist()
drug_names    = compound_df["InchiKey"].tolist()
protein_seqs  = protein_df["Sequence"].tolist()
protein_names = protein_df["UniProtId"].tolist()

N_drugs    = len(drug_smiles)
N_proteins = len(protein_seqs)
print(f"{N_drugs} drugs × {N_proteins} proteins = {N_drugs*N_proteins} pairs")

# Build flat all-pairs lists
revised_drug_smiles, revised_drug_names, revised_protein_seqs, revised_protein_names = [], [], [], []
for i in range(N_drugs):
    for j in range(N_proteins):
        revised_drug_smiles.append(drug_smiles[i])
        revised_drug_names.append(drug_names[i])
        revised_protein_seqs.append(protein_seqs[j])
        revised_protein_names.append(protein_names[j])

model_names = [
    "CNN_CNN_BindingDB_IC50", "MPNN_CNN_DAVIS", "Morgan_CNN_BindingDB_IC50",
    "Morgan_AAC_BindingDB_IC50", "Daylight_AAC_BindingDB_IC50", "CNN_CNN_DAVIS",
    "CNN_CNN_BindingDB", "Morgan_CNN_BindingDB", "Morgan_CNN_DAVIS",
    "MPNN_CNN_BindingDB", "MPNN_CNN_KIBA", "Transformer_CNN_BindingDB",
    "Daylight_AAC_DAVIS", "Daylight_AAC_KIBA", "Daylight_AAC_BindingDB",
    "Morgan_AAC_BindingDB", "Morgan_AAC_KIBA", "Morgan_AAC_DAVIS",
]

# ── Phase 1: Predictions (skip models already done) ──────────────────────────
for model_name in model_names:
    final_csv = os.path.join(RESULTS_DIR, model_name, f"{model_name}_virtual_screening.csv")
    if os.path.exists(final_csv):
        print(f"[skip] {model_name}")
        continue

    print(f"\n>>> {model_name}")
    t0 = time.time()

    result_folder = os.path.join(RESULTS_DIR, model_name)
    os.makedirs(result_folder, exist_ok=True)

    model = models.model_pretrained(os.path.join(MODELS_DIR, model_name))
    model.config['batch_size'] = 1024

    models.virtual_screening(
        revised_drug_smiles, revised_protein_seqs, model,
        drug_names=revised_drug_names, target_names=revised_protein_names,
        result_folder=result_folder, convert_y=False,
        output_num_max=len(revised_drug_smiles), verbose=True,
    )

    # Explicitly unload model from GPU before moving on
    del model
    torch.cuda.empty_cache()
    gc.collect()

    # Rename DeepPurpose output to our convention
    raw = os.path.join(result_folder, "virtual_screening.txt")
    if os.path.exists(raw):
        os.rename(raw, final_csv)
    print(f"    done in {time.time()-t0:.1f}s")

# ── Phase 2: Rank aggregation ─────────────────────────────────────────────────
all_model_df = {}
for model_name in model_names:
    path = os.path.join(RESULTS_DIR, model_name, f"{model_name}_virtual_screening.csv")
    if not os.path.exists(path):
        print(f"[missing] {model_name}")
        continue

    df = pd.read_csv(path, header="infer", sep="|")
    df = df.drop(columns=[c for c in df.columns if "Unnamed" in str(c)], errors="ignore")
    # Columns after dropping index/unnamed: Rank, Drug Name, Target Name, Binding Score
    df.columns = ["Rank", "Drug Name", "Target Name", "Binding Score"]
    df = df.drop(columns=["Rank"])
    df["Drug Name"]   = df["Drug Name"].str.strip()
    df["Target Name"] = df["Target Name"].str.strip()
    df = df.dropna(how="all")
    df["Binding Score"] = df["Binding Score"].astype(float)

    # Pivot to protein × drug matrix
    # Rank drugs per protein (axis=1): each row = one protein, rank all 2216 drugs (1 = best binder)
    pivoted = df.pivot(index="Target Name", columns="Drug Name", values="Binding Score")
    ranked  = pivoted.apply(lambda x: ss.rankdata(-x, method="min", nan_policy="omit"), axis=1,
                            result_type="broadcast")
    ranked.index   = pivoted.index
    ranked.columns = pivoted.columns
    all_model_df[model_name] = ranked
    print(f"loaded {model_name}: {ranked.shape}")

# Borda rank aggregation across models
# For each drug, collect its per-protein drug-ranks from every model,
# then take the median (Borda) → final drug rank per protein (1–N_drugs)
available_models = list(all_model_df.keys())
ref_df       = all_model_df[available_models[0]]
all_drug_ids = ref_df.columns.tolist()
all_prot_ids = ref_df.index.tolist()

agg_rows = {}
for k, drug in enumerate(all_drug_ids):
    # temp_df: proteins × models, values = rank of this drug for each protein in each model (1–N_drugs)
    temp_df = pd.DataFrame({mn: all_model_df[mn][drug]
                            for mn in available_models if drug in all_model_df[mn].columns},
                           index=all_prot_ids)
    # Borda median: median drug rank across models per protein
    agg_rows[drug] = borda_median(temp_df)
    if (k + 1) % 500 == 0:
        print(f"  aggregated {k+1}/{len(all_drug_ids)} drugs")

# Build drug × protein borda-score matrix
borda_matrix = pd.DataFrame(agg_rows).T  # drugs × proteins, values = median rank (not yet 1–N_drugs)

# Re-rank drugs within each protein column → true 1–N_drugs per protein
drug_protein_matrix = borda_matrix.apply(
    lambda col: ss.rankdata(col, method="min", nan_policy="omit"), axis=0
).astype(int)
drug_protein_matrix.index   = borda_matrix.index
drug_protein_matrix.columns = borda_matrix.columns
drug_protein_matrix.index.name   = "Drug_InchiKey"
drug_protein_matrix.columns.name = "UniProtId"

# Aggregate rank: mean of per-protein ranks across all proteins, then rank drugs 1–N_drugs
drug_protein_matrix["Mean_Rank"]      = drug_protein_matrix[all_prot_ids].mean(axis=1)
drug_protein_matrix["Aggregate_Rank"] = ss.rankdata(
    drug_protein_matrix["Mean_Rank"], method="min"
).astype(int)
drug_protein_matrix = drug_protein_matrix.sort_values("Aggregate_Rank")

out_path = os.path.join(BASE_DIR, "Results/Chirag_Analysis", "PanViral_drug_protein_ranks.csv")
drug_protein_matrix.to_csv(out_path, index=True)
print(f"\nSaved {drug_protein_matrix.shape} matrix → {out_path}")
print(drug_protein_matrix[["Aggregate_Rank", "Mean_Rank"] + all_prot_ids].head(10).to_string())
