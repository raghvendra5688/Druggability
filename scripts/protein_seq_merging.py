# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: ProtSol
#     language: python
#     name: protsol
# ---

# ## Merge the protein sequence with UniProt Id and gene names for AlphaFold Database proteins

import pandas as pd
import numpy as np
from Bio import SeqIO

#Load the training set
input_file = "../Data/UP000005640_9606.fasta"
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
seq_names, seq_gene_names, seq_seqs = [],[],[]
for fasta in fasta_sequences:
    fasta_id = fasta.id.split("|")
    name, gene_name, sequence = fasta_id[1], fasta_id[2], str(fasta.seq)
    seq_names.append(name)
    seq_gene_names.append(gene_name)
    seq_seqs.append(sequence)
protein_df = pd.DataFrame([seq_names, seq_gene_names, seq_seqs]).transpose()
protein_df.columns = ["UniProtId","Gene_Name","Protein_Sequence"]
protein_df.head()

# +
#Get the filtered AlphaFold Database
alphafold_protein_df = pd.read_csv("../Data/Filtered_AlphaFold_UniProtKB.csv",header="infer")
print(alphafold_protein_df.shape)

#Perform innner join
merged_protein_df = pd.merge(alphafold_protein_df, protein_df, on=['UniProtId','Gene_Name'], how='left')
merged_protein_df.to_csv("../Data/Filtered_AlphaFold_UniProtKB_with_Sequence.csv",index=None)
# -

missing_protein_ids = np.where(merged_protein_df["Protein_Sequence"].isna())
missing_protein_ids = list(missing_protein_ids[0])
merged_protein_df.iloc[missing_protein_ids,:]


