import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Load your data into a DataFrame (assuming 'data.csv' contains your data)
df = pd.read_csv('../Data/drug_protein_interactions.csv', sep="|")
df.columns

# Define indices or identifiers of molecules to remove due to bad conformer IDs
indices_to_remove = []

# Process each molecule
for idx, row in df.iterrows():
    try:
        drug_id = row['drug_id']  #Assuming 'drug_id' is the column name
        
        # Create RDKit molecule object from SMILES
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol is None:
            print(f"Failed to create molecule from SMILES: {row['smiles']}")
            indices_to_remove.append(idx)
            continue
        
        # Add hydrogens to the molecule
        mol = Chem.AddHs(mol)
        
        # Generate conformers
        AllChem.EmbedMultipleConfs(mol, numConfs=1)  # Generates 1 conformer
        
        # Perform UFF optimization on each conformer
        for conf in mol.GetConformers():
            try:
                AllChem.UFFOptimizeMolecule(mol, confId=conf.GetId(), maxIters=1000)
            except Exception as e:
                print(f"UFF optimization failed for molecule with DrugbankID {drug_id}: {row['smiles']}")
                print(e)
        
        # Save the optimized molecule to an SDF file
        sdf_filename = f'../Data/drug_sdfs/mol_{drug_id}.sdf'  # Naming based on DrugbankID
        writer = Chem.SDWriter(sdf_filename)
        writer.write(mol)
        writer.close()
        
    except Exception as e:
        print(f"Processing failed for molecule with DrugbankID {drug_id}: {row['smiles']}")
        print(e)
        indices_to_remove.append(idx)

print(len(indices_to_remove))


# +
## Optionally, save the processed DataFrame back to a CSV file without the removed molecules
#df_cleaned = df.drop(indices_to_remove)
#df_cleaned.to_csv('../data/cleaned_drug_protein_interactions.csv', index=False, sep="|")

