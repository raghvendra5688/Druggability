import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Load your data into a DataFrame (assuming 'SMILES.csv' contains your data)
df = pd.read_csv('SMILES.csv')

# Define indices or identifiers of molecules to remove due to bad conformer IDs
indices_to_remove = []

# Process each molecule
for idx, smiles in enumerate(df['SMILES']):
    try:
        # Create RDKit molecule object from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Failed to create molecule from SMILES: {smiles}")
            continue
        
        # Add hydrogens to the molecule
        mol = Chem.AddHs(mol)
        
        # Generate conformers
        AllChem.EmbedMultipleConfs(mol, numConfs=1)  # Generates 1 conformer
        
        # Perform UFF optimization on each conformer
        for conf in mol.GetConformers():
            AllChem.UFFOptimizeMolecule(mol, confId=conf.GetId(), maxIters=1000)
        
        # Save the optimized molecule to an SDF file
        writer = Chem.SDWriter(f'molecule_{idx}.sdf')
        writer.write(mol)
        writer.close()
        
    except Exception as e:
        print(f"Optimization failed for molecule with SMILES: {smiles}")
        print(e)
        indices_to_remove.append(idx)

# Optionally, save the processed DataFrame back to a CSV file without the removed molecules
df_cleaned = df.drop(indices_to_remove)
df_cleaned.to_csv('cleaned_data.csv', index=False)

