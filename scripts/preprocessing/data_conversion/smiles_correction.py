import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def correct_valency(smiles):
    try:
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: {smiles}")
            return None
        
        # Attempt to sanitize the molecule
        try:
            Chem.SanitizeMol(mol)
            # Convert sanitized molecule back to SMILES
            corrected_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            return corrected_smiles
        except Chem.AtomValenceException as e:
            print(f"Valency error with SMILES {smiles}: {e}")
            return None
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {e}")
        return None

# Load the CSV file with incorrect SMILES
df = pd.read_csv('All_858.csv')

# Apply the correction function to the SMILES column
df['Corrected_SMILES'] = df['SMILES'].apply(correct_valency)

# Log the number of corrected and invalid SMILES
num_corrected = df['Corrected_SMILES'].notna().sum()
num_invalid = df['Corrected_SMILES'].isna().sum()
print(f"Number of corrected SMILES: {num_corrected}")
print(f"Number of invalid SMILES: {num_invalid}")

# Save the results to a new CSV file
df.to_csv('corrected_smiles.csv', index=False)

# Check the first few rows of the output file
print(df.head())
