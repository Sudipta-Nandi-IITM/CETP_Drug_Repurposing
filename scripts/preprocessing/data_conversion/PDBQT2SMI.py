import os
import pandas as pd
from rdkit import Chem

# Define the input and output folders
input_folder = "/path/to/the/file/with/pdbqts"
output_csv = "Output.csv"

# List all pdbqt files in the input folder
pdbqt_files = [f for f in os.listdir(input_folder) if f.endswith('.pdbqt')]

# Initialize an empty list to store the results
results = []

# Function to convert PDBQT to SMILES
def pdbqt_to_smiles(pdbqt_file):
    try:
        # Read the PDBQT file content
        with open(pdbqt_file, 'r') as f:
            pdbqt_content = f.read()

        # Remove lines that cause issues (e.g., those containing 'A' as element)
        pdbqt_lines = pdbqt_content.split('\n')
        valid_lines = [line for line in pdbqt_lines if not line.startswith('ATOM') or ' A ' not in line]
        pdbqt_content = '\n'.join(valid_lines)

        # Convert PDBQT to Mol
        mol = Chem.MolFromPDBBlock(pdbqt_content, removeHs=False)

        if mol is None:
            return None

        # Convert Mol to SMILES
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception as e:
        print(f"Error processing {pdbqt_file}: {e}")
        return None

# Process each pdbqt file
for pdbqt_file in pdbqt_files:
    input_filepath = os.path.join(input_folder, pdbqt_file)
    smiles = pdbqt_to_smiles(input_filepath)
    if smiles:
        results.append([pdbqt_file, smiles])
    else:
        results.append([pdbqt_file, "Conversion failed"])

# Create a DataFrame and save to CSV
df = pd.DataFrame(results, columns=["Filename", "SMILES"])
df.to_csv(output_csv, index=False)

print(f"Conversion completed. Output saved to {output_csv}")

