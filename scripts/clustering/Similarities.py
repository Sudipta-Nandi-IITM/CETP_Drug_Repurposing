#Calculate similaritites from fingerprints for clustering
#This script will show step wise how to calculate similarities from fingerprints for clustering to a number of defined clusters. The clustering metrics used in this study are given as follows. If you dont want to use any metric, skip the corresponding section.

1.
!pip install rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools, DataStructs
import pandas as pd

# Step 1: Read the CSV file into a Pandas DataFrame
df = pd.read_csv('SMILES.csv')

# Step 2: Convert SMILES strings to RDKit molecules
df['mols'] = df['SMILES'].apply(Chem.MolFromSmiles)

# Calculate the fingerprint from the molecule
df['fingerprints'] = df['mols'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2) if x is not None else None)
df['fingerprints']
df

2.
df.to_csv('molecules_with_fingerprints.csv', index=False)

3.
target = '1.pdbqt'
df.set_index(df['Ligand_name'].apply(lambda x: target if x == target else x), inplace=True)

4. For checking if similarity is calculating from fingerprints:
df.fingerprints['1.pdbqt']
DataStructs.TanimotoSimilarity(df.fingerprints['1.pdbqt'], df.fingerprints['1.pdbqt'])
df.fingerprints['1']
DataStructs.TanimotoSimilarity(df.fingerprints['1.pdbqt'], df.fingerprints['1.pdbqt'])

5. For Tanimoto similarity:
similarity_dict = {}
for compound in df.Ligand_name:
  sim_val = DataStructs.TanimotoSimilarity(df.fingerprints[target],
                                           df.fingerprints[compound]
                                           )
  similarity_dict[compound] = sim_val
pd.Series(similarity_dict)
pd.Series(similarity_dict).sort_values()
pd.Series(similarity_dict).sort_values().plot.bar()
pd.Series(similarity_dict).sort_values().to_csv('Tanimoto_similarity.csv', index = False)


6. For USR similarity:
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

# Assuming you have a DataFrame 'df' with 'Ligand_name' as index and 'mols' column
# Also, 'target' is the name of the first compound you want to use as the target

# Extract the target compound name (the first one in the 'Ligand_name' column)
target = df.index[0]

# Convert the target molecule to RDKit molecule object
target_molecule = df['mols'][target]

# Generate 3D coordinates for the target molecule
AllChem.EmbedMolecule(target_molecule)

# Calculate the USR descriptor for the target molecule
target_usr = AllChem.GetUSR(target_molecule)

# Initialize a dictionary to store similarity values
similarity_dict_usr = {}

# Calculate USR similarity for each compound
for compound in df.index:
    if compound != target:
        # Convert the compound molecule to RDKit molecule object
        compound_molecule = df['mols'][compound]

        # Generate 3D coordinates for the compound molecule
        AllChem.EmbedMolecule(compound_molecule)

        # Calculate the USR descriptor for the compound
        compound_usr = AllChem.GetUSR(compound_molecule)

        # Calculate the USR similarity
        usr_similarity = AllChem.GetUSRScore(target_usr, compound_usr)
        similarity_dict_usr[compound] = usr_similarity

# Convert the similarity dictionary to a Pandas Series
similarity_series_usr = pd.Series(similarity_dict_usr)

# Print the similarity series
print(similarity_series_usr)
similarity_series_usr.sort_values()
similarity_series_usr.sort_values().plot.bar()

7. For VAM similarity:
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
from openeye import oechem, oeshape

# Assuming you have a DataFrame 'df' with 'Ligand_name' as index and 'mols' column
# Also, 'target' is the name of the first compound you want to use as the target

# Extract the target compound name (the first one in the 'Ligand_name' column)
target = df.index[0]

# Convert the target molecule to OpenEye molecule object
target_molecule = oechem.OEMol(df['mols'][target])

# Generate 3D coordinates for the target molecule
oechem.OE3DToInternalStereo(target_molecule)

# Initialize a dictionary to store similarity values
similarity_dict_shape = {}

# Initialize ROCS shape overlay options
rocs_options = oeshape.OEROCSOptions()
rocs_options.SetOverlayMethod(oeshape.OEROCSOverlayShape)

# Perform shape similarity comparison for each compound
for compound in df.index:
    if compound != target:
        # Convert the compound molecule to OpenEye molecule object
        compound_molecule = oechem.OEMol(df['mols'][compound])

        # Generate 3D coordinates for the compound molecule
        oechem.OE3DToInternalStereo(compound_molecule)

        # Perform ROCS shape overlay
        rocs_result = oeshape.OEROCSOverlay(rocs_options, target_molecule, compound_molecule)

        # Get the TanimotoCombo score
        tanimoto_similarity = rocs_result.GetTanimotoCombo()
        similarity_dict_shape[compound] = tanimoto_similarity

# Convert the similarity dictionary to a Pandas Series
similarity_series_shape = pd.Series(similarity_dict_shape)

# Print the similarity series
print(similarity_series_shape)

8. For MEP similarity:
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import pandas as pd

# Assuming you have a DataFrame 'df' with 'Ligand_name' as index and 'mols' column
# Also, 'target' is the name of the first compound you want to use as the target

# Extract the target compound name (the first one in the 'Ligand_name' column)
target = df.index[0]

# Convert the target molecule to RDKit molecule object
target_molecule = df['mols'][target]

# Ensure the target molecule is not None
if target_molecule is not None:
    # Generate 3D coordinates for the target molecule
    AllChem.EmbedMolecule(target_molecule)
    AllChem.MMFFOptimizeMolecule(target_molecule)

    # Calculate molecular descriptors based on electrostatic potential
    def calculate_mep_descriptor(mol):
        # Use AM1-BCC charges
        AllChem.ComputeGasteigerCharges(mol)
        AllChem.Compute2DCoords(mol)
        AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=False)
        return mol

    # Calculate the MEP descriptors for the target molecule
    target_molecule_with_charges = calculate_mep_descriptor(target_molecule)
    target_fingerprint = AllChem.GetMorganFingerprintAsBitVect(target_molecule_with_charges, 2)

    # Initialize a dictionary to store similarity values
    similarity_dict_mep = {}

    # Calculate Morgan fingerprint similarity for each compound
    for compound in df.index:
        if compound != target:
            compound_molecule = df['mols'][compound]

            # Ensure the compound molecule is not None
            if compound_molecule is not None:
                # Generate 3D coordinates for the compound molecule
                AllChem.EmbedMolecule(compound_molecule)
                AllChem.MMFFOptimizeMolecule(compound_molecule)

                # Calculate the MEP descriptors for the compound
                compound_molecule_with_charges = calculate_mep_descriptor(compound_molecule)

                # Calculate Morgan fingerprint for the compound
                compound_fingerprint = AllChem.GetMorganFingerprintAsBitVect(compound_molecule_with_charges, 2)

                # Calculate the Tanimoto similarity
                similarity = DataStructs.TanimotoSimilarity(target_fingerprint, compound_fingerprint)
                similarity_dict_mep[compound] = similarity

    # Convert the similarity dictionary to a Pandas Series
    similarity_series_mep = pd.Series(similarity_dict_mep)

    # Print the similarity series
    print(similarity_series_mep)
else:
    print("Target molecule is None.")
similarity_series_mep.sort_values()
similarity_series_mep.sort_values().plot.bar()
pd.Series(similarity_series_mep).sort_values().to_csv('MEP_similarity.csv', index = False)



9. For Tversky similarity:
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd
import os

# Assuming df has columns 'Filename' and 'SMILES'
input_csv = 'SMILES.csv'
df = pd.read_csv(input_csv)

# Generate fingerprints for all molecules in the DataFrame
df['Molecule'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x))
df['Fingerprint'] = df['Molecule'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2))

# Initialize similarity dictionary
similarity_dict6 = {}

# Target fingerprint (example, you can change this to your target molecule's fingerprint)
target_molecule = Chem.MolFromSmiles('CCO')  # Example: ethanol, replace with your target SMILES
target_fingerprint = AllChem.GetMorganFingerprintAsBitVect(target_molecule, 2)

# Compute Tversky similarity
alpha = 0.5
beta = 1 - alpha

for idx, row in df.iterrows():
    compound = row['Filename']
    fingerprint = row['Fingerprint']
    sim_val = DataStructs.TverskySimilarity(target_fingerprint, fingerprint, alpha, beta)
    similarity_dict6[compound] = sim_val

# Convert similarity dictionary to Series and process
similarity_series = pd.Series(similarity_dict6)
similarity_series_sorted = similarity_series.sort_values()

# Plot and save results
similarity_series_sorted.plot.bar()
similarity_series_sorted.to_csv('Tversky_similarity.csv', index=True)

print("Tversky similarity calculation completed. Results saved to Tversky_similarity.csv")

10. For Dice similarity:
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

# Step 1: Read the CSV file into a Pandas DataFrame
df = pd.read_csv('SMILES.csv')

# Convert SMILES strings to RDKit molecules
df['mols'] = df['SMILES'].apply(Chem.MolFromSmiles)

# Calculate the fingerprint from the molecule
df['fingerprints'] = df['mols'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2) if x is not None else None)

# Set the target molecule
target = '1.pdbqt'

# Ensure the target molecule is in the DataFrame
if target not in df['Filename'].values:
    raise ValueError(f"Target molecule '{target}' not found in DataFrame.")

# Set the target molecule as the index
df.set_index('Filename', inplace=True)

# Check if the target fingerprint is available
if target in df.index:
    target_fingerprint = df.loc[target, 'fingerprints']
else:
    raise ValueError(f"Target '{target}' not found in DataFrame index.")

# Initialize the similarity dictionary
similarity_dict3 = {'Filename': [], 'DiceSimilarity': []}

# Calculate similarities
for compound in df.index:
    if compound != target and df.loc[compound, 'fingerprints'] is not None:
        sim_val = DataStructs.DiceSimilarity(target_fingerprint, df.loc[compound, 'fingerprints'])
        similarity_dict3['Filename'].append(compound)
        similarity_dict3['DiceSimilarity'].append(sim_val)

# Convert dictionary to DataFrame
similarity_df3 = pd.DataFrame(similarity_dict3)

# Check if the similarity DataFrame is empty
if similarity_df3.empty:
    print("No similarity calculations were performed.")
else:
    # Sort the DataFrame by similarity values
    similarity_df3 = similarity_df3.sort_values(by='DiceSimilarity', ascending=False)
    
    # Save to CSV
    similarity_df3.to_csv('Dice_similarity_new.csv', index=False)
    
    print("Dice similarity calculation completed. Results saved to Dice_similarity_new.csv")
    
11. For Cosine similarity:
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

# Step 1: Read the CSV file into a Pandas DataFrame
df = pd.read_csv('SMILES.csv')

# Convert SMILES strings to RDKit molecules
df['mols'] = df['SMILES'].apply(Chem.MolFromSmiles)

# Calculate the fingerprint from the molecule
df['fingerprints'] = df['mols'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2) if x is not None else None)

# Set the target molecule
target = '1.pdbqt'

# Ensure the target molecule is in the DataFrame
if target not in df['Filename'].values:
    raise ValueError(f"Target molecule '{target}' not found in DataFrame.")

# Set the target molecule as the index
df.set_index('Filename', inplace=True)

# Check if the target fingerprint is available
if target in df.index:
    target_fingerprint = df.loc[target, 'fingerprints']
else:
    raise ValueError(f"Target '{target}' not found in DataFrame index.")

# Initialize the similarity dictionary
similarity_dict2 = {'Filename': [], 'CosineSimilarity': []}

# Calculate similarities
for compound in df.index:
    if compound != target and df.loc[compound, 'fingerprints'] is not None:
        sim_val = DataStructs.CosineSimilarity(target_fingerprint, df.loc[compound, 'fingerprints'])
        similarity_dict2['Filename'].append(compound)
        similarity_dict2['CosineSimilarity'].append(sim_val)

# Convert dictionary to DataFrame
similarity_df2 = pd.DataFrame(similarity_dict2)

# Check if the similarity DataFrame is empty
if similarity_df2.empty:
    print("No similarity calculations were performed.")
else:
    # Sort the DataFrame by similarity values
    similarity_df2 = similarity_df2.sort_values(by='CosineSimilarity', ascending=False)
    
    # Save to CSV
    similarity_df2.to_csv('Cosine_similarity_new.csv', index=False)
    
    print("Cosine similarity calculation completed. Results saved to Cosine_similarity_new.csv")
    
12. For Kulczynski similarity:
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

# Define Kulczynski similarity function
def KulczynskiSimilarity(fp1, fp2):
    intersection = DataStructs.BitVectToText(fp1 & fp2).count('1')
    a = DataStructs.BitVectToText(fp1).count('1')
    b = DataStructs.BitVectToText(fp2).count('1')
    if a == 0 or b == 0:
        return 0.0
    return (intersection / a + intersection / b) / 2

# Step 1: Read the CSV file into a Pandas DataFrame
df = pd.read_csv('SMILES.csv')

# Convert SMILES strings to RDKit molecules
df['mols'] = df['SMILES'].apply(Chem.MolFromSmiles)

# Calculate the fingerprint from the molecule
df['fingerprints'] = df['mols'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2) if x is not None else None)

# Set the target molecule
target = '1.pdbqt'

# Ensure the target molecule is in the DataFrame
if target not in df['Filename'].values:
    raise ValueError(f"Target molecule '{target}' not found in DataFrame.")

# Set the target molecule as the index
df.set_index('Filename', inplace=True)

# Check if the target fingerprint is available
if target in df.index:
    target_fingerprint = df.loc[target, 'fingerprints']
else:
    raise ValueError(f"Target '{target}' not found in DataFrame index.")

# Initialize the similarity dictionary
similarity_dict = {'Filename': [], 'KulczynskiSimilarity': []}

# Calculate similarities
for compound in df.index:
    if compound != target and df.loc[compound, 'fingerprints'] is not None:
        sim_val = KulczynskiSimilarity(target_fingerprint, df.loc[compound, 'fingerprints'])
        similarity_dict['Filename'].append(compound)
        similarity_dict['KulczynskiSimilarity'].append(sim_val)

# Convert dictionary to DataFrame
similarity_df = pd.DataFrame(similarity_dict)

# Check if the similarity DataFrame is empty
if similarity_df.empty:
    print("No similarity calculations were performed.")
else:
    # Sort the DataFrame by similarity values
    similarity_df = similarity_df.sort_values(by='KulczynskiSimilarity', ascending=False)
    
    # Save to CSV
    similarity_df.to_csv('Kulczynski_similarity.csv', index=False)
    
    print("Kulczynski similarity calculation completed. Results saved to Kulczynski_similarity.csv")
