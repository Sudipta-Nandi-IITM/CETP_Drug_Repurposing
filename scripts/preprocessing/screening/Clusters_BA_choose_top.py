import pandas as pd
import os

# Directory containing the cluster files
cluster_dir = '/media/sudipta/e48a0503-e8e3-4ac5-92e7-7c98a3d9f7f3/SUDIPTA/vs/CETP/Similarities/individual_clusters'  # Change this to your actual path where cluster files are stored

# File containing the binding affinities
affinity_file = '/media/sudipta/e48a0503-e8e3-4ac5-92e7-7c98a3d9f7f3/SUDIPTA/vs/CETP/Similarities/output.csv'  # Change this to your actual path

# Read the binding affinities into a DataFrame
affinity_df = pd.read_csv(affinity_file)

# Strip leading/trailing spaces from compound names and values
affinity_df['Compound'] = affinity_df['Compound'].str.strip()
affinity_df['Binding_Affinity'] = affinity_df['Binding_Affinity'].astype(str).str.strip()

# Initialize a list to hold the results
results = []

# Iterate through the cluster files
for i in range(1, 10):
    cluster_file = os.path.join(cluster_dir, f'c{i}.txt')
    
    # Read the cluster file
    with open(cluster_file, 'r') as file:
        compounds = file.read().splitlines()
    
    # Strip leading/trailing spaces from compound names
    compounds = [compound.strip() for compound in compounds]
    
    # Find the binding affinities for the compounds in this cluster
    for compound in compounds:
        affinity = affinity_df[affinity_df['Compound'].str.lower() == compound.lower()]['Binding_Affinity'].values
        if affinity.size > 0:
            results.append({'Cluster': f'c{i}', 'Compound': compound, 'Binding_Affinity': affinity[0]})
        else:
            print(f"Compound {compound} not found in affinity file.")

# Convert the results to a DataFrame
results_df = pd.DataFrame(results)

# Save the results to a CSV file
results_df.to_csv('cluster_binding_affinities_basedon1.csv', index=False)

# Print the results to verify
print(results_df)
