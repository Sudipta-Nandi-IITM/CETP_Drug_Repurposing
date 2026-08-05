import os
import pandas as pd

# Input and output directories
input_folder = '/media/sudipta/One Touch/Sudipta_phd/VS/output_folder'  # Replace with the path to your input folder
output_excel_file = 'Screen1_BP'  # Output Excel file name

# Initialize a list to store the results
results = []

for filename in os.listdir(input_folder):
    input_file = os.path.join(input_folder, filename)
    if os.path.isfile(input_file):
        try:
            # Read the file into a DataFrame
            df = pd.read_csv(input_file, delimiter='\t', header=None)

            # Extract the third column
            third_column = df.iloc[:, 2]

            # Count values less than 5
            count_gt_5 = (third_column < 5).sum()

            # Calculate the fraction of values less than 5
            total_values = len(third_column)
            fraction_gt_5 = count_gt_5 / total_values

            # Append results to the list
            results.append([filename, count_gt_5, total_values, fraction_gt_5])
        except Exception as e:
            print(f"Error processing file {filename}: {e}")

# Create a DataFrame for the results
results_df = pd.DataFrame(results, columns=['File', 'Count_GT_5', 'Total_Values', 'Fraction_GT_5'])

# Save the results to an Excel file
results_df.to_excel(output_excel_file, index=False)

print("Data written to Excel file:", output_excel_file)

