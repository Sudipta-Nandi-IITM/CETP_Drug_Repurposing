import os
import shutil

# Define the source folder (folder A) and destination folder (folder B)
source_folder = '/media/sudipta/One Touch/Sudipta_phd/VS/Screen1_inputs_26179_out'
destination_folder = '/media/sudipta/One Touch/Sudipta_phd/VS/Screen1_outputs_6296_out2'

# Read the list of file names from a text file
with open('Outputs.txt', 'r') as file:
    file_names = file.read().splitlines()

# Create the destination folder if it doesn't exist
if not os.path.exists(destination_folder):
    os.makedirs(destination_folder)

# Iterate through the file names and move the corresponding files
for file_name in file_names:
    source_file_path = os.path.join(source_folder, file_name)
    destination_file_path = os.path.join(destination_folder, file_name)

    # Check if the file exists in the source folder
    if os.path.exists(source_file_path):
        # Move the file to the destination folder
        shutil.move(source_file_path, destination_file_path)
        print(f"Moved '{file_name}' to '{destination_folder}'")
    else:
        print(f"File '{file_name}' not found in '{source_folder}'")

print("File separation complete.")

