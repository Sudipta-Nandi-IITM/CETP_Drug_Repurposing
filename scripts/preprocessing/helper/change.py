import os

# Define the path to the folder containing the files
folder_path = 'compounds'

# Function to process each file
def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    file_name = os.path.basename(file_path)
    updated_lines = []
    
    for line in lines:
        if "REMARK" in line and "Name" in line:
            parts = line.split("=")
            if len(parts) > 1:
                updated_line = f"{parts[0].strip()} = {file_name}\n"
            else:
                updated_line = f"{line.strip()} = {file_name}\n"
            updated_lines.append(updated_line)
        else:
            updated_lines.append(line)
    
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

# Get a list of all files in the folder
file_list = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]

# Process each file
for file_name in file_list:
    file_path = os.path.join(folder_path, file_name)
    process_file(file_path)

print("Files have been updated.")

