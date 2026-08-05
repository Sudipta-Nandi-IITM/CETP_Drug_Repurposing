#!/bin/bash

input_folder="/path/to/the/input/files"
output_folder="/path/to/the/output/files"
script_path="/path/to/the/file/with/process_pdbqt.sh"  # Replace with the actual path to your script

for file in "$input_folder"/*.pdbqt; do
    if [ -e "$file" ]; then
        output_file="$output_folder/$(basename "$file" .pdbqt)_processed.pdbqt"

        # Execute the script file
        bash "$script_path" "$file" "$output_file"
    fi
done


