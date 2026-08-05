import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO

# List of file paths
file_paths = ["numcontacts_1.xvg", "numcontacts_2.xvg", "numcontacts_3.xvg", "numcontacts_4.xvg", "numcontacts_5.xvg", "numcontacts_6.xvg"]  # Add more file paths as needed

# Define colors for smoothed plots
smoothed_colors = ["black", "darkblue", "darkgreen", "darkred", "purple", "orange"]  # Extend if more files
raw_color = "gray"  # Same color for raw data

# Initialize the plot
plt.figure(figsize=(14, 7))

# Loop through each file and process it
for i, file_path in enumerate(file_paths):
    # Pre-process the file to filter out lines starting with # or @
    with open(file_path, 'r') as file:
        lines = [line for line in file if not (line.startswith('#') or line.startswith('@'))]
    
    # Load the cleaned data into a DataFrame
    data = pd.read_csv(StringIO(''.join(lines)), delim_whitespace=True, header=None)
    data.columns = ["Time (ps)", "Number of Contacts"]
    
    # Convert time from ps to ns
    data["Time (ns)"] = data["Time (ps)"] / 1000

 # Calculate the average number of contacts
    avg_contacts = data["Number of Contacts"].mean()


    # Smooth the data using a moving average
    window_size = 50  # Adjust based on desired smoothness
    data["Smoothed Contacts"] = data["Number of Contacts"].rolling(window=window_size, center=True).mean()

    # Plot raw data
   # plt.plot(data["Time (ns)"], data["Number of Contacts"], label=f"Raw Contacts ({file_path})", alpha=0.3, color=raw_color, linewidth=1)
    
     # Plot smoothed data
    plt.plot(
        data["Time (ns)"], 
        data["Smoothed Contacts"], 
        label=f"{file_path} (Avg: {avg_contacts:.1f})" # Add average contacts to the label
        color=smoothed_colors[i % len(smoothed_colors)], 
        linewidth=2
    )
# Title and labels
plt.rcParams['font.family'] = 'Serif'
plt.title("Protein-Ligand Contact Stability Analysis", fontsize=16, fontweight='bold')
plt.xlabel("Time (ns)", fontsize=14, fontweight='bold')
plt.ylabel("Number of Contacts", fontsize=14, fontweight='bold')

# Legend
plt.legend(
    fontsize=12, 
    loc="upper left",  # Initial corner location
    bbox_to_anchor=(1.05, 1),  # Move the legend outside the plot area
    borderaxespad=0,  # Remove padding
    frameon=False  # Make the legend background transparent
)
# Customization
plt.grid(False)
plt.tight_layout()

plt.savefig("Contacts_all.png", dpi=300, bbox_inches="tight")

# Show the plot
plt.show()
