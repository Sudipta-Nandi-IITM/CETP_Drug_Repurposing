import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  # For kernel density estimation
import matplotlib.pyplot as plt
from matplotlib import font_manager

# List of your rg.xvg file paths (replace with actual paths)
file_paths = [
    "rg_1.xvg",
    "rg_2.xvg",
    "rg_3.xvg",
    "rg_4.xvg",
    "rg_5.xvg",
    "rg_6.xvg"
]

# Set up seaborn for nice plotting style, using a dark color palette
sns.set(style="white", palette="dark")

plt.figure(figsize=(10, 6))
font_manager.fontManager.addfont("timesbold.ttf")
plt.rcParams["font.family"] = font_manager.FontProperties(fname="timesbold.ttf").get_name()
plt.rcParams["font.weight"] = "bold"
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

# Loop through each file to read and plot the Rg data
for file_path in file_paths:
    try:
        # Read the file, skip comment lines starting with # and @
        data = pd.read_csv(file_path, comment='@', sep='\s+', header=None, usecols=[0, 1])

        # Ensure the data is numeric, ignore rows where 'Radius of Gyration' is non-numeric
        data.columns = ["Time (ps)", "Radius of Gyration (nm)"]
        
        # Convert 'Radius of Gyration (nm)' to numeric, coercing errors to NaN, and drop NaNs
        data["Radius of Gyration (nm)"] = pd.to_numeric(data["Radius of Gyration (nm)"], errors='coerce')
        data = data.dropna(subset=["Radius of Gyration (nm)"])  # Drop rows with NaN values
        
        # Extract the Rg values
        rg_values = data["Radius of Gyration (nm)"]
        
        # Plot the probability density function (KDE) for each Rg dataset
        sns.kdeplot(rg_values, label=file_path.split('.')[0], fill=False, alpha=0.8, linewidth=2.5)

    except Exception as e:
        print(f"Error reading {file_path}: {e}")

plt.title("", fontsize=14, fontweight='bold')
plt.xlabel("Radius of Gyration (nm)", fontsize=16, fontweight='bold')
plt.ylabel("Probability Density", fontsize=16, fontweight='bold')
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, MaxNLocator

ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(0.05))   
ax.yaxis.set_major_locator(MultipleLocator(2.0))    


ax.xaxis.set_minor_locator(AutoMinorLocator(5))   
ax.yaxis.set_minor_locator(AutoMinorLocator(4))

ax.tick_params(axis='both', which='major', length=8, width=1.2, direction='in',
               bottom=True, left=True, labelsize=14)   # show ticks on all sides

for spine in ["bottom", "left"]:
    ax.spines[spine].set_linewidth(2.0)   # thick axis box lines
for spine in ["top", "right"]:
    ax.spines[spine].set_linewidth(2.0)   # hide opposite box lines

plt.legend(title="Simulations", fontsize=14, title_fontsize=14, loc='best', frameon=False, prop={'weight': 'bold'})
plt.xticks()
plt.yticks()
plt.grid(False)

plt.tight_layout()

plt.savefig("Rg_Probability_Distribution2.png")
plt.show()
