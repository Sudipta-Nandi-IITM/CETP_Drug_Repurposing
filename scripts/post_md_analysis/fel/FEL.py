import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Set font and style
plt.rcParams['font.family'] = ['Times New Roman', 'serif']
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

# Load data (Gibbs.txt or equivalent FES_PC1_PC2.dat)
x, y, z = np.loadtxt('Gibbs.txt').T  # Replace with the path to your data file

# Define grid for interpolation
xi = np.linspace(min(x), max(x), 100)
yi = np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')

# Multiplot layout
fig, axes = plt.subplots(1, 3, figsize=(24, 8))

# --- Plot 1: Free Energy Landscape ---
ax = axes[0]
cset = ax.contourf(xi, yi, zi, levels=np.arange(0, 13, 1), cmap='jet', extend='both')
ax.set_title("Free Energy Landscape", fontsize=16, fontweight='bold')
ax.set_xlabel("PC1", fontweight='bold')
ax.set_ylabel("PC2", fontweight='bold')
cbar = fig.colorbar(cset, ax=ax, shrink=0.8, aspect=20)
cbar.set_label("Free Energy (kcal/mol)", fontsize=12)

# --- Plot 2: Energy Profile along PC1 ---
pc1_data = np.loadtxt('energy_PC1.xvg')  # Replace with your PC1 energy file
pc1_x = pc1_data[:, 0]
pc1_y = pc1_data[:, 1]
axes[1].plot(pc1_x, pc1_y, color='red', linewidth=2)
axes[1].set_title("Energy Profile along PC1", fontsize=16)
axes[1].set_xlabel("PC1")
axes[1].set_ylabel(r'$\Delta G$ (kcal/mol)')
axes[1].grid(True)

# --- Plot 3: Energy Profile along PC2 ---
pc2_data = np.loadtxt('energy_PC2.xvg')  # Replace with your PC2 energy file
pc2_x = pc2_data[:, 0]
pc2_y = pc2_data[:, 1]
axes[2].plot(pc2_x, pc2_y, color='red', linewidth=2)
axes[2].set_title("Energy Profile along PC2", fontsize=16)
axes[2].set_xlabel("PC2")
axes[2].set_ylabel(r'$\Delta G$ (kcal/mol)')
axes[2].grid(True)

# Adjust layout
plt.tight_layout()

# Save the figure
plt.savefig('FEL_C_alpha.png', dpi=300, bbox_inches='tight')

# Show the plots
plt.show()
