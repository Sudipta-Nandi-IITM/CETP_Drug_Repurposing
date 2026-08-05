import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
from matplotlib import font_manager

FEL_FILE = "Gibbs.txt"
PC1_FILE = "energy_PC1.xvg"
PC2_FILE = "energy_PC2.xvg"
FONT_FILE = "timesbold.ttf"
OUT_FILE = "FEL.png"

SYSTEM_NAME = "Control_Protein"

font_manager.fontManager.addfont(FONT_FILE)
timesbold_name = font_manager.FontProperties(fname=FONT_FILE).get_name()

plt.rcParams["font.family"] = timesbold_name
plt.rcParams["font.size"] = 16
plt.rcParams["axes.titlesize"] = 18
plt.rcParams["axes.labelsize"] = 16
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16
plt.rcParams["axes.linewidth"] = 2.0

def load_xvg(filename):
    data = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("@") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)

x, y, z = np.loadtxt(FEL_FILE).T

pc1_data = load_xvg(PC1_FILE)
pc2_data = load_xvg(PC2_FILE)

pc1_x, pc1_y = pc1_data[:, 0], pc1_data[:, 1]
pc2_x, pc2_y = pc2_data[:, 0], pc2_data[:, 1]

NGRID = 200

xi = np.linspace(np.min(x), np.max(x), NGRID)
yi = np.linspace(np.min(y), np.max(y), NGRID)
XI, YI = np.meshgrid(xi, yi)
ZI = griddata((x, y), z, (XI, YI), method="linear")
ZI_linear = griddata((x, y), z, (XI, YI), method="linear")
ZI_nearest = griddata((x, y), z, (XI, YI), method="nearest")

ZI[np.isnan(ZI)] = ZI_linear[np.isnan(ZI)]
ZI[np.isnan(ZI)] = ZI_nearest[np.isnan(ZI)]

fig = plt.figure(figsize=(13, 7))
gs = GridSpec(2, 2, width_ratios=[1.0, 1.35], height_ratios=[1, 1],
              wspace=0.30, hspace=0.30)

ax_pc1 = fig.add_subplot(gs[0, 0])
ax_pc2 = fig.add_subplot(gs[1, 0])
ax_fel = fig.add_subplot(gs[:, 1])

line_color = "#B13A1B"   
cmap = plt.get_cmap("RdYlBu_r").copy()
cmap.set_bad(cmap(0.0))

ax_pc1.plot(pc1_x, pc1_y, color=line_color, linewidth=2.2)
ax_pc1.set_xlabel("PC1")
ax_pc1.set_ylabel("ΔG (kcal/mol)")

ax_pc1.grid(False)

ax_pc1.text(-0.12, 1.05, "(a)", transform=ax_pc1.transAxes,
            fontsize=16, fontweight="bold", va="bottom")
ax_pc2.plot(pc2_x, pc2_y, color=line_color, linewidth=2.2)
ax_pc2.set_xlabel("PC2")
ax_pc2.set_ylabel("ΔG (kcal/mol)")
ax_pc2.grid(False)

ax_pc2.text(-0.12, 1.05, "(b)", transform=ax_pc2.transAxes,
            fontsize=16, fontweight="bold", va="bottom")

for ax in [ax_pc1, ax_pc2]:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.2)  
        spine.set_color("black")   
levels = np.arange(0, 13, 1)

cf = ax_fel.contourf(XI, YI, ZI, levels=levels, cmap=cmap, extend="max")
ax_fel.contour(XI, YI, ZI, levels=levels, colors="k", linewidths=0.2, alpha=0.25)

ax_fel.set_xlabel("PC1")
ax_fel.set_ylabel("PC2")

ax_fel.text(-0.08, 1.02, "(c)", transform=ax_fel.transAxes,
            fontsize=18, fontweight="bold", va="bottom")

cbar = fig.colorbar(cf, ax=ax_fel, fraction=0.046, pad=0.04)
cbar.set_label("Free energy (kcal/mol)", fontsize=16)

plt.tight_layout()
plt.savefig(OUT_FILE, dpi=600, bbox_inches="tight")
plt.show()
