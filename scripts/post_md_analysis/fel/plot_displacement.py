import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm

INPUT_CSV = "pc1_displacement.csv"
OUTPUT_PLOT = "PC1_displacement_plot.png"

residues = []
values = []

with open(INPUT_CSV, newline="") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        residues.append(int(row["residue_number"]))
        values.append(float(row["displacement_A"]))

cmap = cm.get_cmap("viridis")
max_v = max(values)
colors = [cmap(v / max_v) for v in values]

fig, ax = plt.subplots(figsize=(6, 3))
ax.bar(residues, values, color=colors, width=1.0)
ax.set_xlabel("Residue number")
ax.set_ylabel(r"C$\alpha$ displacement along PC1 (Å)")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(OUTPUT_PLOT, dpi=300)
print(f"Saved {OUTPUT_PLOT}")
