import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path

# ===================== SETTINGS =====================
CSV = "flex_strain_merged.csv"

ID_COL   = "base_id"      # compound identifier
AFF_COL  = "affinity"
STR_COL  = "strain"
MODE_COL = "mode"

OUT = "Fig4B_boxplot_minTEU_per_compound.png"
DPI = 600

BINS = 10
STRAIN_THRESHOLD = 5.0

FONT_PATH = "timesbold.ttf"

FIGSIZE = (12.0, 6.8)   # slightly taller to make room for n= row below axis

FS_AXIS   = 20
FS_TICK   = 16
FS_ANNOT  = 16
FS_N      = 12

BOX_COLOR    = "#AFC6DA"   # muted publication teal
EDGE_COLOR   = "#1A3A5C"
MEDIAN_COLOR = "#1A1A1A"
SHADE_COLOR  = "#D3D3D3"
THRESH_COLOR = "#2A2A2A"
CAP_COLOR    = "#1A1A1A"   # neutral -- same family as median/edge, not an
                           # accent color, so it doesn't compete visually

SHOW_FLIERS = False
# ====================================================

# ---------- font ----------
if Path(FONT_PATH).is_file():
    font_manager.fontManager.addfont(FONT_PATH)
    plt.rcParams["font.family"] = font_manager.FontProperties(fname=FONT_PATH).get_name()

plt.rcParams.update({
    "font.weight": "bold",
    "axes.labelsize": FS_AXIS,
    "xtick.labelsize": FS_TICK,
    "ytick.labelsize": FS_TICK,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.transparent": False,
})

# ---------- load ----------
df = pd.read_csv(CSV).dropna(subset=[AFF_COL, STR_COL]).copy()

# ====================================================
# STEP 1: choose one representative pose per compound
# ====================================================
df_rep = (
    df.sort_values([ID_COL, STR_COL, AFF_COL], ascending=[True, True, True])
      .groupby(ID_COL, as_index=False)
      .first()
)

print(f"Total compounds after min-TEU selection: {len(df_rep)}")

# ---------- affinity binning ----------
df_rep["aff_bin"] = pd.qcut(df_rep[AFF_COL], BINS, duplicates="drop")

order = sorted(df_rep["aff_bin"].unique(), key=lambda x: x.left)
labels = [f"{b.left:.1f} to {b.right:.1f}" for b in order]

plot_data = [df_rep.loc[df_rep["aff_bin"] == b, STR_COL].values for b in order]

# ---------- pass percentage + n per bin ----------
pass_pct = []
n_per_bin = []
for b in order:
    sub = df_rep.loc[df_rep["aff_bin"] == b, STR_COL]
    pct = 100.0 * (sub <= STRAIN_THRESHOLD).sum() / len(sub)
    pass_pct.append(pct)
    n_per_bin.append(len(sub))

# ---------- y-axis limits ----------
ymax = np.percentile(df_rep[STR_COL], 98)
ymax = max(ymax, STRAIN_THRESHOLD + 5)
ymax = np.ceil(ymax)

# Small padding below the true data floor (strain >= 0) so whisker caps
# landing at/near 0 don't visually fuse with the axis spine.
y_bottom_pad = 0.05 * ymax
ymin_display = -y_bottom_pad

# ---------- plot ----------
fig, ax = plt.subplots(figsize=FIGSIZE)

ax.axhspan(0, STRAIN_THRESHOLD, color=SHADE_COLOR, alpha=1.0, zorder=0)
ax.axhline(STRAIN_THRESHOLD, linestyle="--", linewidth=1.5, color=THRESH_COLOR)
ax.axhline(0, color="0.6", linewidth=0.8, zorder=0.5)

bp = ax.boxplot(
    plot_data,
    patch_artist=True,
    showfliers=SHOW_FLIERS,
    widths=0.8
)

for patch in bp["boxes"]:
    patch.set_facecolor(BOX_COLOR)
    patch.set_edgecolor(EDGE_COLOR)
    patch.set_linewidth(1.2)

for whisk in bp["whiskers"]:
    whisk.set_color(EDGE_COLOR)
    whisk.set_linewidth(2.0)

# Caps: neutral color, just thicker than the whisker itself so the "end"
# is unambiguous without introducing a competing accent color
for cap in bp["caps"]:
    cap.set_color(CAP_COLOR)
    cap.set_linewidth(2.4)
    cap.set_solid_capstyle("butt")

for med in bp["medians"]:
    med.set_color(MEDIAN_COLOR)
    med.set_linewidth(2.0)

# annotate % below threshold (inside the boxes, unchanged)
for i, pct in enumerate(pass_pct, start=1):
    ax.text(i, STRAIN_THRESHOLD + 0.7, f"{pct:.0f}%", ha="center", va="bottom", fontsize=FS_ANNOT)

# axes
ax.set_xticks(range(1, len(order) + 1))
ax.set_xticklabels(labels, rotation=24, ha="right")
ax.set_xlabel("Affinity bin (kcal/mol)")
ax.set_ylabel("Ligand strain (TEU)")
ax.set_ylim(ymin_display, ymax)

ax.grid(False)
ax.spines["top"].set_visible(True)
ax.spines["right"].set_visible(True)

# n= labels: placed BELOW the axis (outside the data area entirely) using
# the x-axis transform, so they never compete with whisker caps or the
# shaded threshold region
#for i, n in enumerate(n_per_bin, start=1):
#    ax.text(i, -0.04, f"n={n}", transform=ax.get_xaxis_transform(),
#            ha="center", va="top", fontsize=FS_N, color="0.35", clip_on=False)

plt.subplots_adjust(bottom=0.32, left=0.11, right=0.98, top=0.96)

plt.savefig(OUT, dpi=DPI, bbox_inches="tight")
plt.savefig(OUT.replace(".png", ".pdf"), bbox_inches="tight")
print("Saved", OUT)
