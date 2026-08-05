# figure3a_scatter_with_marginals_rug.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# ---------- settings ----------
SCREEN1_CSV = "Screen1_all.csv"
SCREEN2_CSV = "Screen2_all.csv"
ID_COL      = "compound_id"
S1_COL_RAW  = "score"
S2_COL_RAW  = "score"

S1_CUTOFF = -9.0
S2_CUTOFF = -9.0

CUTOFF_COLOR    = "black"
POINT_COLOR     = "#2f5f8f"   # overlap scatter points
HIST_ALL_COLOR  = "#a1a1a1"   # full Screen1 distribution (light gray, background)
HIST_ADV_COLOR  = "#2f5f8f"   # subset that advanced to Screen2
HIST_S2_COLOR   = "#a1a1a1"   # right-hand Screen2 histogram
RUG_COLOR       = "#9a9a9a"   # Screen1-only rug ticks (real x, no fabricated y)

OUT_PNG = "Figure3A_scatter_marginals_rug.png"
DPI     = 600
# -----------------------------

try:
    font_manager.fontManager.addfont("timesbold.ttf")
    plt.rcParams["font.family"] = font_manager.FontProperties(fname="timesbold.ttf").get_name()
    plt.rcParams["font.weight"] = "bold"
except Exception as e:
    print(f"Font not loaded ({e}); using default font.")

# ---------- load & normalize IDs ----------
s1 = pd.read_csv(SCREEN1_CSV, usecols=[ID_COL, S1_COL_RAW]).rename(columns={S1_COL_RAW: "S1"})
s2 = pd.read_csv(SCREEN2_CSV, usecols=[ID_COL, S2_COL_RAW]).rename(columns={S2_COL_RAW: "S2"})
for d in (s1, s2):
    d[ID_COL] = d[ID_COL].astype(str).str.strip().str.upper()
s1 = s1.drop_duplicates(subset=ID_COL, keep="first")
s2 = s2.drop_duplicates(subset=ID_COL, keep="first")

df = s1.merge(s2, on=ID_COL, how="inner")  # overlap only, for the scatter + stats
print(f"n(S1)={len(s1):,}  n(S2)={len(s2):,}  overlap={len(df):,}")

# Screen1 compounds with no Screen2 value -- only their real S1 (x) value is
# known, so they get a rug (x-position only), never a fabricated y
s1_only = s1[~s1[ID_COL].isin(df[ID_COL])].copy()
print(f"Screen1-only (no S2 value) = {len(s1_only):,}")

pear  = df[["S1", "S2"]].corr(method="pearson").iloc[0, 1]
spear = df[["S1", "S2"]].corr(method="spearman").iloc[0, 1]

# ---------- layout ----------
fig = plt.figure(figsize=(7.0, 6.8), constrained_layout=True)
gs  = GridSpec(4, 4, figure=fig, wspace=0.05, hspace=0.05)

ax_scatter = fig.add_subplot(gs[1:4, 0:3])
ax_top     = fig.add_subplot(gs[0,    0:3], sharex=ax_scatter)
ax_right   = fig.add_subplot(gs[1:4, 3],   sharey=ax_scatter)

# ---- Real data scatter (overlap only) ----
ax_scatter.scatter(df["S1"], df["S2"], s=7, c=POINT_COLOR, alpha=0.65, edgecolors="none",
                    rasterized=True, zorder=3)

# ---- Extend the y-axis downward to make room for the rug band, so it never
#      overlaps real data points ----
y_data_min = df["S2"].min()
y_data_max = df["S2"].max()
y_range    = y_data_max - y_data_min
rug_band_h = 0.07 * y_range          # rug tick height
rug_gap    = 0.02 * y_range          # small gap between real data and rug
rug_base   = y_data_min - rug_gap - rug_band_h
new_ylim_bottom = rug_base - 0.01 * y_range

ax_scatter.set_ylim(new_ylim_bottom, y_data_max + 0.04 * y_range)

# ---- Rug: thin vertical ticks at each Screen1-only compound's real S1 value ----
ax_scatter.eventplot(s1_only["S1"].values, orientation="horizontal",
                      lineoffsets=rug_base + rug_band_h / 2, linelengths=rug_band_h,
                      colors=RUG_COLOR, alpha=0.35, linewidths=0.4, zorder=1)

# faint separator line between the real data region and the rug band
ax_scatter.axhline(y_data_min - rug_gap / 2, color="0.85", linewidth=0.8, zorder=2)

ax_scatter.axvline(S1_CUTOFF, linestyle="--", linewidth=1.2, color=CUTOFF_COLOR, zorder=4)
ax_scatter.axhline(S2_CUTOFF, linestyle="--", linewidth=1.2, color=CUTOFF_COLOR, zorder=4)
ax_scatter.set_xlabel("Docking score — Screen 1 (kcal/mol)", fontsize=16)
ax_scatter.set_ylabel("Docking score — Screen 2 (kcal/mol)", fontsize=16)
ax_scatter.tick_params(axis="both", labelsize=14)
ax_scatter.grid(False)

legend_handles = [
    Line2D([0], [0], marker="o", linestyle="", color=POINT_COLOR, markersize=6,
           alpha=0.8, label="In both screens"),
    Line2D([0], [0], marker="|", linestyle="", color=RUG_COLOR, markersize=10,
           markeredgewidth=1.5, label="Screen 1 only (S1 value shown, no S2)"),
]
ax_scatter.legend(handles=legend_handles, fontsize=8, loc="lower left", frameon=True)

# ---- Top marginal: FULL S1 distribution (gray) with the S2-advanced subset
#      overlaid in the scatter color on identical bins ----
bins_s1 = np.histogram_bin_edges(s1["S1"].dropna().values,
                                  bins=max(40, int(np.sqrt(len(s1)))))
ax_top.hist(s1["S1"].dropna().values, bins=bins_s1, color=HIST_ALL_COLOR,
            alpha=1.0, label="All Screen 1")
ax_top.hist(df["S1"].values, bins=bins_s1, color=HIST_ADV_COLOR,
            alpha=0.85, label="Advanced to Screen 2")
ax_top.axvline(S1_CUTOFF, linestyle="--", linewidth=1.0, color=CUTOFF_COLOR)
ax_top.axis("off")
ax_top.legend(handles=[Patch(fc=HIST_ALL_COLOR, label=f"All Screen 1 (n={len(s1):,})"),
                       Patch(fc=HIST_ADV_COLOR, label=f"Advanced to Screen 2 (n={len(df):,})")],
              fontsize=8, loc="upper left", frameon=False, handlelength=1.2)

# ---- Right marginal: Screen 2 distribution (unchanged) ----
bins_s2 = max(30, int(np.sqrt(len(s2))))
ax_right.hist(s2["S2"].dropna().values, bins=bins_s2, orientation="horizontal",
              color=HIST_S2_COLOR, alpha=0.9)
ax_right.axhline(S2_CUTOFF, linestyle="--", linewidth=1.0, color=CUTOFF_COLOR)
ax_right.axis("off")

# ---- Stats box ----
txt = (f"Pearson r = {pear:.2f}\nSpearman ρ = {spear:.2f}\n"
       f"n(S1)={len(s1):,}, n(S2)={len(s2):,}, overlap={len(df):,}")
ax_scatter.text(0.45, 0.98, txt, transform=ax_scatter.transAxes, va="top", ha="left", fontsize=10,
                bbox=dict(boxstyle="round,pad=0.35", fc="#e9f2ff", ec="none"))

plt.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
print(f"Saved {OUT_PNG}")
