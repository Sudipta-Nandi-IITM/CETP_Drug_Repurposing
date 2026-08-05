# figure3b_final_violin.py
# Output: Figure3B_final_violin.png

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
import matplotlib as mpl
from matplotlib.patches import Rectangle

# =========================
# SETTINGS
# =========================
SCREEN1_CSV = "Screen1_all.csv"   # columns: compound_id, score
SCREEN2_CSV = "Screen2_all.csv"   # columns: compound_id, score

ID_COL      = "compound_id"
S1_COL_RAW  = "score"
S2_COL_RAW  = "score"

LOWER_IS_BETTER = True
S1_CUTOFF = -9.0
S2_CUTOFF = -9.0

OUT_PNG = "Figure3B_final_violin.png"
DPI     = 600

# colors
GREY  = "#cfcfcf"
BLUE  = "#8DA0CB"
BLACK = "black"

# y-axis range similar to your original style
YLIM = (-16.2, -0.2)


# =========================
# FONT / STYLE
# =========================
mpl.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 18,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

# use timesbold.ttf if present
font_file = "timesbold.ttf"
if os.path.exists(font_file):
    font_manager.fontManager.addfont(font_file)
    font_name = font_manager.FontProperties(fname=font_file).get_name()
    plt.rcParams["font.family"] = font_name
    plt.rcParams["font.weight"] = "bold"


# =========================
# HELPERS
# =========================
def pass_mask(values, cutoff, lower_is_better=True):
    if lower_is_better:
        return values <= cutoff
    else:
        return values >= cutoff


def draw_continuous_colored_violin(ax, all_vals, cutoff, title,
                                   xtick_label_all, xtick_label_sel):
    """
    Draw one continuous violin using the full dataset shape.
    Grey part is only above cutoff.
    Blue part is only below cutoff.
    """

    all_vals = np.asarray(all_vals)
    pass_vals = all_vals[pass_mask(all_vals, cutoff, LOWER_IS_BETTER)]

    y_min, y_max = YLIM

    # ---- full violin in grey ----
    v_grey = ax.violinplot(
        [all_vals],
        positions=[1],
        widths=0.85,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )

    # clip grey to only ABOVE cutoff
    if LOWER_IS_BETTER:
        grey_clip = Rectangle(
            (0.0, cutoff),   # x, y
            2.0,             # width
            y_max - cutoff,  # height
            transform=ax.transData
        )
    else:
        grey_clip = Rectangle(
            (0.0, y_min),
            2.0,
            cutoff - y_min,
            transform=ax.transData
        )

    for body in v_grey["bodies"]:
        body.set_facecolor(GREY)
        body.set_edgecolor("none")
        body.set_alpha(1.0)
        body.set_clip_path(grey_clip)
        body.set_zorder(1)

    # ---- same full violin again in blue ----
    v_blue = ax.violinplot(
        [all_vals],
        positions=[1],
        widths=0.85,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )

    # clip blue to only BELOW cutoff
    if LOWER_IS_BETTER:
        blue_clip = Rectangle(
            (0.0, y_min),
            2.0,
            cutoff - y_min,
            transform=ax.transData
        )
    else:
        blue_clip = Rectangle(
            (0.0, cutoff),
            2.0,
            y_max - cutoff,
            transform=ax.transData
        )

    for body in v_blue["bodies"]:
        body.set_facecolor(BLUE)
        body.set_edgecolor("none")
        body.set_alpha(1.0)
        body.set_clip_path(blue_clip)
        body.set_zorder(2)

    # median of full dataset
    if len(all_vals):
        ax.scatter([1], [np.median(all_vals)], s=24, zorder=4, c=BLACK)

    # cutoff line
    ax.axhline(cutoff, linestyle="--", linewidth=1, color=BLACK, zorder=3)

    # x tick label
    ax.set_xticks([1], [f"{xtick_label_all}\n{xtick_label_sel}"])

    # title / axes
    ax.set_title(title)
    ax.set_xlim(0.5, 1.5)
    ax.set_ylim(YLIM)
    ax.grid(False)


# =========================
# LOAD DATA
# =========================
s1 = pd.read_csv(
    SCREEN1_CSV,
    usecols=[ID_COL, S1_COL_RAW]
).rename(columns={S1_COL_RAW: "S1"})

s2 = pd.read_csv(
    SCREEN2_CSV,
    usecols=[ID_COL, S2_COL_RAW]
).rename(columns={S2_COL_RAW: "S2"})

# compounds advanced to screen-2
adv = s1.merge(s2, on=ID_COL, how="inner")

# arrays
S1_all = s1["S1"].dropna().values
S2_all = adv["S2"].dropna().values

# counts for labels
S1_sel = S1_all[pass_mask(S1_all, S1_CUTOFF, LOWER_IS_BETTER)]
S2_sel = S2_all[pass_mask(S2_all, S2_CUTOFF, LOWER_IS_BETTER)]


# =========================
# PLOT
# =========================
fig, axes = plt.subplots(
    1, 2,
    figsize=(8.2, 4.4),   # like your 2nd figure
    sharey=True
)

draw_continuous_colored_violin(
    axes[0],
    S1_all,
    S1_CUTOFF,
    title="Screen-1 distribution",
    xtick_label_all=f"S1 all  (n={len(S1_all):,})",
    xtick_label_sel=f"pass cutoff (n={len(S1_sel):,})"
)

draw_continuous_colored_violin(
    axes[1],
    S2_all,
    S2_CUTOFF,
    title="Screen-2 distribution",
    xtick_label_all=f"S2 all  (n={len(S2_all):,})",
    xtick_label_sel=f"pass cutoff (n={len(S2_sel):,})"
)

axes[0].set_ylabel("Docking score (kcal/mol)")

# small gap between panels like your 2nd figure
plt.subplots_adjust(wspace=0.04, left=0.08, right=0.985, top=0.90, bottom=0.16)

# save
plt.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
plt.show()

print(f"Saved {OUT_PNG}")
