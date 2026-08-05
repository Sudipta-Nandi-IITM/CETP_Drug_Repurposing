# residue_decomposition_heatmap_simple.py

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import TwoSlopeNorm

# =========================
# SETTINGS
# =========================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

DATA_GLOB = os.path.join(SCRIPT_DIR, "energyMapIn_*.dat")
FONT_PATH = os.path.join(SCRIPT_DIR, "timesbold.ttf")

OUT_PNG = os.path.join(SCRIPT_DIR, "Residue_energy_decomposition_heatmap.png")
OUT_PDF = os.path.join(SCRIPT_DIR, "Residue_energy_decomposition_heatmap.pdf")
OUT_CSV = os.path.join(SCRIPT_DIR, "Residue_energy_decomposition_matrix.csv")

TOP_N = 20                 # number of residues to display
FIG_DPI = 600
ANNOTATE = True
ANNOT_FONTSIZE = 8

# choose how residues are selected
# "absmean" = largest absolute contributions across ligands
# "most_negative" = most favorable residues
RESIDUE_SELECTION_MODE = "absmean"

# ligand name mapping from filename suffix
LIGAND_NAME_MAP = {
    "1": "Control",
    "2": "Compound1",
    "3": "Compound2",
    "4": "Compound3",
    "5": "Compound4",
    "6": "Compound5",
}

# preferred ligand order in plot
CUSTOM_LIGAND_ORDER = [
    "Control",
    "Compound1",
    "Compound2",
    "Compound3",
    "Compound4",
    "Compound5"
]

CMAP = "RdYlBu_r"   # negative favorable = blue, positive unfavorable = red
# =========================


def apply_font():
    if os.path.isfile(FONT_PATH):
        font_manager.fontManager.addfont(FONT_PATH)
        font_name = font_manager.FontProperties(fname=FONT_PATH).get_name()
        plt.rcParams["font.family"] = font_name

    plt.rcParams.update({
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 15,
        "xtick.labelsize": 10,
        "ytick.labelsize": 12,
        "axes.linewidth": 1.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def ligand_name_from_file(fp):
    base = os.path.basename(fp)
    # energyMapin_Tor.dat -> Tor
    short = base.replace("energyMapin_", "").replace(".dat", "")
    return LIGAND_NAME_MAP.get(short, short)


def parse_simple_two_column_file(fp):
    """
    Reads files of form:
    residue_number   contribution
    """
    data = np.loadtxt(fp)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 2:
        raise ValueError(f"{fp} does not have at least 2 columns.")

    residue_numbers = data[:, 0].astype(int)
    values = data[:, 1].astype(float)

    df = pd.DataFrame({
        "Residue": residue_numbers,
        "Energy": values
    })

    # if same residue number appears multiple times, sum as fallback
    df = df.groupby("Residue", as_index=False)["Energy"].sum()
    return df


def build_matrix(files):
    tables = []
    for fp in files:
        lig = ligand_name_from_file(fp)
        df = parse_simple_two_column_file(fp)
        df = df.rename(columns={"Energy": lig})
        tables.append(df.set_index("Residue"))

    mat = pd.concat(tables, axis=1).fillna(0.0)
    mat = mat.T  # rows = ligands, cols = residues

    # reorder ligands
    ordered = [x for x in CUSTOM_LIGAND_ORDER if x in mat.index]
    remaining = [x for x in mat.index if x not in ordered]
    mat = mat.loc[ordered + remaining]

    return mat


def select_residues(mat):
    if RESIDUE_SELECTION_MODE == "most_negative":
        score = mat.mean(axis=0).sort_values(ascending=True)
        residues = score.head(TOP_N).index.tolist()
    else:
        score = mat.abs().mean(axis=0).sort_values(ascending=False)
        residues = score.head(TOP_N).index.tolist()
    return residues


def plot_heatmap(mat, residues):
    plot_df = mat[residues].copy()

    # sort residues numerically for cleaner display
    plot_df = plot_df[sorted(plot_df.columns)]

    nrows, ncols = plot_df.shape
    fig_w = max(8.0, 0.45 * ncols + 3.5)
    fig_h = max(4.5, 0.70 * nrows + 2.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmax = np.nanmax(np.abs(plot_df.values))
    if vmax == 0:
        vmax = 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax.imshow(plot_df.values, cmap=CMAP, norm=norm, aspect="auto")

    ax.set_xticks(np.arange(ncols))
    ax.set_xticklabels([str(c) for c in plot_df.columns], rotation=45, ha="right", fontweight="bold")
    ax.set_yticks(np.arange(nrows))
    ax.set_yticklabels(plot_df.index, fontweight="bold")

    ax.set_xlabel("Residue number")
    ax.set_ylabel("Ligands")
    ax.set_title("Residue-wise binding free-energy decomposition", fontweight="bold", pad=12)

    # cell borders
    ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nrows, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    if ANNOTATE:
        for i in range(nrows):
            for j in range(ncols):
                val = plot_df.iloc[i, j]
                txt_color = "white" if abs(val) > 0.55 * vmax else "black"
                ax.text(j, i, f"{val:.2f}",
                        ha="center", va="center",
                        fontsize=ANNOT_FONTSIZE,
                        color=txt_color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label(r"Per-residue contribution (kcal/mol)", fontweight="bold")

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.0)

    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=FIG_DPI, bbox_inches="tight")
    plt.savefig(OUT_PDF, dpi=FIG_DPI, bbox_inches="tight")
    print(f"Saved: {OUT_PNG}")
    print(f"Saved: {OUT_PDF}")


def main():
    apply_font()

    files = sorted(glob.glob(DATA_GLOB))
    print("Searching pattern:", DATA_GLOB)
    print("Files found:")
    for f in files:
        print("  ", os.path.basename(f))

    if not files:
        raise FileNotFoundError(f"No files found with pattern: {DATA_GLOB}")

    mat = build_matrix(files)
    residues = select_residues(mat)

    mat[residues].to_csv(OUT_CSV)
    print(f"Saved matrix: {OUT_CSV}")

    plot_heatmap(mat, residues)


if __name__ == "__main__":
    main()
