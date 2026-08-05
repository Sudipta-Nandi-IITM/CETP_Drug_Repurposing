import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path

# =========================
# USER SETTINGS
# =========================
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()

# Input files
DATA_GLOB = os.path.join(SCRIPT_DIR, "energyMapIn*.dat")
# If your files have different capitalization, use:
# DATA_GLOB = os.path.join(SCRIPT_DIR, "energyMap*.dat")

# PDB file
PDB_FILE = os.path.join(SCRIPT_DIR, "Protein.pdb")

# Output files
OUT_PNG = os.path.join(SCRIPT_DIR, "Residue_Decomposition_Heatmap.png")
OUT_PDF = os.path.join(SCRIPT_DIR, "Residue_Decomposition_Heatmap.pdf")
OUT_SVG = os.path.join(SCRIPT_DIR, "Residue_Decomposition_Heatmap.svg")

FONT_PATH = os.path.join(SCRIPT_DIR, "timesbold.ttf")

KEY_RESIDUES = [
    13, 15, 129, 136, 195, 197, 198, 199, 202,
    205, 215, 217, 221, 228, 232, 261, 263, 441, 468
]

AUTO_SELECT = False
TOP_N_AUTO = 20

# If your PDB has multiple chains, set the chain here.
# Example: CHAIN_ID = "A"
# If None, the first matching residue number is used.
CHAIN_ID = None

CUSTOM_LIGAND_ORDER = [
    "Control",
    "Compound1",
    "Compound2",
    "Compound3",
    "Compound4",
    "Compound5"
]

# Appearance
FIGSIZE = (17, 7)
DPI = 600

FS_BASE = 14
FS_AXIS = 18
FS_TICK = 14
FS_TITLE = 18
FS_CBAR = 18
FS_CELL = 12

CMAP = "RdYlBu_r"

SHOW_CELL_VALUES = True
VALUE_FMT = "{:.2f}"

ROUND_DECIMALS = 1
# =========================


def apply_font():
    if Path(FONT_PATH).is_file():
        font_manager.fontManager.addfont(FONT_PATH)
        family = font_manager.FontProperties(fname=FONT_PATH).get_name()
        plt.rcParams["font.family"] = family

    plt.rcParams.update({
        "font.size": FS_BASE,
        "axes.labelsize": FS_AXIS,
        "axes.titlesize": FS_TITLE,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "axes.linewidth": 1.5,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.transparent": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def ligand_from_filename(filepath):
    base = os.path.basename(filepath).lower()

    filename_map = {
        "1": "Control",
        "2": "Compound1",
        "3": "Compound2",
        "4": "Compound3",
        "5": "Compound4",
        "6": "Compound5",
    }

    for key, ligand in filename_map.items():
        if key in base:
            return ligand

    raise ValueError(f"Could not determine ligand name from filename: {filepath}")


def load_energy_file(filepath):

    df = pd.read_csv(
        filepath,
        sep=r"\s+",
        header=None,
        names=["residue", "energy"],
        engine="python",
        comment="#"
    )

    df["residue"] = df["residue"].astype(int)
    df["energy"] = df["energy"].astype(float)

    return df


def standardize_resname(resname):
    resname = resname.strip().upper()

    residue_fix = {
        "HID": "HIS",
        "HIE": "HIS",
        "HIP": "HIS",
        "HSD": "HIS",
        "HSE": "HIS",
        "HSP": "HIS",
        "CYX": "CYS",
        "CYM": "CYS",
        "ASH": "ASP",
        "GLH": "GLU",
        "LYN": "LYS",
        "MSE": "MET",
    }

    return residue_fix.get(resname, resname)


def parse_pdb_line_fixed_width(line):
    """
    Standard PDB fixed-width parser.
    """
    record = line[0:6].strip()

    if record not in ["ATOM", "HETATM"]:
        return None

    resname = line[17:20].strip()
    chain = line[21].strip()
    resnum_text = line[22:26].strip()

    if not resname or not resnum_text:
        return None

    try:
        resnum = int(resnum_text)
    except ValueError:
        return None

    return chain, resnum, standardize_resname(resname)


def parse_pdb_line_split(line):
   
    parts = line.split()

    if len(parts) < 5:
        return None

    if parts[0] not in ["ATOM", "HETATM"]:
        return None

    # Common form:
    # ATOM serial atom resname chain resnum ...
    if len(parts) >= 6:
        resname = parts[3]

        try:
            resnum = int(parts[5])
            chain = parts[4]
            return chain, resnum, standardize_resname(resname)
        except ValueError:
            pass

    # Form without chain:
    # ATOM serial atom resname resnum ...
    try:
        resname = parts[3]
        resnum = int(parts[4])
        chain = ""
        return chain, resnum, standardize_resname(resname)
    except ValueError:
        return None


def read_pdb_residue_map(pdb_file, chain_id=None):
    
    if not os.path.isfile(pdb_file):
        print(f"\nWARNING: protein.pdb not found at: {pdb_file}")
        print("Trying to find any .pdb file in the same folder...")

        pdb_candidates = sorted(glob.glob(os.path.join(SCRIPT_DIR, "*.pdb")))

        if len(pdb_candidates) == 0:
            print("No PDB file found. Residue labels will remain numbers only.")
            return {}

        pdb_file = pdb_candidates[0]
        print(f"Using PDB file: {os.path.basename(pdb_file)}")

    residue_map = {}

    with open(pdb_file, "r") as f:
        for line in f:
            parsed = parse_pdb_line_fixed_width(line)

            if parsed is None:
                parsed = parse_pdb_line_split(line)

            if parsed is None:
                continue

            chain, resnum, resname = parsed

            if chain_id is not None:
                if chain != chain_id:
                    continue

            if resnum not in residue_map:
                residue_map[resnum] = f"{resname}{resnum}"

    print(f"\nPDB file used: {os.path.basename(pdb_file)}")
    print(f"Total residues read from PDB: {len(residue_map)}")

    return residue_map


def build_matrix(files):
    all_data = {}

    for fp in files:
        ligand = ligand_from_filename(fp)
        print(f"Parsed file {os.path.basename(fp)} -> ligand='{ligand}'")

        df = load_energy_file(fp)
        all_data[ligand] = df.set_index("residue")["energy"]

    mat = pd.DataFrame(all_data)
    return mat


def auto_select_residues(mat, top_n=20):
    mean_abs = mat.abs().mean(axis=1).sort_values(ascending=False)
    return sorted(mean_abs.head(top_n).index.tolist())


def get_residue_labels(residue_numbers, residue_map):
    labels = []

    for r in residue_numbers:
        r = int(r)
        labels.append(residue_map.get(r, str(r)))

    return labels


def round_down(x, decimals=1):
    factor = 10 ** decimals
    return np.floor(x * factor) / factor


def round_up(x, decimals=1):
    factor = 10 ** decimals
    return np.ceil(x * factor) / factor


def choose_text_color(value, vmin, vmax):
    """
    White text for dark cells, black text for light cells.
    """
    span = max(abs(vmin), abs(vmax))

    if abs(value) >= 0.45 * span:
        return "white"
    else:
        return "black"


def make_colorbar_ticks(vmin, vmax):

    ticks = []

    # Negative ticks
    if vmin < 0:
        # if vmin is exactly on a 2-unit grid like -10.0, -8.0 etc.
        if np.isclose(vmin % 2, 0, atol=1e-8) or np.isclose((vmin % 2), 2, atol=1e-8):
            start = int(vmin)
        else:
            # move to next 2-step tick above vmin
            start = int(np.ceil(vmin / 2.0) * 2 + 2)

        neg_ticks = list(np.arange(start, 0, 2))
        ticks.extend(neg_ticks)

    # Zero
    if vmin < 0 < vmax:
        ticks.append(0)

    # Positive ticks
    if vmax > 0:
        if vmax <= 4:
            pos_ticks = list(np.arange(1, int(np.floor(vmax)) + 1, 1))
        else:
            pos_ticks = list(np.arange(2, int(np.floor(vmax)) + 1, 2))
        ticks.extend(pos_ticks)

    # Add actual limits
    ticks = [t for t in ticks if vmin <= t <= vmax]

    if not any(np.isclose(t, vmin) for t in ticks):
        ticks.insert(0, vmin)

    if not any(np.isclose(t, vmax) for t in ticks):
        ticks.append(vmax)

    return ticks

def main():
    apply_font()

    residue_map = read_pdb_residue_map(PDB_FILE, CHAIN_ID)

    print("\nResidue labels for selected residues:")
    for r in KEY_RESIDUES:
        print(f"  {r} -> {residue_map.get(r, str(r))}")

    files = sorted(glob.glob(DATA_GLOB))

    print("\nSearching pattern:", DATA_GLOB)
    print("Files found:")
    for f in files:
        print("  ", os.path.basename(f))

    if not files:
        raise FileNotFoundError(f"No files found matching: {DATA_GLOB}")

    mat = build_matrix(files)

    print("\nMatrix columns before ordering:")
    print(list(mat.columns))

    if AUTO_SELECT:
        residues = auto_select_residues(mat, TOP_N_AUTO)
    else:
        residues = [r for r in KEY_RESIDUES if r in mat.index]

    if len(residues) == 0:
        raise ValueError("None of the selected residues were found in the data files.")

    plot_df = mat.loc[residues].copy()

    ligand_order = [x for x in CUSTOM_LIGAND_ORDER if x in plot_df.columns]
    remaining = [x for x in plot_df.columns if x not in ligand_order]
    ligand_order = ligand_order + remaining

    if len(ligand_order) == 0:
        raise ValueError("No ligand columns available after filename parsing.")

    plot_df = plot_df[ligand_order]

    print("\nLigand order used in plot:")
    print(ligand_order)

    # Transpose:
    # rows = ligands
    # columns = residues
    plot_df = plot_df.T

    vals = plot_df.values.astype(float)

    actual_min = np.nanmin(vals)
    actual_max = np.nanmax(vals)

    vmin = round_down(actual_min, ROUND_DECIMALS)
    vmax = round_up(actual_max, ROUND_DECIMALS)

    print(f"\nActual data minimum: {actual_min:.3f}")
    print(f"Actual data maximum: {actual_max:.3f}")
    print(f"Colorbar limits used: {vmin} to {vmax}")

    if not (vmin < 0 < vmax):
        raise ValueError(
            "TwoSlopeNorm requires vmin < 0 < vmax. "
            "Your selected data may not contain both negative and positive values."
        )

    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    # Dedicated colorbar axis fixes blank / missing colorbar issue
    fig = plt.figure(figsize=FIGSIZE)
    gs = fig.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=[40, 1],
        wspace=0.05
    )

    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    im = ax.imshow(
        vals,
        cmap=CMAP,
        norm=norm,
        aspect="auto",
        interpolation="nearest"
    )

    x_labels = get_residue_labels(plot_df.columns, residue_map)
    y_labels = list(plot_df.index)

    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_xticklabels(
        x_labels,
        rotation=45,
        ha="right",
        fontweight="bold"
    )

    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_yticklabels(
        y_labels,
        fontweight="bold"
    )

    ax.set_xlabel("Residue", fontweight="bold")
    ax.set_ylabel("Ligands", fontweight="bold")
    ax.set_title(
        "Residue-wise binding free-energy decomposition",
        fontweight="bold",
        pad=16
    )

    # Grid between heatmap cells
    ax.set_xticks(np.arange(-0.5, len(x_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(y_labels), 1), minor=True)

    ax.grid(
        which="minor",
        color="white",
        linestyle="-",
        linewidth=1.0,
        alpha=0.85
    )

    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="both", direction="out", length=5, width=1.1)

    # Add numerical values
    if SHOW_CELL_VALUES:
        for i in range(vals.shape[0]):
            for j in range(vals.shape[1]):
                value = vals[i, j]
                ax.text(
                    j,
                    i,
                    VALUE_FMT.format(value),
                    ha="center",
                    va="center",
                    fontsize=FS_CELL,
                    color=choose_text_color(value, vmin, vmax),
                    fontweight="bold"
                )

    # Colorbar
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(
        "Per-residue contribution (kcal/mol)",
        fontsize=FS_CBAR,
        fontweight="bold",
        rotation=90,
        labelpad=18
    )

    cbar_ticks = make_colorbar_ticks(vmin, vmax)
    cbar.set_ticks(cbar_ticks)
    cbar.ax.tick_params(
        labelsize=12,
        width=1.1,
        length=5
    )

    for label in cbar.ax.get_yticklabels():
        label.set_fontweight("bold")

    # This helps avoid white gaps in PDF/SVG colorbars
    cbar.solids.set_edgecolor("face")

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    for spine in cax.spines.values():
        spine.set_linewidth(1.2)

    fig.subplots_adjust(
        left=0.12,
        right=0.94,
        bottom=0.20,
        top=0.88
    )

    plt.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
    plt.savefig(OUT_PDF, bbox_inches="tight")
    plt.savefig(OUT_SVG, bbox_inches="tight")

    print(f"\nSaved: {OUT_PNG}")
    print(f"Saved: {OUT_PDF}")
    print(f"Saved: {OUT_SVG}")

    plt.show()


if __name__ == "__main__":
    main()
