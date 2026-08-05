from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
from matplotlib import font_manager
from matplotlib.ticker import MaxNLocator

SYSTEMS = [
    {
        "name": "Ligand 1–Protein",
        "folder": Path("ligand1"),
        "output": "FEL_Ligand1.png",
    },
    {
        "name": "Ligand 2–Protein",
        "folder": Path("ligand2"),
        "output": "FEL_Ligand2.png",
    },
]

FONT_FILE = Path("timesbold.ttf")

NGRID = 200
FEL_LEVELS = np.arange(0, 13, 1)
NUMBER_OF_PC_TICKS = 6
PC_RANGE_PADDING = 0.00
REMAP_XVG_BIN_AXIS = True

if not FONT_FILE.exists():
    raise FileNotFoundError(
        f"Font file was not found: {FONT_FILE.resolve()}"
    )

font_manager.fontManager.addfont(str(FONT_FILE))

timesbold_name = font_manager.FontProperties(
    fname=str(FONT_FILE)
).get_name()

plt.rcParams["font.family"] = timesbold_name
plt.rcParams["font.size"] = 20
plt.rcParams["axes.titlesize"] = 22
plt.rcParams["axes.labelsize"] = 20
plt.rcParams["xtick.labelsize"] = 20
plt.rcParams["ytick.labelsize"] = 20
plt.rcParams["axes.linewidth"] = 2.0


def load_xvg(filename):
    filename = Path(filename)

    if not filename.exists():
        raise FileNotFoundError(
            f"XVG file was not found: {filename.resolve()}"
        )

    data = []

    with filename.open("r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith("@") or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])

    if not data:
        raise ValueError(
            f"No numerical data were found in {filename.resolve()}"
        )

    return np.asarray(data, dtype=float)


def load_gibbs(filename):
    filename = Path(filename)

    if not filename.exists():
        raise FileNotFoundError(
            f"Gibbs file was not found: {filename.resolve()}"
        )

    data = np.loadtxt(filename)

    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(
            f"{filename.resolve()} must contain at least three columns: "
            "PC1, PC2 and free energy."
        )

    pc1 = data[:, 0]
    pc2 = data[:, 1]
    energy = data[:, 2]

    valid = np.isfinite(pc1) & np.isfinite(pc2) & np.isfinite(energy)

    pc1 = pc1[valid]
    pc2 = pc2[valid]
    energy = energy[valid]

    if pc1.size == 0:
        raise ValueError(
            f"No valid numerical data were found in {filename.resolve()}"
        )

    return pc1, pc2, energy


def calculate_common_range(arrays, padding_fraction=0.0):
    combined = np.concatenate(
        [np.asarray(array, dtype=float).ravel() for array in arrays]
    )

    combined = combined[np.isfinite(combined)]

    if combined.size == 0:
        raise ValueError("No finite values were available for range calculation.")

    value_min = np.min(combined)
    value_max = np.max(combined)

    value_width = value_max - value_min

    if value_width == 0:
        padding = 0.5
    else:
        padding = value_width * padding_fraction

    return value_min - padding, value_max + padding


def map_profile_axis_to_fel(profile_x, fel_coordinate, remap_bin_axis=True):
    profile_x = np.asarray(profile_x, dtype=float)
    fel_coordinate = np.asarray(fel_coordinate, dtype=float)

    if not remap_bin_axis:
        return profile_x.copy()

    if not np.all(np.isfinite(profile_x)):
        raise ValueError("The XVG profile x-axis contains non-finite values.")

    fel_coordinate = fel_coordinate[np.isfinite(fel_coordinate)]

    if fel_coordinate.size == 0:
        raise ValueError("The FEL coordinate array contains no finite values.")

    unique_fel_coordinates = np.unique(fel_coordinate)
    unique_fel_coordinates.sort()

    if len(profile_x) == len(unique_fel_coordinates):
        return unique_fel_coordinates

    source_min = np.min(profile_x)
    source_max = np.max(profile_x)

    target_min = np.min(fel_coordinate)
    target_max = np.max(fel_coordinate)

    if source_max == source_min:
        return np.full_like(
            profile_x,
            fill_value=(target_min + target_max) / 2.0,
            dtype=float,
        )

    normalized_position = (profile_x - source_min) / (source_max - source_min)

    mapped_axis = target_min + normalized_position * (target_max - target_min)

    return mapped_axis


def create_common_ticks(axis_limits, number_of_ticks=6):
    locator = MaxNLocator(nbins=number_of_ticks)

    ticks = locator.tick_values(axis_limits[0], axis_limits[1])

    tolerance = abs(axis_limits[1] - axis_limits[0]) * 1.0e-9

    ticks = ticks[
        (ticks >= axis_limits[0] - tolerance)
        & (ticks <= axis_limits[1] + tolerance)
    ]

    return ticks


for system in SYSTEMS:
    folder = system["folder"]

    system["gibbs_file"] = folder / "Gibbs.txt"
    system["pc1_file"] = folder / "energy_PC1.xvg"
    system["pc2_file"] = folder / "energy_PC2.xvg"

    system["fel_pc1"], system["fel_pc2"], system["fel_energy"] = load_gibbs(
        system["gibbs_file"]
    )

    system["pc1_profile"] = load_xvg(system["pc1_file"])
    system["pc2_profile"] = load_xvg(system["pc2_file"])

    pc1_profile_bin = system["pc1_profile"][:, 0]
    pc2_profile_bin = system["pc2_profile"][:, 0]

    system["pc1_profile_coordinate"] = map_profile_axis_to_fel(
        profile_x=pc1_profile_bin,
        fel_coordinate=system["fel_pc1"],
        remap_bin_axis=REMAP_XVG_BIN_AXIS,
    )

    system["pc2_profile_coordinate"] = map_profile_axis_to_fel(
        profile_x=pc2_profile_bin,
        fel_coordinate=system["fel_pc2"],
        remap_bin_axis=REMAP_XVG_BIN_AXIS,
    )


GLOBAL_PC1_LIMITS = calculate_common_range(
    arrays=[system["fel_pc1"] for system in SYSTEMS],
    padding_fraction=PC_RANGE_PADDING,
)

GLOBAL_PC2_LIMITS = calculate_common_range(
    arrays=[system["fel_pc2"] for system in SYSTEMS],
    padding_fraction=PC_RANGE_PADDING,
)

GLOBAL_PC1_TICKS = create_common_ticks(
    GLOBAL_PC1_LIMITS,
    number_of_ticks=NUMBER_OF_PC_TICKS,
)

GLOBAL_PC2_TICKS = create_common_ticks(
    GLOBAL_PC2_LIMITS,
    number_of_ticks=NUMBER_OF_PC_TICKS,
)

print("\nCommon ranges used for every ligand:")
print(f"PC1: {GLOBAL_PC1_LIMITS[0]:.4f} to {GLOBAL_PC1_LIMITS[1]:.4f}")
print(f"PC2: {GLOBAL_PC2_LIMITS[0]:.4f} to {GLOBAL_PC2_LIMITS[1]:.4f}")


def create_figure(system):
    system_name = system["name"]
    output_file = system["output"]

    x = system["fel_pc1"]
    y = system["fel_pc2"]
    z = system["fel_energy"]

    pc1_x = system["pc1_profile_coordinate"]
    pc1_y = system["pc1_profile"][:, 1]

    pc2_x = system["pc2_profile_coordinate"]
    pc2_y = system["pc2_profile"][:, 1]

   
    xi = np.linspace(np.min(x), np.max(x), NGRID)
    yi = np.linspace(np.min(y), np.max(y), NGRID)
    XI, YI = np.meshgrid(xi, yi)

    ZI = griddata((x, y), z, (XI, YI), method="linear")

    ZI_linear = griddata((x, y), z, (XI, YI), method="linear")
    ZI_nearest = griddata((x, y), z, (XI, YI), method="nearest")

    nan_mask = np.isnan(ZI)
    ZI[nan_mask] = ZI_linear[nan_mask]

    nan_mask = np.isnan(ZI)
    ZI[nan_mask] = ZI_nearest[nan_mask]

   
    fig = plt.figure(figsize=(13, 7))

    gs = GridSpec(
        2, 2,
        width_ratios=[1.0, 1.35],
        height_ratios=[1, 1],
        wspace=0.30,
        hspace=0.30,
    )

    ax_pc1 = fig.add_subplot(gs[0, 0])
    ax_pc2 = fig.add_subplot(gs[1, 0])
    ax_fel = fig.add_subplot(gs[:, 1])
    
    line_color = "#B13A1B"
    cmap = plt.get_cmap("RdYlBu_r").copy()
    cmap.set_bad(cmap(1.0))

    ax_fel.set_facecolor(cmap(1.0))
    
    ax_pc1.plot(pc1_x, pc1_y, color=line_color, linewidth=2.2)
    ax_pc1.set_xlabel("PC1")
    ax_pc1.set_ylabel("ΔG (kJ/mol)")
    ax_pc1.set_xlim(GLOBAL_PC1_LIMITS)
    ax_pc1.set_xticks(GLOBAL_PC1_TICKS)
    ax_pc1.grid(False)
    ax_pc2.plot(pc2_x, pc2_y, color=line_color, linewidth=2.2)
    ax_pc2.set_xlabel("PC2")
    ax_pc2.set_ylabel("ΔG (kJ/mol)")
    ax_pc2.set_xlim(GLOBAL_PC2_LIMITS)
    ax_pc2.set_xticks(GLOBAL_PC2_TICKS)
    ax_pc2.grid(False)

    for axis in [ax_pc1, ax_pc2]:
        for spine in axis.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.2)
            spine.set_color("black")

    cf = ax_fel.contourf(
        XI,
        YI,
        ZI,
        levels=FEL_LEVELS,
        cmap=cmap,
        extend="max",
    )

    ax_fel.contour(
        XI,
        YI,
        ZI,
        levels=FEL_LEVELS,
        colors="k",
        linewidths=0.2,
        alpha=0.25,
    )

    ax_fel.set_xlabel("PC1")
    ax_fel.set_ylabel("PC2")
    ax_fel.set_xlim(GLOBAL_PC1_LIMITS)
    ax_fel.set_xticks(GLOBAL_PC1_TICKS)
    ax_fel.set_ylim(GLOBAL_PC2_LIMITS)
    ax_fel.set_yticks(GLOBAL_PC2_TICKS)

    cbar = fig.colorbar(
        cf,
        ax=ax_fel,
        fraction=0.046,
        pad=0.04,
    )
    cbar.set_label("Free energy (kJ/mol)", fontsize=16)
    
    fig.suptitle(
        system_name,
        fontsize=18,
        fontweight="bold",
        y=0.99,
    )

    plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])

    plt.savefig(
        output_file,
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {output_file}")
    
for system in SYSTEMS:
    create_figure(system)
