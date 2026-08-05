#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.interpolate import griddata
from scipy.ndimage import minimum_filter

FEL_FILE = "Gibbs.txt"
PROJECTION_FILE = "pc12_time.dat"

SYSTEM_NAME = "Control-Protein"
BASIN_PREFIX = "A"  
NUMBER_OF_MINIMA = 2
ENERGY_WINDOW = 3.0
MINIMUM_NORMALIZED_SEPARATION = 0.18
MINIMUM_FILTER_SIZE = 9
NGRID = 300
PLOT_MAX_ENERGY = 12.0
ENERGY_UNIT = "kcal/mol"
FONT_FILE = "timesbold.ttf"

OUTPUT_FIGURE = "FEL_with_representative_frames.png"
OUTPUT_CSV = "representative_frames.csv"
OUTPUT_SHELL_SCRIPT = "extract_representative_frames.sh"


TPR_FILE = "md.tpr"
TRAJECTORY_FILE = "md_fitted.xtc"
INDEX_FILE = "index.ndx"

TIME_UNIT = "ps"

OUTPUT_GROUP_NUMBER: int | None = None

def load_numeric_file(
    filename: str,
    minimum_columns: int,
    description: str,
) -> np.ndarray:
    """Load a whitespace-separated numeric file."""

    path = Path(filename)

    if not path.exists():
        raise FileNotFoundError(
            f"{description} file was not found: {filename}"
        )

    try:
        data = np.loadtxt(path, dtype=float)
    except ValueError as exc:
        raise ValueError(
            f"Could not read numeric data from {filename}. "
            "Check for non-numeric lines or malformed columns."
        ) from exc

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < minimum_columns:
        raise ValueError(
            f"{filename} must contain at least {minimum_columns} columns, "
            f"but {data.shape[1]} column(s) were detected."
        )

    if data.shape[0] == 0:
        raise ValueError(f"No data rows were found in {filename}.")

    return data


def check_filter_size(value: int) -> int:
    """Ensure that the minimum-filter size is a positive odd integer."""

    if value < 1:
        raise ValueError("MINIMUM_FILTER_SIZE must be at least 1.")

    if value % 2 == 0:
        value += 1
        print(
            f"MINIMUM_FILTER_SIZE changed to {value} "
            "because an odd value is required."
        )

    return value


fel_data = load_numeric_file(
    FEL_FILE,
    minimum_columns=3,
    description="FEL",
)

fel_pc1 = fel_data[:, 0]
fel_pc2 = fel_data[:, 1]
fel_energy = fel_data[:, 2]

valid_fel = (
    np.isfinite(fel_pc1)
    & np.isfinite(fel_pc2)
    & np.isfinite(fel_energy)
)

fel_pc1 = fel_pc1[valid_fel]
fel_pc2 = fel_pc2[valid_fel]
fel_energy = fel_energy[valid_fel]

if len(fel_energy) == 0:
    raise ValueError("No finite FEL data points remain after filtering.")

fel_energy = fel_energy - np.min(fel_energy)

pc1_min = float(np.min(fel_pc1))
pc1_max = float(np.max(fel_pc1))
pc2_min = float(np.min(fel_pc2))
pc2_max = float(np.max(fel_pc2))

pc1_range = pc1_max - pc1_min
pc2_range = pc2_max - pc2_min

if pc1_range <= 0 or pc2_range <= 0:
    raise ValueError(
        "The FEL must span nonzero ranges along both PC1 and PC2."
    )



projection_data = load_numeric_file(
    PROJECTION_FILE,
    minimum_columns=3,
    description="PC projection",
)

trajectory_time = projection_data[:, 0]
trajectory_pc1 = projection_data[:, 1]
trajectory_pc2 = projection_data[:, 2]

valid_projection = (
    np.isfinite(trajectory_time)
    & np.isfinite(trajectory_pc1)
    & np.isfinite(trajectory_pc2)
)

trajectory_time = trajectory_time[valid_projection]
trajectory_pc1 = trajectory_pc1[valid_projection]
trajectory_pc2 = trajectory_pc2[valid_projection]

if len(trajectory_time) == 0:
    raise ValueError(
        "No valid time-resolved projection data were found."
    )

if np.any(np.diff(trajectory_time) < 0):
    raise ValueError(
        "The times in pc12_time.dat are not in increasing order."
    )

print(f"FEL points loaded: {len(fel_energy)}")
print(f"Trajectory frames loaded: {len(trajectory_time)}")

print(
    f"FEL PC1 range: {pc1_min:.4f} to {pc1_max:.4f}"
)
print(
    f"Trajectory PC1 range: "
    f"{np.min(trajectory_pc1):.4f} to "
    f"{np.max(trajectory_pc1):.4f}"
)

print(
    f"FEL PC2 range: {pc2_min:.4f} to {pc2_max:.4f}"
)
print(
    f"Trajectory PC2 range: "
    f"{np.min(trajectory_pc2):.4f} to "
    f"{np.max(trajectory_pc2):.4f}"
)


grid_pc1 = np.linspace(pc1_min, pc1_max, NGRID)
grid_pc2 = np.linspace(pc2_min, pc2_max, NGRID)

GRID_PC1, GRID_PC2 = np.meshgrid(grid_pc1, grid_pc2)

GRID_ENERGY_LINEAR = griddata(
    points=(fel_pc1, fel_pc2),
    values=fel_energy,
    xi=(GRID_PC1, GRID_PC2),
    method="linear",
)
GRID_ENERGY_NEAREST = griddata(
    points=(fel_pc1, fel_pc2),
    values=fel_energy,
    xi=(GRID_PC1, GRID_PC2),
    method="nearest",
)

GRID_ENERGY_PLOT = GRID_ENERGY_LINEAR.copy()

missing_plot_values = ~np.isfinite(GRID_ENERGY_PLOT)
GRID_ENERGY_PLOT[missing_plot_values] = (
    GRID_ENERGY_NEAREST[missing_plot_values]
)


sampled_grid_mask = np.isfinite(GRID_ENERGY_LINEAR)

working_energy = np.where(
    sampled_grid_mask,
    GRID_ENERGY_LINEAR,
    np.inf,
)

neighborhood_minimum = minimum_filter(
    working_energy,
    size=filter_size,
    mode="constant",
    cval=np.inf,
)

local_minimum_mask = (
    sampled_grid_mask
    & np.isclose(
        working_energy,
        neighborhood_minimum,
        rtol=1.0e-7,
        atol=1.0e-8,
    )
    & (working_energy <= ENERGY_WINDOW)
)

candidate_grid_indices = np.argwhere(local_minimum_mask)

candidate_minima: list[dict[str, float]] = []

for row_index, column_index in candidate_grid_indices:
    candidate_minima.append(
        {
            "pc1": float(GRID_PC1[row_index, column_index]),
            "pc2": float(GRID_PC2[row_index, column_index]),
            "energy": float(
                GRID_ENERGY_LINEAR[row_index, column_index]
            ),
        }
    )

candidate_minima.sort(key=lambda item: item["energy"])

selected_minima: list[dict[str, float]] = []

for candidate in candidate_minima:
    keep_candidate = True

    for selected in selected_minima:
        normalized_distance = math.sqrt(
            (
                (candidate["pc1"] - selected["pc1"])
                / pc1_range
            ) ** 2
            + (
                (candidate["pc2"] - selected["pc2"])
                / pc2_range
            ) ** 2
        )

        if (
            normalized_distance
            < MINIMUM_NORMALIZED_SEPARATION
        ):
            keep_candidate = False
            break

    if keep_candidate:
        selected_minima.append(candidate)

    if len(selected_minima) == NUMBER_OF_MINIMA:
        break

if len(selected_minima) < NUMBER_OF_MINIMA:
    print(
        "\nWarning: the local-minimum detector identified only "
        f"{len(selected_minima)} distinct minimum/minima."
    )
    print(
        "Applying a fallback search using the original FEL points."
    )

    sorted_original_indices = np.argsort(fel_energy)

    for original_index in sorted_original_indices:
        candidate_energy = float(fel_energy[original_index])

        if candidate_energy > ENERGY_WINDOW:
            continue

        candidate = {
            "pc1": float(fel_pc1[original_index]),
            "pc2": float(fel_pc2[original_index]),
            "energy": candidate_energy,
        }

        keep_candidate = True

        for selected in selected_minima:
            normalized_distance = math.sqrt(
                (
                    (candidate["pc1"] - selected["pc1"])
                    / pc1_range
                ) ** 2
                + (
                    (candidate["pc2"] - selected["pc2"])
                    / pc2_range
                ) ** 2
            )

            if (
                normalized_distance
                < MINIMUM_NORMALIZED_SEPARATION
            ):
                keep_candidate = False
                break

        if keep_candidate:
            selected_minima.append(candidate)

        if len(selected_minima) == NUMBER_OF_MINIMA:
            break

if not selected_minima:
    raise RuntimeError(
        "No minima were identified. Increase ENERGY_WINDOW or "
        "decrease MINIMUM_NORMALIZED_SEPARATION."
    )

if len(selected_minima) < NUMBER_OF_MINIMA:
    print(
        f"\nOnly {len(selected_minima)} distinct low-energy "
        "basin(s) were identified."
    )

trajectory_pc1_scale = max(
    float(np.ptp(trajectory_pc1)),
    1.0e-12,
)

trajectory_pc2_scale = max(
    float(np.ptp(trajectory_pc2)),
    1.0e-12,
)

representative_frames: list[dict[str, float | int | str]] = []

for basin_number, minimum in enumerate(
    selected_minima,
    start=1,
):
    basin_label = f"{BASIN_PREFIX}{basin_number}"

    scaled_squared_distance = (
        (
            (trajectory_pc1 - minimum["pc1"])
            / trajectory_pc1_scale
        ) ** 2
        + (
            (trajectory_pc2 - minimum["pc2"])
            / trajectory_pc2_scale
        ) ** 2
    )

    nearest_index = int(np.argmin(scaled_squared_distance))

    representative_frames.append(
        {
            "basin": basin_label,
            "minimum_PC1": minimum["pc1"],
            "minimum_PC2": minimum["pc2"],
            f"relative_energy_{ENERGY_UNIT.replace('/', '_per_')}":
                minimum["energy"],
            "trajectory_frame_index_zero_based":
                nearest_index,
            f"trajectory_time_{TIME_UNIT}":
                float(trajectory_time[nearest_index]),
            "frame_PC1":
                float(trajectory_pc1[nearest_index]),
            "frame_PC2":
                float(trajectory_pc2[nearest_index]),
            "normalized_PC_distance":
                float(
                    math.sqrt(
                        scaled_squared_distance[nearest_index]
                    )
                ),
        }
    )


csv_fieldnames = list(representative_frames[0].keys())

with open(
    OUTPUT_CSV,
    "w",
    newline="",
    encoding="utf-8",
) as csv_handle:
    writer = csv.DictWriter(
        csv_handle,
        fieldnames=csv_fieldnames,
    )
    writer.writeheader()
    writer.writerows(representative_frames)


print("\nRepresentative FEL structures")
print("=" * 90)

for frame in representative_frames:
    energy_key = (
        f"relative_energy_{ENERGY_UNIT.replace('/', '_per_')}"
    )
    time_key = f"trajectory_time_{TIME_UNIT}"

    print(
        f'{frame["basin"]}: '
        f'minimum = '
        f'({frame["minimum_PC1"]:.5f}, '
        f'{frame["minimum_PC2"]:.5f}); '
        f'ΔG = {frame[energy_key]:.4f} {ENERGY_UNIT}; '
        f'nearest actual frame = '
        f'{frame[time_key]:.4f} {TIME_UNIT}; '
        f'frame PC coordinates = '
        f'({frame["frame_PC1"]:.5f}, '
        f'{frame["frame_PC2"]:.5f})'
    )

print(f"\nFrame information saved to: {OUTPUT_CSV}")

font_path = Path(FONT_FILE)

if font_path.exists():
    font_manager.fontManager.addfont(str(font_path))
    custom_font_name = font_manager.FontProperties(
        fname=str(font_path)
    ).get_name()
    plt.rcParams["font.family"] = custom_font_name
else:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = [
        "Times New Roman",
        "DejaVu Serif",
    ]
    print(
        f"\nFont file {FONT_FILE} was not found. "
        "Using the available serif font."
    )

plt.rcParams["font.size"] = 15
plt.rcParams["axes.labelsize"] = 16
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["axes.linewidth"] = 1.5

figure, axis = plt.subplots(figsize=(8.4, 7.0))

if PLOT_MAX_ENERGY is None:
    finite_plot_energies = GRID_ENERGY_PLOT[
        np.isfinite(GRID_ENERGY_PLOT)
    ]
    plot_maximum = float(
        np.ceil(np.percentile(finite_plot_energies, 98))
    )
else:
    plot_maximum = float(PLOT_MAX_ENERGY)

if plot_maximum <= 0:
    raise ValueError("PLOT_MAX_ENERGY must be greater than zero.")

contour_levels = np.linspace(
    0.0,
    plot_maximum,
    13,
)

filled_contour = axis.contourf(
    GRID_PC1,
    GRID_PC2,
    GRID_ENERGY_PLOT,
    levels=contour_levels,
    cmap="viridis",
    extend="max",
)

axis.contour(
    GRID_PC1,
    GRID_PC2,
    GRID_ENERGY_PLOT,
    levels=contour_levels,
    linewidths=0.35,
    alpha=0.35,
)

for frame in representative_frames:
    minimum_pc1 = float(frame["minimum_PC1"])
    minimum_pc2 = float(frame["minimum_PC2"])
    actual_pc1 = float(frame["frame_PC1"])
    actual_pc2 = float(frame["frame_PC2"])
    basin_label = str(frame["basin"])

    # FEL minimum
    axis.scatter(
        minimum_pc1,
        minimum_pc2,
        marker="*",
        s=300,
        edgecolors="black",
        linewidths=1.2,
        zorder=10,
        label=(
            "FEL minimum"
            if basin_label == representative_frames[0]["basin"]
            else None
        ),
    )

    # Nearest actual trajectory frame
    axis.scatter(
        actual_pc1,
        actual_pc2,
        marker="o",
        s=90,
        facecolors="none",
        edgecolors="black",
        linewidths=1.4,
        zorder=11,
        label=(
            "Selected MD frame"
            if basin_label == representative_frames[0]["basin"]
            else None
        ),
    )

    # Connect minimum and selected frame
    axis.plot(
        [minimum_pc1, actual_pc1],
        [minimum_pc2, actual_pc2],
        linestyle="--",
        linewidth=0.9,
        alpha=0.8,
        zorder=9,
    )

    axis.annotate(
        basin_label,
        xy=(minimum_pc1, minimum_pc2),
        xytext=(7, 8),
        textcoords="offset points",
        fontsize=15,
        fontweight="bold",
        zorder=12,
    )

axis.set_xlabel("PC1")
axis.set_ylabel("PC2")
axis.set_title(
    f"Free-energy landscape: {SYSTEM_NAME}",
    pad=12,
)

colorbar = figure.colorbar(
    filled_contour,
    ax=axis,
    fraction=0.050,
    pad=0.04,
)

colorbar.set_label(
    f"Relative free energy ({ENERGY_UNIT})"
)

axis.legend(
    loc="best",
    frameon=True,
    fontsize=11,
)

figure.tight_layout()

figure.savefig(
    OUTPUT_FIGURE,
    dpi=600,
    bbox_inches="tight",
)

plt.show()

print(f"Annotated FEL saved to: {OUTPUT_FIGURE}")

with open(
    OUTPUT_SHELL_SCRIPT,
    "w",
    encoding="utf-8",
) as shell_handle:
    shell_handle.write("#!/usr/bin/env bash\n")
    shell_handle.write("set -euo pipefail\n\n")

    shell_handle.write(
        "# Representative structures selected from FEL minima\n"
    )
    shell_handle.write(
        "# Confirm that the trajectory contains CETP and the ligand.\n"
    )
    shell_handle.write(
        "# The selected output group should contain CETP + ligand.\n\n"
    )

    for frame in representative_frames:
        basin_label = str(frame["basin"])
        time_value = float(
            frame[f"trajectory_time_{TIME_UNIT}"]
        )

        output_pdb = (
            f"{SYSTEM_NAME.replace('–', '_').replace('-', '_')}"
            f"_{basin_label}_representative.pdb"
        )

        command = (
            f"gmx trjconv "
            f"-s {TPR_FILE} "
            f"-f {TRAJECTORY_FILE} "
            f"-n {INDEX_FILE} "
            f"-tu {TIME_UNIT} "
            f"-dump {time_value:.6f} "
            f"-o {output_pdb}"
        )

        shell_handle.write(
            f'echo "Extracting {basin_label} at '
            f'{time_value:.6f} {TIME_UNIT}"\n'
        )

        if OUTPUT_GROUP_NUMBER is None:
            shell_handle.write(
                "# Select the CETP + ligand group when prompted.\n"
            )
            shell_handle.write(f"{command}\n\n")
        else:
            shell_handle.write(
                f'echo "{OUTPUT_GROUP_NUMBER}" | {command}\n\n'
            )

os.chmod(OUTPUT_SHELL_SCRIPT, 0o755)

print(
    "GROMACS extraction commands saved to: "
    f"{OUTPUT_SHELL_SCRIPT}"
)
