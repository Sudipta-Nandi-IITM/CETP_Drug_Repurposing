from __future__ import annotations

import csv
import math
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from pymol.cgo import CONE, CYLINDER
import matplotlib.cm as cm
import os

import numpy as np
from pymol import cmd

import matplotlib
matplotlib.use("Agg")  

import matplotlib.pyplot as plt
from matplotlib import font_manager

FONT_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "timesbold.ttf",
)

if not os.path.isfile(FONT_FILE):
    raise FileNotFoundError(
        f"Font file was not found: {FONT_FILE}"
    )

font_manager.fontManager.addfont(FONT_FILE)

timesbold_prop = font_manager.FontProperties(
    fname=FONT_FILE
)

TIMESBOLD_NAME = timesbold_prop.get_name()

print(f"Colour-bar font: {TIMESBOLD_NAME}")


PALETTE_NAME = "RdYlBu_r"


def register_rdylbu(n=64, name="rdylbu_ramp"):
    """
    Register discrete RdYlBu_r colours in PyMOL.

    Low values  -> blue
    High values -> red
    """
    cmap = cm.get_cmap(PALETTE_NAME, n)

    color_names = []

    for index in range(n):
        r, g, b, _ = cmap(index)

        color_name = f"{name}_{index:03d}"

        cmd.set_color(
            color_name,
            [float(r), float(g), float(b)],
        )

        color_names.append(color_name)

    return " ".join(color_names)


def magnitude_rgb_rdylbu(value, minimum, maximum):
    
    if maximum <= minimum:
        return 0.5, 0.5, 0.5

    fraction = (value - minimum) / (maximum - minimum)
    fraction = max(0.0, min(1.0, fraction))

    cmap = cm.get_cmap(PALETTE_NAME)

    r, g, b, _ = cmap(fraction)

    return float(r), float(g), float(b)
    
def register_viridis(n=64, name="viridis_ramp"):
    
    cmap = cm.get_cmap("viridis", n)
    names = []
    for i in range(n):
        r, g, b, _ = cmap(i)
        color_name = f"{name}_{i:03d}"
        cmd.set_color(color_name, [r, g, b])
        names.append(color_name)
    return " ".join(names)

def magnitude_rgb_viridis(value, minimum, maximum):
    if maximum <= minimum:
        return 0.5, 0.5, 0.5
    fraction = (value - minimum) / (maximum - minimum)
    fraction = max(0.0, min(1.0, fraction))
    cmap = cm.get_cmap("viridis")
    r, g, b, _ = cmap(fraction)
    return r, g, b

def pc1_colored_motion(
    start_object="pc1_neg",
    end_object="pc1_pos",
    display_object="pc1_neg",
    arrow_object="pc1_vectors",
    every=14,
    arrow_scale=0.25,
    cutoff=3.0,
    minimum=0.0,
    maximum=0.0,
    shaft_radius=0.055,
    head_radius=0.14,
    color_arrows=0,
    csv_file="pc1_displacement.csv",
):

    every = max(1, int(every))
    arrow_scale = float(arrow_scale)
    cutoff = float(cutoff)
    minimum = float(minimum)
    maximum = float(maximum)
    shaft_radius = float(shaft_radius)
    head_radius = float(head_radius)
    color_arrows = int(color_arrows)

    start_model = cmd.get_model(
        f"({start_object}) and name CA",
        state=1,
    )

    end_model = cmd.get_model(
        f"({end_object}) and name CA",
        state=1,
    )

    start_atoms = start_model.atom
    end_atoms = end_model.atom

    print(f"Starting C-alpha atoms: {len(start_atoms)}")
    print(f"Ending C-alpha atoms:   {len(end_atoms)}")

    if not start_atoms or not end_atoms:
        raise ValueError(
            "No C-alpha atoms were found in one or both endpoint objects."
        )

    if len(start_atoms) != len(end_atoms):
        raise ValueError(
            "The endpoint structures contain different numbers "
            "of C-alpha atoms."
        )

    displacement_data = []

    for start_atom, end_atom in zip(start_atoms, end_atoms):
        vector = [
            float(end_atom.coord[index] - start_atom.coord[index])
            for index in range(3)
        ]

        magnitude = math.sqrt(
            sum(component * component for component in vector)
        )

        displacement_data.append(
            {
                "chain": start_atom.chain,
                "residue_number": start_atom.resi,
                "residue_name": start_atom.resn,
                "atom_name": start_atom.name,
                "dx_A": vector[0],
                "dy_A": vector[1],
                "dz_A": vector[2],
                "displacement_A": magnitude,
            }
        )

    displacement_values = [
        row["displacement_A"]
        for row in displacement_data
    ]

    mean_value = (
        sum(displacement_values) / len(displacement_values)
    )

    rms_value = math.sqrt(
        sum(value * value for value in displacement_values)
        / len(displacement_values)
    )

    observed_maximum = max(displacement_values)

    if maximum <= minimum:
        maximum = observed_maximum

    print(f"Mean displacement: {mean_value:.4f} Å")
    print(f"RMS displacement:  {rms_value:.4f} Å")
    print(f"Maximum observed:  {observed_maximum:.4f} Å")
    print(
        f"Color scale:       {minimum:.4f}–{maximum:.4f} Å"
    )

    with open(csv_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(displacement_data[0].keys()),
        )
        writer.writeheader()
        writer.writerows(displacement_data)

    print(f"Displacement table written to: {csv_file}")

   
    for row in displacement_data:
        chain = row["chain"]
        residue_number = row["residue_number"]
        value = row["displacement_A"]

        selection = (
            f"({display_object}) and resi {residue_number}"
        )

        if chain:
            selection += f" and chain {chain}"

        cmd.alter(selection, f"b={value}")

    cmd.rebuild(display_object)

   
    palette_str = register_rdylbu(64)
    cmd.spectrum(
    	expression="b",
    	palette=palette_str,
    	selection=display_object,
    	minimum=minimum,
    	maximum=maximum,
    	byres=1,
    )
   
    arrow_geometry = []
    arrows_drawn = 0

    uniform_rgb = (0.15, 0.15, 0.15)

    for atom_index in range(
        0,
        len(start_atoms),
        every,
    ):
        start_atom = start_atoms[atom_index]
        end_atom = end_atoms[atom_index]

        start = [
            float(value)
            for value in start_atom.coord
        ]

        displacement = [
            float(end_atom.coord[index] - start_atom.coord[index])
            for index in range(3)
        ]

        magnitude = math.sqrt(
            sum(value * value for value in displacement)
        )

        if magnitude < cutoff:
            continue

        arrow_tip = [
            start[index] + arrow_scale * displacement[index]
            for index in range(3)
        ]

        displayed_vector = [
            arrow_tip[index] - start[index]
            for index in range(3)
        ]

        displayed_length = math.sqrt(
            sum(value * value for value in displayed_vector)
        )

        if displayed_length < 1.0e-8:
            continue

        unit_vector = [
            value / displayed_length
            for value in displayed_vector
        ]

        head_length = min(
            max(0.40, 2.5 * head_radius),
            0.35 * displayed_length,
        )

        shaft_end = [
            arrow_tip[index]
            - head_length * unit_vector[index]
            for index in range(3)
        ]

        if color_arrows:
            rgb = magnitude_rgb_rdylbu(
                magnitude,
                minimum,
                maximum,
            )
        else:
            rgb = uniform_rgb

        arrow_geometry.extend(
            [
                CYLINDER,
                start[0], start[1], start[2],
                shaft_end[0], shaft_end[1], shaft_end[2],
                shaft_radius,
                rgb[0], rgb[1], rgb[2],
                rgb[0], rgb[1], rgb[2],

                CONE,
                shaft_end[0], shaft_end[1], shaft_end[2],
                arrow_tip[0], arrow_tip[1], arrow_tip[2],
                head_radius,
                0.0,
                rgb[0], rgb[1], rgb[2],
                rgb[0], rgb[1], rgb[2],
                1.0,
                0.0,
            ]
        )

        arrows_drawn += 1

    cmd.delete(arrow_object)

    if arrow_geometry:
        cmd.load_cgo(arrow_geometry, arrow_object)

    print(
        f"Generated {arrows_drawn} arrows in "
        f"'{arrow_object}'."
    )


cmd.extend("pc1_colored_motion", pc1_colored_motion)

def render_pc1_with_colorbar(
    raw_png="pc1_motion_raw.png",
    output_png="pc1_motion_with_colorbar.png",
    minimum=0.0,
    maximum=15.0,
    width=2400,
    height=1200,
    dpi=600,
):

    minimum = float(minimum)
    maximum = float(maximum)
    width = int(width)
    height = int(height)
    dpi = int(dpi)

    if maximum <= minimum:
        raise ValueError(
            "maximum must be greater than minimum"
        )

    raw_png = os.path.abspath(raw_png)
    output_png = os.path.abspath(output_png)

    print(f"Raw PyMOL image: {raw_png}")
    print(f"Final image:      {output_png}")

    # --------------------------------------------------------
    # Render current PyMOL view
    # --------------------------------------------------------
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)

    cmd.ray(width, height)

    cmd.png(
        raw_png,
        width=width,
        height=height,
        dpi=300,
        ray=0,
        quiet=0,
    )

    cmd.sync()

    if not os.path.isfile(raw_png):
        raise FileNotFoundError(
            f"PyMOL did not create the raw image: {raw_png}"
        )

    structure_image = plt.imread(raw_png)
    fig = plt.figure(
        figsize=(15, 7),
        facecolor="white",
    )
    image_axis = fig.add_axes(
        [0.01, 0.02, 0.84, 0.96]
    )

    image_axis.imshow(structure_image)
    image_axis.axis("off")
    
    colorbar_axis = fig.add_axes(
        [0.875, 0.16, 0.025, 0.68]
    )

    cmap = plt.get_cmap(PALETTE_NAME)

    norm = Normalize(
        vmin=minimum,
        vmax=maximum,
    )

    colorbar = ColorbarBase(
        colorbar_axis,
        cmap=cmap,
        norm=norm,
        orientation="vertical",
        ticks=np.arange(
            minimum,
            maximum + 0.001,
            3.0,
        ),
    )

    colorbar.set_label(
        r"C$\alpha$ displacement ($\AA$)",
        fontsize=18,
        fontproperties=timesbold_prop,
        labelpad=14,
    )

    colorbar.ax.tick_params(
        axis="y",
        labelsize=16,
        width=1.5,
        length=5,
        direction="out",
    )

    for tick_label in colorbar.ax.get_yticklabels():
    	tick_label.set_fontproperties(timesbold_prop)
    	tick_label.set_fontsize(16)

    colorbar.outline.set_linewidth(1.5)

    fig.savefig(
        output_png,
        dpi=dpi,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.05,
    )

    plt.close(fig)

    if not os.path.isfile(output_png):
        raise FileNotFoundError(
            f"The combined image was not created: {output_png}"
        )

    print("Colour-bar image created successfully:")
    print(output_png)


cmd.extend(
    "render_pc1_with_colorbar",
    render_pc1_with_colorbar,
)

