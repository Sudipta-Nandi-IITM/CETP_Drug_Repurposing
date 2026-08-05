from __future__ import annotations

import csv
import math
from pymol import cmd, stored
from pymol.cgo import CONE, CYLINDER
import matplotlib.cm as cm

TRAJECTORY_FILE   = "pc1_extreme_20.pdb"
STABLE_SELECTION  = "resi 1-200 and name CA"   
ARROW_EVERY       = 6
ARROW_CUTOFF      = 1.5
ARROW_SCALE       = 1.0
SHAFT_RADIUS      = 0.06
HEAD_RADIUS       = 0.16
COLOR_ARROWS      = 1
OUTPUT_CSV        = "pc1_displacement.csv"


def register_viridis(n=64, name="viridis_ramp"):
    cmap = cm.get_cmap("viridis", n)
    names = []
    for i in range(n):
        r, g, b, _ = cmap(i)
        color_name = f"{name}_{i:03d}"
        cmd.set_color(color_name, [r, g, b])
        names.append(color_name)
    return names


def magnitude_rgb_viridis(value, minimum, maximum):
    if maximum <= minimum:
        return 0.5, 0.5, 0.5
    fraction = (value - minimum) / (maximum - minimum)
    fraction = max(0.0, min(1.0, fraction))
    cmap = cm.get_cmap("viridis")
    r, g, b, _ = cmap(fraction)
    return r, g, b
    
def fix_model_numbering(input_file, output_file):
    with open(input_file) as handle:
        lines = handle.readlines()

    model_count = 0
    fixed_lines = []
    for line in lines:
        if line.startswith("MODEL"):
            model_count += 1
            fixed_lines.append(f"MODEL     {model_count:4d}\n")
        else:
            fixed_lines.append(line)

    with open(output_file, "w") as handle:
        handle.writelines(fixed_lines)

    print(f"[fix_model_numbering] wrote {model_count} models -> {output_file}")
    return model_count

def load_and_prepare(trajectory_file, stable_selection):
    fixed_file = trajectory_file.replace(".pdb", "_fixed.pdb")
    n_models = fix_model_numbering(trajectory_file, fixed_file)

    cmd.delete("pc1_traj")
    cmd.load(fixed_file, "pc1_traj")
    n_states = cmd.count_states("pc1_traj")
    print(f"[load_and_prepare] loaded states: {n_states} (expected {n_models})")

    if n_states < 2:
        raise RuntimeError("Only 1 state loaded — falling back to manual split.")

    cmd.intra_fit(f"pc1_traj and {stable_selection}")
    cmd.smooth("pc1_traj", window=3, first=0, last=0, ends=1)

    cmd.split_states("pc1_traj", prefix="frame")
    cmd.delete("pc1_traj")

    frame_objects = sorted(name for name in cmd.get_names() if name.startswith("frame"))
    print(f"[load_and_prepare] split into {len(frame_objects)} objects")
    return frame_objects


def load_and_prepare_fallback(trajectory_file, stable_selection):
    with open(trajectory_file) as handle:
        content = handle.read()

    blocks = content.split("MODEL")[1:]
    frame_objects = []
    for i, block in enumerate(blocks):
        fname = f"pc1_frame_{i:02d}.pdb"
        with open(fname, "w") as out:
            out.write("MODEL" + block)
        obj_name = f"frame{i:04d}"
        cmd.load(fname, obj_name)
        frame_objects.append(obj_name)

    ref = frame_objects[0]
    for obj in frame_objects[1:]:
        cmd.align(f"{obj} and {stable_selection}", f"{ref} and {stable_selection}")

    print(f"[load_and_prepare_fallback] loaded {len(frame_objects)} objects")
    return frame_objects

def compute_displacement(frame_objects, csv_file):
    first_obj = frame_objects[0]
    last_obj = frame_objects[-1]

    first_atoms = cmd.get_model(f"{first_obj} and name CA", state=1).atom
    last_atoms = cmd.get_model(f"{last_obj} and name CA", state=1).atom

    if len(first_atoms) != len(last_atoms):
        raise ValueError("First and last frame have different CA counts.")

    displacement_data = []
    for start_atom, end_atom in zip(first_atoms, last_atoms):
        vector = [float(end_atom.coord[i] - start_atom.coord[i]) for i in range(3)]
        magnitude = math.sqrt(sum(v * v for v in vector))
        displacement_data.append(
            {
                "chain": start_atom.chain,
                "residue_number": start_atom.resi,
                "residue_name": start_atom.resn,
                "dx_A": vector[0],
                "dy_A": vector[1],
                "dz_A": vector[2],
                "displacement_A": magnitude,
            }
        )

    with open(csv_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(displacement_data[0].keys()))
        writer.writeheader()
        writer.writerows(displacement_data)

    values = [row["displacement_A"] for row in displacement_data]
    print(f"[compute_displacement] mean={sum(values)/len(values):.4f} A, "
          f"max={max(values):.4f} A -> {csv_file}")

    return displacement_data

def color_trail_viridis(frame_objects):
    n = len(frame_objects)
    color_names = register_viridis(n, name="trail")

    cmd.set("cartoon_trace_atoms", 1)
    cmd.set("cartoon_flat_sheets", 0)
    cmd.set("cartoon_loop_radius", 0.2)
    cmd.set("cartoon_tube_radius", 0.2)

    for i, obj in enumerate(frame_objects):
        cmd.hide("everything", obj)
        cmd.cartoon("tube", obj)
        cmd.show("cartoon", obj)
        cmd.color(color_names[i], obj)

        if i == 0 or i == n - 1:
            cmd.set("cartoon_transparency", 0.0, obj)
        else:
            cmd.set("cartoon_transparency", 0.65, obj)

    print(f"[color_trail_viridis] colored {n} frames, viridis(0) -> viridis(1)")

def draw_arrows_viridis(
    frame_objects,
    displacement_data,
    every=ARROW_EVERY,
    cutoff=ARROW_CUTOFF,
    arrow_scale=ARROW_SCALE,
    shaft_radius=SHAFT_RADIUS,
    head_radius=HEAD_RADIUS,
    color_arrows=COLOR_ARROWS,
    arrow_object="pc1_vectors",
):
    first_obj = frame_objects[0]
    last_obj = frame_objects[-1]

    first_atoms = cmd.get_model(f"{first_obj} and name CA", state=1).atom
    last_atoms = cmd.get_model(f"{last_obj} and name CA", state=1).atom

    max_disp = max(row["displacement_A"] for row in displacement_data)
    min_disp = 0.0
    uniform_rgb = (0.15, 0.15, 0.15)

    arrow_geometry = []
    arrows_drawn = 0

    for idx in range(0, len(first_atoms), every):
        start_atom = first_atoms[idx]
        end_atom = last_atoms[idx]

        start = [float(v) for v in start_atom.coord]
        vector = [float(end_atom.coord[i] - start_atom.coord[i]) for i in range(3)]
        magnitude = math.sqrt(sum(v * v for v in vector))

        if magnitude < cutoff:
            continue

        arrow_tip = [start[i] + arrow_scale * vector[i] for i in range(3)]
        displayed_vector = [arrow_tip[i] - start[i] for i in range(3)]
        displayed_length = math.sqrt(sum(v * v for v in displayed_vector))

        if displayed_length < 1e-8:
            continue

        unit_vector = [v / displayed_length for v in displayed_vector]
        head_length = min(max(0.4, 3 * head_radius), 0.35 * displayed_length)
        shaft_end = [arrow_tip[i] - head_length * unit_vector[i] for i in range(3)]

        rgb = magnitude_rgb_viridis(magnitude, min_disp, max_disp) if color_arrows else uniform_rgb

        arrow_geometry.extend([
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
        ])
        arrows_drawn += 1

    cmd.delete(arrow_object)
    if arrow_geometry:
        cmd.load_cgo(arrow_geometry, arrow_object)

    print(f"[draw_arrows_viridis] drew {arrows_drawn} arrows in '{arrow_object}'")

def run_pipeline():
    try:
        frame_objects = load_and_prepare(TRAJECTORY_FILE, STABLE_SELECTION)
    except RuntimeError as exc:
        print(f"[run_pipeline] {exc}")
        frame_objects = load_and_prepare_fallback(TRAJECTORY_FILE, STABLE_SELECTION)

    displacement_data = compute_displacement(frame_objects, OUTPUT_CSV)
    color_trail_viridis(frame_objects)
    draw_arrows_viridis(frame_objects, displacement_data)

    cmd.bg_color("white")
    cmd.orient()
    cmd.zoom(buffer=3)

    print("[run_pipeline] scene built — no render performed. "
          "Inspect in the PyMOL viewer, then render manually when ready.")


run_pipeline()
