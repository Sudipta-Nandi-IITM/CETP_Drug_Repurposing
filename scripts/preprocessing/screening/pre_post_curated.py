#!/usr/bin/env python3
"""
figure2a_source_composition.py

Generate a grouped bar plot comparing pre- and post-curation source
composition for FDA, ChEMBL, and DrugBank collections.

Pre-curation values
-------------------
Counted directly from PDB files in the four source-specific folders.

Post-curation values
--------------------
Each occurrence in the post-curation SMILES file is matched against the
source libraries and assigned to one source using an explicit priority:

    FDA > ChEMBL > DrugBank Approved > DrugBank Investigational

Post-curation entries are not deduplicated. Therefore, repeated SMILES
remain repeated occurrences.

Outputs
-------
1. Fig2a_source_composition.png
2. Fig2a_source_composition.pdf
3. Fig2a_source_counts.csv
4. Fig2a_post_provenance_audit.csv
5. Fig2a_unmatched_post_entries.csv, if unmatched entries are detected

Requirements
------------
pandas
numpy
matplotlib
rdkit
"""

from __future__ import annotations

import argparse
import logging
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------
# Source configuration
# ---------------------------------------------------------------------

SOURCE_FOLDERS: dict[str, list[str]] = {
    "FDA": [
        "FDA",
    ],
    "ChEMBL": [
        "CHEMBL",
        "ChEMBL",
    ],
    "DrugBank_Approved": [
        "Drugbank_Approved",
        "DrugBank_Approved",
    ],
    "DrugBank_Investigational": [
        "Drugbank_Investigational",
        "DrugBank_Investigational",
    ],
}

ASSIGNMENT_PRIORITY = [
    "FDA",
    "ChEMBL",
    "DrugBank_Approved",
    "DrugBank_Investigational",
]

PLOT_GROUPS = {
    "FDA": ["FDA"],
    "ChEMBL": ["ChEMBL"],
    "DrugBank": [
        "DrugBank_Approved",
        "DrugBank_Investigational",
    ],
}


# ---------------------------------------------------------------------
# Publication-style plot settings
# ---------------------------------------------------------------------

PRE_COLOR = "#78B7A5"
POST_COLOR = "#8E7DB3"
EDGE_COLOR = "#333333"

FIGURE_SIZE = (8.0, 5.3)
BAR_WIDTH = 0.36

FONT_SIZE_AXIS = 15
FONT_SIZE_TICK = 12
FONT_SIZE_LEGEND = 11
FONT_SIZE_VALUE = 10
FONT_SIZE_TOTAL = 10

OUTPUT_DPI = 600


# ---------------------------------------------------------------------
# General utilities
# ---------------------------------------------------------------------

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Figure 2a: pre- versus post-curation source composition."
    )

    parser.add_argument(
        "--root",
        required=True,
        type=Path,
        help="Root directory containing the four source-specific folders.",
    )

    parser.add_argument(
        "--post-smiles",
        required=True,
        type=Path,
        help="Post-curation SMILES file; the first token of each line is treated as SMILES.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figure2a_outputs"),
        help="Directory for the figure and audit tables.",
    )

    parser.add_argument(
        "--font",
        type=Path,
        default=None,
        help="Optional custom font file, for example timesbold.ttf.",
    )

    parser.add_argument(
        "--pdb-pattern",
        default="*.pdb",
        help="Pattern used to identify pre-curation compound files.",
    )

    parser.add_argument(
        "--match-mode",
        choices=["connectivity", "full_inchikey", "smiles"],
        default="connectivity",
        help=(
            "Molecular key used for source matching. "
            "'connectivity' uses the first InChIKey block and is tolerant "
            "of stereochemical differences; 'full_inchikey' is stricter."
        ),
    )

    parser.add_argument(
        "--allow-unmatched",
        action="store_true",
        help="Generate the figure even when post-curation entries remain unmatched.",
    )

    return parser.parse_args()


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
    )


def configure_font(font_path: Path | None) -> None:
    if font_path is not None:
        if not font_path.is_file():
            warnings.warn(f"Font file was not found: {font_path}")
        else:
            font_manager.fontManager.addfont(str(font_path))
            family = font_manager.FontProperties(
                fname=str(font_path)
            ).get_name()
            matplotlib.rcParams["font.family"] = family

    matplotlib.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": FONT_SIZE_AXIS,
            "xtick.labelsize": FONT_SIZE_TICK,
            "ytick.labelsize": FONT_SIZE_TICK,
            "legend.fontsize": FONT_SIZE_LEGEND,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
        }
    )


# ---------------------------------------------------------------------
# Molecular parsing and matching
# ---------------------------------------------------------------------

def load_pdb_molecule(path: Path) -> Chem.Mol | None:
    """
    Read a PDB file using RDKit.

    Bond perception from PDB coordinates may occasionally fail. Files that
    cannot be converted to a molecular key are reported in the audit output.
    """
    try:
        mol = Chem.MolFromPDBFile(
            str(path),
            sanitize=False,
            removeHs=False,
        )
    except Exception:
        return None

    if mol is None:
        return None

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        # A partially sanitized molecule can still be usable for generating
        # a non-stereochemical structural key.
        try:
            operations = (
                Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            )
            Chem.SanitizeMol(mol, sanitizeOps=operations)
        except Exception:
            pass

    return mol


def molecular_key(
    mol: Chem.Mol,
    match_mode: str,
) -> str | None:
    """
    Convert an RDKit molecule to the selected matching key.
    """
    if mol is None:
        return None

    if match_mode in {"connectivity", "full_inchikey"}:
        try:
            key = inchi.MolToInchiKey(mol)
            if key:
                if match_mode == "connectivity":
                    return key.split("-")[0]
                return key
        except Exception:
            pass

    try:
        return Chem.MolToSmiles(
            mol,
            canonical=True,
            isomericSmiles=(match_mode != "connectivity"),
        )
    except Exception:
        return None


def smiles_to_key(
    smiles: str,
    match_mode: str,
) -> str | None:
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None

    if mol is None:
        return None

    return molecular_key(mol, match_mode)


# ---------------------------------------------------------------------
# Pre-curation processing
# ---------------------------------------------------------------------

def locate_source_folder(
    root: Path,
    aliases: list[str],
) -> Path | None:
    """
    Locate a source folder using case-insensitive matching.
    """
    available = {
        path.name.lower(): path
        for path in root.iterdir()
        if path.is_dir()
    }

    for alias in aliases:
        candidate = available.get(alias.lower())
        if candidate is not None:
            return candidate

    return None


def process_precuration_sources(
    root: Path,
    pdb_pattern: str,
    match_mode: str,
) -> tuple[dict[str, int], dict[str, set[str]], pd.DataFrame]:
    """
    Count source files and generate molecular-key sets for provenance matching.
    """
    source_counts: dict[str, int] = {}
    source_keys: dict[str, set[str]] = {}
    audit_rows: list[dict[str, object]] = []

    for source, aliases in SOURCE_FOLDERS.items():
        folder = locate_source_folder(root, aliases)

        if folder is None:
            raise FileNotFoundError(
                f"No folder found for {source}. Tried aliases: {aliases}"
            )

        pdb_files = sorted(folder.glob(pdb_pattern))
        source_counts[source] = len(pdb_files)

        keys: set[str] = set()
        failed = 0

        for pdb_file in pdb_files:
            mol = load_pdb_molecule(pdb_file)
            key = molecular_key(mol, match_mode) if mol else None

            if key is None:
                failed += 1
            else:
                keys.add(key)

            audit_rows.append(
                {
                    "source": source,
                    "folder": str(folder),
                    "file": pdb_file.name,
                    "match_key": key,
                    "key_generated": key is not None,
                }
            )

        source_keys[source] = keys

        logging.info(
            "%-28s files=%5d | molecular keys=%5d | failed=%d",
            source,
            len(pdb_files),
            len(keys),
            failed,
        )

    return source_counts, source_keys, pd.DataFrame(audit_rows)


# ---------------------------------------------------------------------
# Post-curation processing
# ---------------------------------------------------------------------

def read_postcuration_smiles(
    path: Path,
    match_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read occurrence-level SMILES.

    No deduplication is performed.
    """
    valid_rows: list[dict[str, object]] = []
    invalid_rows: list[dict[str, object]] = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.strip().split()

            if not fields:
                continue

            smiles = fields[0]
            key = smiles_to_key(smiles, match_mode)

            row = {
                "occurrence_id": len(valid_rows) + len(invalid_rows) + 1,
                "line_number": line_number,
                "input_smiles": smiles,
                "match_key": key,
            }

            if key is None:
                invalid_rows.append(row)
            else:
                valid_rows.append(row)

    valid_df = pd.DataFrame(valid_rows)
    invalid_df = pd.DataFrame(invalid_rows)

    logging.info(
        "Post-curation entries: valid=%d | invalid=%d",
        len(valid_df),
        len(invalid_df),
    )

    return valid_df, invalid_df


def assign_postcuration_sources(
    post_df: pd.DataFrame,
    source_keys: dict[str, set[str]],
) -> tuple[pd.DataFrame, dict[str, int]]:
    """
    Assign every post-curation occurrence to one source.

    When an entry matches multiple sources, the first source in
    ASSIGNMENT_PRIORITY is used. All matching sources are retained in the
    audit table.
    """
    assigned_counts = {
        source: 0
        for source in SOURCE_FOLDERS
    }

    audit_rows: list[dict[str, object]] = []

    for row in post_df.itertuples(index=False):
        matches = [
            source
            for source in ASSIGNMENT_PRIORITY
            if row.match_key in source_keys.get(source, set())
        ]

        assigned_source = matches[0] if matches else None

        if assigned_source is not None:
            assigned_counts[assigned_source] += 1

        audit_rows.append(
            {
                "occurrence_id": row.occurrence_id,
                "line_number": row.line_number,
                "input_smiles": row.input_smiles,
                "match_key": row.match_key,
                "matched_sources": ";".join(matches),
                "number_of_source_matches": len(matches),
                "assigned_source": assigned_source,
                "assignment_rule": (
                    "first match in priority order"
                    if len(matches) > 1
                    else "single source match"
                    if len(matches) == 1
                    else "unmatched"
                ),
            }
        )

    return pd.DataFrame(audit_rows), assigned_counts


# ---------------------------------------------------------------------
# Count aggregation
# ---------------------------------------------------------------------

def combine_plot_groups(
    detailed_counts: dict[str, int],
) -> dict[str, int]:
    combined: dict[str, int] = {}

    for plot_label, detailed_sources in PLOT_GROUPS.items():
        combined[plot_label] = sum(
            detailed_counts.get(source, 0)
            for source in detailed_sources
        )

    return combined


def build_summary_table(
    pre_grouped: dict[str, int],
    post_grouped: dict[str, int],
    valid_post_count: int,
    unmatched_count: int,
) -> pd.DataFrame:
    rows = []

    for source in PLOT_GROUPS:
        rows.append(
            {
                "source": source,
                "pre_curation_count": pre_grouped[source],
                "post_curation_count": post_grouped[source],
            }
        )

    rows.append(
        {
            "source": "TOTAL_MATCHED",
            "pre_curation_count": sum(pre_grouped.values()),
            "post_curation_count": sum(post_grouped.values()),
        }
    )

    rows.append(
        {
            "source": "TOTAL_VALID_POST_ENTRIES",
            "pre_curation_count": np.nan,
            "post_curation_count": valid_post_count,
        }
    )

    rows.append(
        {
            "source": "UNMATCHED_POST_ENTRIES",
            "pre_curation_count": np.nan,
            "post_curation_count": unmatched_count,
        }
    )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------

def annotate_bars(
    ax: plt.Axes,
    bars,
    vertical_padding: float,
) -> None:
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + vertical_padding,
            f"{int(height):,}",
            ha="center",
            va="bottom",
            fontsize=FONT_SIZE_VALUE,
        )


def plot_source_composition(
    pre_grouped: dict[str, int],
    post_grouped: dict[str, int],
    output_directory: Path,
    valid_post_count: int,
    unmatched_count: int,
) -> None:
    labels = list(PLOT_GROUPS.keys())
    pre_values = [pre_grouped[label] for label in labels]
    post_values = [post_grouped[label] for label in labels]

    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    pre_bars = ax.bar(
        x - BAR_WIDTH / 2,
        pre_values,
        width=BAR_WIDTH,
        label="Pre-curation",
        color=PRE_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=0.7,
    )

    post_bars = ax.bar(
        x + BAR_WIDTH / 2,
        post_values,
        width=BAR_WIDTH,
        label="Post-curation",
        color=POST_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=0.7,
    )

    maximum = max(pre_values + post_values)
    annotation_padding = maximum * 0.018

    annotate_bars(ax, pre_bars, annotation_padding)
    annotate_bars(ax, post_bars, annotation_padding)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Compound count")

    ax.legend(
        frameon=False,
        loc="upper right",
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    ax.set_ylim(0, maximum * 1.17)

    computed_pre_total = sum(pre_values)
    matched_post_total = sum(post_values)

    totals_text = (
        f"Pre-curation total: {computed_pre_total:,}\n"
        f"Post-curation total: {valid_post_count:,}"
    )

    if unmatched_count:
        totals_text += (
            f"\nMatched post entries: {matched_post_total:,}"
            f"\nUnmatched post entries: {unmatched_count:,}"
        )

    ax.text(
        0.02,
        0.98,
        totals_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FONT_SIZE_TOTAL,
    )

    ax.text(
        -0.08,
        1.02,
        "(a)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=14,
        fontweight="bold",
    )

    fig.tight_layout()

    png_path = output_directory / "Fig2a_source_composition.png"
    pdf_path = output_directory / "Fig2a_source_composition.pdf"
    svg_path = output_directory / "Fig2a_source_composition.svg"

    fig.savefig(png_path, dpi=OUTPUT_DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")

    plt.close(fig)

    logging.info("Saved %s", png_path)
    logging.info("Saved %s", pdf_path)
    logging.info("Saved %s", svg_path)


# ---------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------

def main() -> None:
    args = parse_arguments()
    configure_logging()
    configure_font(args.font)

    if not args.root.is_dir():
        raise NotADirectoryError(
            f"Root directory does not exist: {args.root}"
        )

    if not args.post_smiles.is_file():
        raise FileNotFoundError(
            f"Post-curation SMILES file does not exist: {args.post_smiles}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Process pre-curation folders
    pre_counts, source_keys, pre_audit = process_precuration_sources(
        root=args.root,
        pdb_pattern=args.pdb_pattern,
        match_mode=args.match_mode,
    )

    # 2. Read post-curation SMILES without deduplication
    post_valid, post_invalid = read_postcuration_smiles(
        path=args.post_smiles,
        match_mode=args.match_mode,
    )

    # 3. Match post-curation occurrences to source libraries
    provenance_audit, post_counts = assign_postcuration_sources(
        post_df=post_valid,
        source_keys=source_keys,
    )

    unmatched = provenance_audit[
        provenance_audit["assigned_source"].isna()
    ].copy()

    multi_source = provenance_audit[
        provenance_audit["number_of_source_matches"] > 1
    ].copy()

    # 4. Combine approved and investigational DrugBank counts
    pre_grouped = combine_plot_groups(pre_counts)
    post_grouped = combine_plot_groups(post_counts)

    # 5. Verify totals
    matched_post_count = sum(post_grouped.values())
    valid_post_count = len(post_valid)
    unmatched_count = len(unmatched)

    if matched_post_count + unmatched_count != valid_post_count:
        raise RuntimeError(
            "Internal count check failed: matched + unmatched does not "
            "equal the number of valid post-curation entries."
        )

    if len(post_invalid):
        invalid_path = args.output_dir / "Fig2a_invalid_post_entries.csv"
        post_invalid.to_csv(invalid_path, index=False)
        logging.warning(
            "%d post-curation lines were invalid. See %s",
            len(post_invalid),
            invalid_path,
        )

    if unmatched_count and not args.allow_unmatched:
        unmatched_path = (
            args.output_dir / "Fig2a_unmatched_post_entries.csv"
        )
        unmatched.to_csv(unmatched_path, index=False)

        raise RuntimeError(
            f"{unmatched_count} valid post-curation entries could not be "
            f"matched to a source. Inspect {unmatched_path}. "
            "Use --allow-unmatched only after reviewing these entries."
        )

    # 6. Save audit tables
    pre_audit.to_csv(
        args.output_dir / "Fig2a_precuration_key_audit.csv",
        index=False,
    )

    provenance_audit.to_csv(
        args.output_dir / "Fig2a_post_provenance_audit.csv",
        index=False,
    )

    if unmatched_count:
        unmatched.to_csv(
            args.output_dir / "Fig2a_unmatched_post_entries.csv",
            index=False,
        )

    if len(multi_source):
        multi_source.to_csv(
            args.output_dir / "Fig2a_multisource_post_entries.csv",
            index=False,
        )

    summary = build_summary_table(
        pre_grouped=pre_grouped,
        post_grouped=post_grouped,
        valid_post_count=valid_post_count,
        unmatched_count=unmatched_count,
    )

    summary.to_csv(
        args.output_dir / "Fig2a_source_counts.csv",
        index=False,
    )

    # 7. Plot
    plot_source_composition(
        pre_grouped=pre_grouped,
        post_grouped=post_grouped,
        output_directory=args.output_dir,
        valid_post_count=valid_post_count,
        unmatched_count=unmatched_count,
    )

    logging.info("Computed pre-curation total: %d", sum(pre_grouped.values()))
    logging.info("Valid post-curation entries: %d", valid_post_count)
    logging.info("Matched post-curation entries: %d", matched_post_count)
    logging.info("Unmatched post-curation entries: %d", unmatched_count)
    logging.info("Multiple-source matches: %d", len(multi_source))


if __name__ == "__main__":
    main()
