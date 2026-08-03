# CETP_drug_repurposing: Structure-Guided Virtual Screening for CETP Inhibitor Repurposing

A multi-stage, structure-based drug-repurposing pipeline that combines molecular docking, ligand strain filtering, consensus clustering, molecular dynamics (MD), MM-PBSA free-energy analysis, and biochemical validation to identify repurposed small-molecule inhibitors of Cholesteryl Ester Transfer Protein (CETP).

Paper: [DOI link — add once published]
Preprint: [bioRxiv/ChemRxiv link — add if posted]
Data (Zenodo): [DOI badge — add once archived]

![Graphical abstract](assets/graphical_abstract.svg)

## Table of Contents
1. [Abstract](#abstract)
2. [Repository Structure](#repository-structure)
3. [Environment Setup](#environment-setup)
4. [Usage](#usage)
5. [Results](#results)
6. [Citing](#citing)
7. [License](#license)

## Abstract

CETP is a clinically validated but difficult target for cardiovascular drug discovery. This repository contains the full computational pipeline used to identify repurposed small-molecule CETP inhibitors from a curated library of ~6,500 FDA-approved, investigational, and clinical-trial compounds. The workflow proceeds through rigid and flexible docking against the CETP tunnel, ligand-strain filtering, consensus similarity-based clustering to enforce chemotype diversity, ADME/PAINS filtering, all-atom MD simulation of shortlisted complexes, MM-PBSA binding free-energy decomposition, and free-energy-landscape analysis. Two candidates, Adapalene and Buclizine, were experimentally confirmed to inhibit CETP activity in a concentration-dependent manner in a fluorescence-based biochemical assay.

## Repository Structure

Large/raw outputs (docking pose libraries, MD trajectories, MM-PBSA tables, raw assay reads) live on **Zenodo** — see [Results](#results). This repo holds the code needed to reproduce them.

```
CETP-VS/
├── README.md
├── LICENSE
├── environment.yml
├── CITATION.cff
├── requirements.txt
├── CONTRIBUTING.md
├── assets/
│   └── graphical_abstract.svg
├── data/
│   └── README.md                  # pointers to the Zenodo record; no large files here
├── scripts/
│   ├── 01_library_curation/
│   │   ├── merge_fda_drugbank_chembl.py
│   │   └── standardize_dedupe.py
│   ├── 02_conformer_generation/
│   │   └── generate_conformers.py       # Open Babel weighted rotor search
│   ├── 03_docking/
│   │   ├── screen1_rigid_docking.py
│   │   └── screen2_flexible_docking.py
│   ├── 04_strain_filter/
│   │   └── torsional_strain_filter.py
│   ├── 05_clustering/
│   │   └── consensus_clustering.py      # 2D fingerprint + 3D shape/electrostatic
│   ├── 06_admet_pains_filter/
│   │   └── admet_pains_filter.py
│   ├── 07_md_simulation/
│   │   ├── prep_topology.py             # ATB-generated GROMOS topologies
│   │   └── run_md.sh                    # GROMACS 2021.1 production run
│   ├── 08_mmpbsa/
│   │   └── run_mmpbsa.py                # g_mmpbsa wrapper + FEL projection
│   └── 09_assay_analysis/
│       └── plot_inhibition_curves.py
└── notebooks/
    └── figures.ipynb                    # regenerates manuscript figures
```

## Environment Setup

```bash
conda env create --name cetp-vs --file=environment.yml
conda activate cetp-vs
```

`environment.yml` should pin: RDKit, Open Babel, AutoDock Vina, scikit-learn, GROMACS (or a note pointing to a separate cluster install, since GROMACS is usually built from source on HPC), MDAnalysis, g_mmpbsa, pandas, matplotlib.

## Usage

```bash
# 0. Curate and standardize the FDA + DrugBank + ChEMBL library
python scripts/01_library_curation/merge_fda_drugbank_chembl.py
python scripts/01_library_curation/standardize_dedupe.py

# 1. Generate low-energy 3D conformers
python scripts/02_conformer_generation/generate_conformers.py

# 2. Screen-1: rigid docking against the Torcetrapib-defined CETP pocket
python scripts/03_docking/screen1_rigid_docking.py

# 3. Screen-2: flexible docking with local side-chain sampling
python scripts/03_docking/screen2_flexible_docking.py

# 4. Torsional strain filtering (5 TEU cutoff)
python scripts/04_strain_filter/torsional_strain_filter.py

# 5. Consensus clustering into chemotype families
python scripts/05_clustering/consensus_clustering.py

# 6. ADME/PAINS translational filtering
python scripts/06_admet_pains_filter/admet_pains_filter.py

# 7. MD simulation of shortlisted complexes (250 ns, GROMACS)
bash scripts/07_md_simulation/run_md.sh

# 8. MM-PBSA binding free energy + free-energy-landscape analysis
python scripts/08_mmpbsa/run_mmpbsa.py

# 9. Biochemical assay analysis and plotting
python scripts/09_assay_analysis/plot_inhibition_curves.py
```

## Data availability

The complete computational and experimental dataset is available from Zenodo:

**Dataset DOI:** [insert Zenodo DOI]

- **docking/** — rigid and flexible docking poses and scores for all screened conformers
- **clustering/** — similarity matrices, cluster assignments, and consensus scores
- **md_trajectories/** — GROMACS input/output for the five lead complexes plus the Torcetrapib control (250 ns each)
- **mmpbsa/** — per-complex and per-residue binding free-energy tables
- **assay_data/** — raw fluorescence reads and normalized inhibition data for Adapalene, Buclizine, and Mefloquine
- **figures/** — source data and scripts used to generate all manuscript figures

## Installation

```bash
git clone https://github.com/Sudipta-Nandi-IITM/CETP_Drug_Repurposing.git
cd CETP_Drug_Repurposing
conda env create -f environment.yml
conda activate cetp-repurposing

## Citing

If this pipeline is useful for your work, please cite:

```bibtex
@article{nandi2026cetp,
  author  = {Nandi, Sudipta and Senapati, Sanjib},
  title   = {Structure-guided Exploration of Repurposed Small Molecule Inhibitors Targeting Cholesteryl Ester Transfer Protein},
  journal = {[Journal name — add on acceptance]},
  year    = {2026},
  doi     = {[DOI — add on publication]}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Contact

Sudipta Nandi — Computational Biophysics Lab, Dept. of Biotechnology, IIT Madras ([nandi.sudipta5997@gmail.com](mailto:nandi.sudipta5997@gmail.com))
