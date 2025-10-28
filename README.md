# BIOL 499 Thesis

## Project Title
Structural and Functional Analysis of Venom Proteins (H-SN1, Crotoxin) in Inflammation Pathways

## Overview
This project explores the structure and molecular interactions of venom-derived proteins (H-SN1, Crotoxin) with human inflammation-related targets (TNFR1, TNF-α, IL-1β), using computational tools for structure prediction, molecular docking, and visualization.

## Directory Structure
- data/sequences/: FASTA files for all proteins and targets
    - e.g., HSN1.fasta, Crotoxin.fasta, TNFR1.fasta
- data/structures/: PDB files (downloaded or predicted)
    - e.g., TNFR1_AlphaFold2_2025-10-28.pdb
- notebooks/: Jupyter/Colab notebooks for analysis
- scripts/: Custom Python or shell scripts for processing and automation
- docs/: Documentation, protocols, and literature notes
- results/: Docking results, summary files, figures

## Labeling Conventions
- FASTA files: {Protein}.fasta (e.g., HSN1.fasta)
- PDB files: {Protein}_{method}_{YYYY-MM-DD}.pdb (e.g., Crotoxin_AlphaFold2_2025-10-28.pdb)
- Docking results: {Ligand}_{Target}_vina_{YYYY-MM-DD}.pdbqt
- Notebooks: {AnalysisType}_{Date}.ipynb

## Getting Started
1. All software dependencies and installation instructions are listed in /docs/installation.md.
2. Raw sequences and structures are under /data/.
3. Analysis code and protocols are in /notebooks/ and /scripts/.

## Author
Dante Maggiotto  
BIOL 499 Thesis, University of Waterloo

---

For issues or questions about this project, contact Dante at dante.maggiotto@gmail.com
