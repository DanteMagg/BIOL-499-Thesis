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

Research Workflow:

Background
•	IBD affects 3+ million Americans; anti-TNF biologics have 30% non-responder rate
•	Current therapies block both TNFR1 (pro-inflammatory) and TNFR2 (anti-inflammatory)
•	H-SN1: 23-aa sea snake peptide; selectively blocks TNFR1; KD = 32 μM; ~$800/mg synthesis
•	HSN10: 10-aa truncated variant; 11-fold better affinity (KD = 2.8 μM); ~$150/mg—mechanism unknown
Knowledge Gaps
1.	No residue-level structure-activity relationship for H-SN1
2.	HSN10's improvement unexplained mechanistically
3.	No framework integrates binding affinity with manufacturability
4.	TNFR1 polymorphism impact on binding unknown
Aims
1.	Validate computational pipeline; generate baseline H-SN1–TNFR1 complex
2.	Alanine scanning to identify binding hotspots
3.	Design manufacturable variants optimized for affinity and cost
4.	Assess TNFR1 polymorphism impact
Workflow
Step 1: Pipeline Validation
Why: Before we can study H-SN1 binding or design variants, we need a reliable 3D model of the H-SN1–TNFR1 complex. No crystal structure exists, so we must computationally dock the peptide to the receptor and validate our model against published experimental data.
•	Download TNFR1 structure (PDB: 1TNR); prepare for docking
•	Run CABS-dock blind docking for H-SN1 peptide
•	Select highest-density cluster; extract top model
•	Rebuild all-atom structure from Cα trace (CABS outputs coarse-grained only)
•	FlexPepDock refinement (100 low-res + 100 high-res models)
•	Validate: RMSD vs. Zheng et al. (<2 Å); confirm CRD2/CRD3 binding site
•	MM/GBSA binding energy; calibrate against experimental KD = 32 μM
Step 2: Alanine Scanning
Why: We don't know which of H-SN1's 23 residues are critical for binding TNFR1. By mutating each position to alanine and measuring the energy change, we identify "hotspots" (essential residues) vs. "scaffold" (dispensable residues that can be removed). This explains why HSN10's truncation works and guides our variant design.
•	Generate 23 alanine mutants (one per position)
•	Rosetta relax protocol for each mutant complex (5 replicates)
•	Calculate ΔΔG for each position
•	Classify: hotspots (>1.5 kcal/mol), moderate (0.8–1.5), scaffold (<0.8)
•	Generate heat map (Figure 1); explain HSN10 mechanism
Step 3: Variant Design
Why: H-SN1 is too long (expensive to synthesize) and has modest affinity. Using our hotspot map, we can rationally design shorter peptides that keep the critical residues while removing the dispensable ones. We also optimize hotspot residues to improve binding and score each variant for manufacturability.
•	Truncate non-hotspot N/C-terminal residues (target: 10–13 aa)
•	Introduce hotspot mutations to enhance binding
•	Generate 10–12 combined variants
•	Calculate M-scores (manufacturability: length, cysteines, hydrophobicity)
•	AlphaFold for unbound structures; compare to docked conformations
Step 4: Variant Docking & Affinity Prediction
Why: We need to predict how well each designed variant binds TNFR1 and balance this against synthesis cost. Pareto optimization lets us identify variants that are "best of both worlds"—strong binders that are also cheap to make.
•	CABS-dock + FlexPepDock for all variants
•	MM/GBSA binding energy → predicted KD
•	Pareto plot: KD vs. M-score (Figure 4)
•	Select top 3 leads (KD ≤ 5 μM, M-score ≥ 1.6)
Step 5: Polymorphism Analysis
Why: Some patients don't respond to anti-TNF therapy, and genetic variants in TNFR1 may explain this. If our peptide binds differently to polymorphic receptors, we need to either design "promiscuous" variants that work across isoforms or stratify patients by genotype.
•	Literature search: TNFR1 variants linked to IBD/anti-TNF non-response
•	Model 3–5 polymorphic receptors (AlphaFold/homology modeling)
•	Dock H-SN1 WT + leads to each receptor variant
•	Calculate ΔΔGpolymorphism; assess isoform sensitivity
•	Generate polymorphism impact panel (Figure 5)
Step 6: Figure Generation
Why: Clear figures are essential for communicating our findings. Each figure addresses a specific question: which residues matter (Fig 1), where does binding occur (Fig 2), how do variants compare (Fig 3–4), and does receptor variation matter (Fig 5).
•	Figure 1: Alanine scanning heat map
•	Figure 2: H-SN1–TNFR1 docking pose (PyMOL)
•	Figure 3: Variant overlay (WT + leads)
•	Figure 4: Pareto plot
•	Figure 5: Polymorphism impact panel
•	Supplementary Figure S1: AlphaFold validation
Step 7: Thesis Writing
Why: Document all methods, results, and interpretations in a thesis formatted for BIOL 499 requirements.
•	Results (8–10 pages): validation, alanine scan, variants, polymorphisms
•	Methods (6–8 pages): detailed protocols
•	Introduction (4–5 pages): IBD, anti-TNF, H-SN1 background, aims
•	Discussion (5–6 pages): interpretation, limitations, future directions
•	Complete draft: 25–35 pages
Step 8: Revision & Defense
Why: Incorporate feedback from supervisor and second reader to strengthen the thesis, then present and defend the work.
•	Supervisor feedback → revisions
•	Second reader review → final revisions
•	15-slide defense presentation
•	Oral defense + final submission
Expected Outcomes
•	First comprehensive SAR for H-SN1 (5–7 hotspot residues identified)
•	Mechanistic explanation for HSN10's 11-fold affinity improvement
•	2–3 lead variants: KD ≤ 5 μM, 5-fold lower synthesis cost
•	Polymorphism impact assessment for patient stratification

