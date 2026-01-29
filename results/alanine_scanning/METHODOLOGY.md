# Alanine Scanning Methodology

## Purpose
Identify binding hotspots in H-SN1 peptide (23 residues) by computational alanine scanning mutagenesis against TNFR1 receptor.

## Theoretical Background

### MM/GBSA Binding Energy
The Molecular Mechanics/Generalized Born Surface Area (MM/GBSA) method estimates binding free energy:

```
ΔG_bind = E_complex − E_receptor − E_peptide
```

Where each energy term includes:
- Molecular mechanics (bonds, angles, dihedrals, van der Waals, electrostatics)
- Generalized Born solvation (polar contribution)
- Surface area term (nonpolar contribution)

### Alanine Scanning
Alanine is the "null" mutation - it removes the sidechain while preserving backbone geometry. By comparing:

```
ΔΔG = ΔG_bind(mutant) − ΔG_bind(WT)
```

We identify residues critical for binding:
- **Hotspot** (ΔΔG > 1.5 kcal/mol): Essential for binding
- **Moderate** (0.8–1.5 kcal/mol): Contributes to binding
- **Scaffold** (< 0.8 kcal/mol): Dispensable for binding

## Computational Protocol

### Step 1: Structure Preparation
- **Input**: FlexPepDock-refined complex (`model_01.pdb`)
- **Tool**: PDBFixer (OpenMM)
- **Operations**:
  - Add missing atoms
  - Add hydrogens at pH 7.0
  - Replace nonstandard residues

### Step 2: Mutant Generation
- **Tool**: Biopython
- **Method**: For each position 1–23 in chain A:
  - Remove sidechain atoms (except Cβ)
  - Rename residue to ALA
  - Preserve chain R (TNFR1) unchanged
- **Output**: 23 mutant PDB files

### Step 3: Structure Relaxation
- **Tool**: OpenMM (substitute for Rosetta relax)
- **Force Field**: AMBER99SB + GB-OBC2 implicit solvent
- **Protocol** (per structure, 5 independent replicates):
  1. Energy minimization (500 steps)
  2. Heating: 0→300K over 5 ps (NVT, 1 fs timestep)
  3. Equilibration: 300K for 10 ps (NVT)
  4. Final minimization (1000 steps)
- **Purpose**: Remove steric clashes, allow local relaxation
- **Random seeds**: Different for each replicate to sample conformational diversity

### Step 4: MM/GBSA Calculation
- **Force Field**: AMBER99SB
- **Implicit Solvent**: GB-OBC2 (Onufriev-Bashford-Case model II)
- **Salt concentration**: 0.15 M (physiological)
- **Dielectric constants**: solute=1.0, solvent=78.5

For each relaxed structure:
1. Compute E_complex (chains R+A together)
2. Extract chain R → compute E_receptor
3. Extract chain A → compute E_peptide
4. ΔG_bind = E_complex − E_receptor − E_peptide

### Step 5: Scaling Calibration
Raw MM/GBSA overestimates binding energy due to neglected entropy terms.

**Calibration to experimental KD**:
- H-SN1 KD = 32 μM (Zheng et al.)
- ΔG_exp = RT ln(KD) = 0.00199 × 298.15 × ln(32×10⁻⁶) = −6.14 kcal/mol
- Scale factor = ΔG_exp / ΔG_pred(WT)

### Step 6: ΔΔG Calculation
```
ΔΔG_i = ΔG_scaled(mutant_i) − ΔG_scaled(WT)
```

Positive ΔΔG indicates the mutation destabilizes binding (residue is important).

### Step 7: Classification
| Class    | ΔΔG Range      | Interpretation                    |
|----------|----------------|-----------------------------------|
| Hotspot  | > 1.5 kcal/mol | Essential for binding             |
| Moderate | 0.8–1.5        | Contributes to binding            |
| Scaffold | < 0.8          | Dispensable; candidate for removal|

## Software Versions
- Python: 3.12
- OpenMM: 8.4.0
- PDBFixer: 1.12.0
- Biopython: (system version)
- NumPy: (system version)
- Matplotlib: (system version)

## References
1. Onufriev A, Bashford D, Case DA. Exploring protein native states and large-scale conformational changes with a modified generalized born model. Proteins. 2004;55(2):383-394.
2. Kollman PA, et al. Calculating structures and free energies of complex molecules: combining molecular mechanics and continuum models. Acc Chem Res. 2000;33(12):889-897.
3. Zheng Y, et al. H-SN1 peptide selectively inhibits TNFR1 signaling. (KD = 32 μM reference)

## Output Files
- `relaxed_structures/` - All relaxed PDB files
- `energy_log.csv` - Raw energies for each replicate
- `alanine_scan_results.csv` - Final ΔΔG values and classifications
- `alanine_scan_heatmap.png` - Figure 1: Binding hotspot visualization
- `run_log.txt` - Complete execution log with timestamps

