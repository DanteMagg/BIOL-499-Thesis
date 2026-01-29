#!/usr/bin/env python3
"""
MM/GBSA Binding Energy for H-SN1 / TNFR1 Complex
=================================================
- Loads model_01.pdb (chain R = TNFR1, chain A = H-SN1)
- Builds 3 implicit-solvent systems: complex, receptor, peptide
- Uses amber99sb + GBSA (OBC2)
- Minimizes each, computes ΔG_bind
- Calibrates to KD = 32 μM
"""

import numpy as np
from pathlib import Path

from openmm import app, unit, Platform, LangevinMiddleIntegrator
from openmm.app import PDBFile, ForceField, Modeller
import openmm as mm
from pdbfixer import PDBFixer

# Constants
R = 0.001987204  # kcal/(mol·K)
T = 298.15       # K
KD_EXP = 32e-6   # 32 μM in M


def fix_structure(pdb_path: str) -> tuple:
    """Fix PDB and return topology + positions."""
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer.topology, fixer.positions


def extract_chain(topology, positions, chain_id: str) -> tuple:
    """Extract single chain from structure."""
    modeller = Modeller(topology, positions)
    to_delete = [c for c in modeller.topology.chains() if c.id != chain_id]
    modeller.delete(to_delete)
    return modeller.topology, modeller.positions


def create_gbsa_system(topology, forcefield: ForceField):
    """Create system with GB-OBC2 implicit solvent (OpenMM 8.x API)."""
    return forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )


def minimize_energy(topology, positions, forcefield: ForceField, 
                    label: str, max_iter: int = 1000) -> float:
    """Minimize and return energy in kcal/mol."""
    system = create_gbsa_system(topology, forcefield)
    
    # Pick best available platform
    for pname in ['CUDA', 'OpenCL', 'CPU']:
        try:
            platform = Platform.getPlatformByName(pname)
            break
        except:
            continue
    
    integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds)
    sim = app.Simulation(topology, system, integrator, platform)
    sim.context.setPositions(positions)
    
    # Initial energy
    E0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    E0_kcal = E0.value_in_unit(unit.kilocalories_per_mole)
    
    # Minimize
    sim.minimizeEnergy(maxIterations=max_iter)
    
    # Final energy
    state = sim.context.getState(getEnergy=True, getPositions=True)
    E_final = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    print(f"  {label:20s} | Initial: {E0_kcal:10.1f} | Minimized: {E_final:10.1f} kcal/mol")
    return E_final, state.getPositions()


def main():
    # Paths
    project = Path(__file__).parent.parent
    input_pdb = project / "rosie-51675.flex_input_A_R_refinement/output/model_01.pdb"
    output_dir = project / "mmgbsa_results"
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("MM/GBSA BINDING ENERGY CALCULATION")
    print("=" * 70)
    print(f"Input:       {input_pdb.name}")
    print(f"Force field: amber99sb + GB-OBC2")
    print(f"KD (exp):    {KD_EXP*1e6:.0f} μM")
    print("=" * 70)
    
    # Step 1: Fix structure
    print("\n[1] Fixing PDB structure...")
    topology, positions = fix_structure(str(input_pdb))
    
    # Check chains
    chains = {c.id: sum(1 for _ in c.residues()) for c in topology.chains()}
    print(f"    Chains found: {chains}")
    
    # Step 2: Load force field with implicit solvent
    # In OpenMM 8.x, use implicit/*.xml files for GB models
    print("\n[2] Loading amber99sb + implicit/obc2 force field...")
    ff = ForceField('amber99sb.xml', 'implicit/obc2.xml')
    
    # Step 3: Minimize all three systems
    print("\n[3] Minimizing systems...")
    print("-" * 70)
    
    # Complex (R + A)
    E_complex, _ = minimize_energy(topology, positions, ff, "Complex (R+A)")
    
    # Receptor only (chain R)
    topo_R, pos_R = extract_chain(topology, positions, 'R')
    E_receptor, _ = minimize_energy(topo_R, pos_R, ff, "Receptor (R)")
    
    # Peptide only (chain A)
    topo_A, pos_A = extract_chain(topology, positions, 'A')
    E_peptide, _ = minimize_energy(topo_A, pos_A, ff, "Peptide (A)")
    
    # Step 4: Compute binding energy
    print("\n" + "=" * 70)
    print("[4] BINDING ENERGY")
    print("=" * 70)
    
    DG_pred = E_complex - E_receptor - E_peptide
    
    print(f"\n  E(complex)    = {E_complex:12.2f} kcal/mol")
    print(f"  E(receptor)   = {E_receptor:12.2f} kcal/mol")
    print(f"  E(peptide)    = {E_peptide:12.2f} kcal/mol")
    print(f"  {'-'*45}")
    print(f"  ΔG_bind (pred) = {DG_pred:12.2f} kcal/mol")
    
    # Step 5: Experimental ΔG
    print("\n" + "=" * 70)
    print("[5] CALIBRATION TO KD = 32 μM")
    print("=" * 70)
    
    DG_exp = R * T * np.log(KD_EXP)  # ΔG = RT ln(KD)
    
    print(f"\n  KD (experimental) = {KD_EXP*1e6:.0f} μM = {KD_EXP:.2e} M")
    print(f"  ΔG_exp = RT ln(KD) = {DG_exp:.2f} kcal/mol")
    
    # Step 6: Scale factor
    if DG_pred != 0:
        scale_factor = DG_exp / DG_pred
    else:
        scale_factor = float('nan')
    
    print(f"\n  Scale factor (ΔG_exp / ΔG_pred) = {scale_factor:.4f}")
    
    # Interpretation
    print("\n" + "=" * 70)
    print("[6] INTERPRETATION")
    print("=" * 70)
    if DG_pred < 0:
        print("  ✓ Negative ΔG_pred → favorable binding predicted")
    else:
        print("  ⚠ Positive ΔG_pred → unfavorable (check structure)")
    
    print(f"  Typical MM/GBSA scale factors: 0.1–0.3")
    print(f"  Your scale factor: {scale_factor:.4f}")
    
    if 0.05 < abs(scale_factor) < 0.5:
        print("  → Within expected range for MM/GBSA")
    else:
        print("  → Outside typical range; consider structure quality")
    
    # Step 7: Write results
    results_file = output_dir / "mmgbsa_results.txt"
    with open(results_file, 'w') as f:
        f.write("MM/GBSA Binding Energy Results\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input: model_01.pdb\n")
        f.write(f"Chain R: TNFR1 (receptor)\n")
        f.write(f"Chain A: H-SN1 (peptide)\n")
        f.write(f"Force field: amber99sb + GB-OBC2\n")
        f.write(f"Temperature: {T} K\n\n")
        f.write("Component Energies (kcal/mol):\n")
        f.write(f"  E(complex)  = {E_complex:.2f}\n")
        f.write(f"  E(receptor) = {E_receptor:.2f}\n")
        f.write(f"  E(peptide)  = {E_peptide:.2f}\n\n")
        f.write("Predicted Binding Energy:\n")
        f.write(f"  ΔG_bind (MM/GBSA) = {DG_pred:.2f} kcal/mol\n\n")
        f.write(f"Experimental Reference (KD = {KD_EXP*1e6:.0f} μM):\n")
        f.write(f"  ΔG_exp = RT ln(KD) = {DG_exp:.2f} kcal/mol\n\n")
        f.write(f"Scale Factor:\n")
        f.write(f"  ΔG_exp / ΔG_pred = {scale_factor:.4f}\n")
    
    print(f"\n  Results written to: {results_file}")
    print("=" * 70)


if __name__ == "__main__":
    main()
