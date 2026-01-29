#!/usr/bin/env python3
"""
Alanine Scanning ΔΔG Calculation Pipeline
==========================================
BIOL 499 Thesis - Dante Maggiotto

Purpose:
    Calculate binding energy changes (ΔΔG) for each alanine mutant of H-SN1
    peptide bound to TNFR1 receptor to identify binding hotspots.

Protocol:
    1. Load WT and 23 mutant structures
    2. Relax each structure (5 replicates) using OpenMM
    3. Calculate MM/GBSA binding energy for each relaxed structure
    4. Average over replicates
    5. Calibrate to experimental KD = 32 μM
    6. Calculate ΔΔG = ΔG(mutant) - ΔG(WT)
    7. Classify residues as hotspot/moderate/scaffold

Input:
    - Wild-type: model_01.pdb (FlexPepDock refined)
    - Mutants: model_01_Ala01.pdb ... model_01_Ala23.pdb
    
Output:
    - alanine_scan_results.csv: Complete results table
    - alanine_scan_heatmap.png: Figure 1 visualization
    - energy_log.csv: Raw energy data for all replicates
    - run_log.txt: Detailed execution log

Author: Dante Maggiotto
Date: January 2026
"""

import os
import sys
import csv
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple
import numpy as np

# OpenMM imports
from openmm import app, unit, Platform, LangevinMiddleIntegrator
from openmm.app import PDBFile, ForceField, Modeller, Simulation
import openmm as mm
from pdbfixer import PDBFixer

# Plotting
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# =============================================================================
# CONSTANTS
# =============================================================================

# Thermodynamic constants
R = 0.001987204  # Gas constant in kcal/(mol·K)
T = 298.15       # Temperature in Kelvin (25°C)
KD_EXP = 32e-6   # Experimental KD in M (32 μM)

# Experimental binding free energy
DG_EXP = R * T * np.log(KD_EXP)  # = -6.14 kcal/mol

# Classification thresholds (kcal/mol)
HOTSPOT_THRESHOLD = 1.5
MODERATE_THRESHOLD = 0.8

# Relaxation parameters
N_REPLICATES = 5
HEATING_STEPS = 5000      # 5 ps at 1 fs timestep
EQUILIBRATION_STEPS = 10000  # 10 ps
MINIMIZE_STEPS_INITIAL = 500
MINIMIZE_STEPS_FINAL = 1000

# H-SN1 sequence (for reference)
HSN1_SEQUENCE = "DEQHLETELHTHLTSVLTANGFQ"

# =============================================================================
# LOGGING SETUP
# =============================================================================

def setup_logging(output_dir: Path) -> logging.Logger:
    """Configure logging to both file and console."""
    log_file = output_dir / "run_log.txt"
    
    # Create logger
    logger = logging.getLogger("AlanineScan")
    logger.setLevel(logging.DEBUG)
    
    # File handler (detailed)
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.DEBUG)
    fh_format = logging.Formatter('%(asctime)s | %(levelname)-8s | %(message)s',
                                   datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(fh_format)
    
    # Console handler (info only)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch_format = logging.Formatter('%(message)s')
    ch.setFormatter(ch_format)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

# =============================================================================
# STRUCTURE HANDLING
# =============================================================================

def fix_structure(pdb_path: str, logger: logging.Logger) -> Tuple:
    """
    Prepare PDB structure for simulation using PDBFixer.
    
    Operations:
        - Find and add missing residues
        - Replace nonstandard residues
        - Add missing atoms
        - Add hydrogens at pH 7.0
    
    Args:
        pdb_path: Path to input PDB file
        logger: Logger instance
        
    Returns:
        tuple: (topology, positions)
    """
    logger.debug(f"Fixing structure: {pdb_path}")
    
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    
    return fixer.topology, fixer.positions


def extract_chain(topology, positions, chain_id: str) -> Tuple:
    """
    Extract a single chain from a multi-chain structure.
    
    Args:
        topology: OpenMM Topology object
        positions: OpenMM Positions object
        chain_id: Chain identifier ('R' or 'A')
        
    Returns:
        tuple: (new_topology, new_positions)
    """
    modeller = Modeller(topology, positions)
    chains_to_delete = [c for c in modeller.topology.chains() if c.id != chain_id]
    modeller.delete(chains_to_delete)
    return modeller.topology, modeller.positions

# =============================================================================
# RELAXATION PROTOCOL
# =============================================================================

def relax_structure(topology, positions, forcefield: ForceField, 
                    seed: int, logger: logging.Logger) -> Tuple:
    """
    Relax structure using OpenMM MD protocol.
    
    Protocol:
        1. Initial minimization (500 steps)
        2. Heating from 0K to 300K over 5 ps (NVT)
        3. Equilibration at 300K for 10 ps (NVT)
        4. Final minimization (1000 steps)
    
    Args:
        topology: OpenMM Topology
        positions: Initial positions
        forcefield: ForceField object
        seed: Random seed for reproducibility
        logger: Logger instance
        
    Returns:
        tuple: (final_positions, final_energy)
    """
    # Create system with implicit solvent
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )
    
    # Get platform
    platform = Platform.getPlatformByName('CPU')
    
    # Initial minimization
    logger.debug(f"  Initial minimization ({MINIMIZE_STEPS_INITIAL} steps)...")
    integrator = LangevinMiddleIntegrator(10*unit.kelvin, 1/unit.picosecond, 
                                           1*unit.femtoseconds)
    integrator.setRandomNumberSeed(seed)
    sim = Simulation(topology, system, integrator, platform)
    sim.context.setPositions(positions)
    sim.minimizeEnergy(maxIterations=MINIMIZE_STEPS_INITIAL)
    
    # Get minimized positions
    state = sim.context.getState(getPositions=True)
    positions = state.getPositions()
    
    # Heating phase (0K -> 300K)
    logger.debug(f"  Heating 0→300K ({HEATING_STEPS} steps)...")
    del sim, integrator
    
    integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond,
                                           1*unit.femtoseconds)
    integrator.setRandomNumberSeed(seed + 1)
    sim = Simulation(topology, system, integrator, platform)
    sim.context.setPositions(positions)
    sim.context.setVelocitiesToTemperature(10*unit.kelvin, seed + 2)
    
    # Gradually heat
    for i in range(10):
        target_temp = 30 * (i + 1)
        sim.context.setVelocitiesToTemperature(target_temp * unit.kelvin, seed + i)
        sim.step(HEATING_STEPS // 10)
    
    # Equilibration at 300K
    logger.debug(f"  Equilibration at 300K ({EQUILIBRATION_STEPS} steps)...")
    sim.step(EQUILIBRATION_STEPS)
    
    # Get equilibrated positions
    state = sim.context.getState(getPositions=True)
    positions = state.getPositions()
    
    # Final minimization
    logger.debug(f"  Final minimization ({MINIMIZE_STEPS_FINAL} steps)...")
    sim.minimizeEnergy(maxIterations=MINIMIZE_STEPS_FINAL)
    
    # Get final state
    state = sim.context.getState(getPositions=True, getEnergy=True)
    final_positions = state.getPositions()
    final_energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    return final_positions, final_energy

# =============================================================================
# MM/GBSA CALCULATION
# =============================================================================

def calculate_mmgbsa(topology, positions, forcefield: ForceField, 
                     logger: logging.Logger) -> Dict[str, float]:
    """
    Calculate MM/GBSA binding energy components.
    
    Computes:
        E_complex: Energy of receptor-peptide complex
        E_receptor: Energy of receptor alone (chain R)
        E_peptide: Energy of peptide alone (chain A)
        dG_bind: E_complex - E_receptor - E_peptide
    
    Args:
        topology: Complex topology
        positions: Complex positions
        forcefield: ForceField with implicit solvent
        logger: Logger instance
        
    Returns:
        dict: Energy components in kcal/mol
    """
    def get_energy(topo, pos):
        """Get minimized potential energy."""
        system = forcefield.createSystem(topo, nonbondedMethod=app.NoCutoff,
                                         constraints=app.HBonds)
        integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond,
                                               2*unit.femtoseconds)
        platform = Platform.getPlatformByName('CPU')
        sim = Simulation(topo, system, integrator, platform)
        sim.context.setPositions(pos)
        sim.minimizeEnergy(maxIterations=100)
        state = sim.context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    # Complex energy
    E_complex = get_energy(topology, positions)
    
    # Receptor energy (chain R)
    topo_R, pos_R = extract_chain(topology, positions, 'R')
    E_receptor = get_energy(topo_R, pos_R)
    
    # Peptide energy (chain A)
    topo_A, pos_A = extract_chain(topology, positions, 'A')
    E_peptide = get_energy(topo_A, pos_A)
    
    # Binding energy
    dG_bind = E_complex - E_receptor - E_peptide
    
    return {
        'E_complex': E_complex,
        'E_receptor': E_receptor,
        'E_peptide': E_peptide,
        'dG_bind': dG_bind
    }

# =============================================================================
# MAIN PIPELINE
# =============================================================================

def process_structure(name: str, pdb_path: Path, forcefield: ForceField,
                      output_dir: Path, logger: logging.Logger) -> Dict:
    """
    Process a single structure through the full pipeline.
    
    Args:
        name: Structure identifier (e.g., 'WT', 'Ala01')
        pdb_path: Path to PDB file
        forcefield: ForceField object
        output_dir: Directory for relaxed structures
        logger: Logger instance
        
    Returns:
        dict: Results including all replicate energies and average
    """
    logger.info(f"Processing {name}...")
    
    # Fix structure
    topology, positions = fix_structure(str(pdb_path), logger)
    
    # Run replicates
    replicate_results = []
    
    for rep in range(1, N_REPLICATES + 1):
        logger.debug(f"  Replicate {rep}/{N_REPLICATES}")
        
        # Relax with different random seed
        seed = hash(f"{name}_{rep}") % (2**31)
        relaxed_pos, _ = relax_structure(topology, positions, forcefield, seed, logger)
        
        # Save relaxed structure
        relaxed_dir = output_dir / "relaxed_structures"
        relaxed_dir.mkdir(exist_ok=True)
        relaxed_pdb = relaxed_dir / f"{name}_rep{rep}.pdb"
        with open(relaxed_pdb, 'w') as f:
            PDBFile.writeFile(topology, relaxed_pos, f)
        
        # Calculate MM/GBSA
        energies = calculate_mmgbsa(topology, relaxed_pos, forcefield, logger)
        energies['replicate'] = rep
        energies['name'] = name
        replicate_results.append(energies)
        
        logger.debug(f"    ΔG_bind = {energies['dG_bind']:.2f} kcal/mol")
    
    # Calculate statistics
    dG_values = [r['dG_bind'] for r in replicate_results]
    avg_dG = np.mean(dG_values)
    std_dG = np.std(dG_values)
    
    logger.info(f"  {name}: ΔG_bind = {avg_dG:.2f} ± {std_dG:.2f} kcal/mol")
    
    return {
        'name': name,
        'replicates': replicate_results,
        'dG_avg': avg_dG,
        'dG_std': std_dG
    }


def main():
    """Main execution function."""
    
    # Setup paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent
    output_dir = script_dir
    
    # Input paths
    wt_pdb = project_root / "rosie-51675.flex_input_A_R_refinement/output/model_01.pdb"
    mutant_dir = script_dir / "mutants"
    
    # Setup logging
    logger = setup_logging(output_dir)
    
    # Header
    logger.info("=" * 70)
    logger.info("ALANINE SCANNING ΔΔG CALCULATION")
    logger.info("BIOL 499 Thesis - Dante Maggiotto")
    logger.info("=" * 70)
    logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("")
    
    # Log parameters
    logger.info("PARAMETERS:")
    logger.info(f"  Replicates per structure: {N_REPLICATES}")
    logger.info(f"  Heating steps: {HEATING_STEPS}")
    logger.info(f"  Equilibration steps: {EQUILIBRATION_STEPS}")
    logger.info(f"  Experimental KD: {KD_EXP*1e6:.0f} μM")
    logger.info(f"  Experimental ΔG: {DG_EXP:.2f} kcal/mol")
    logger.info(f"  Hotspot threshold: >{HOTSPOT_THRESHOLD} kcal/mol")
    logger.info(f"  Moderate threshold: {MODERATE_THRESHOLD}-{HOTSPOT_THRESHOLD} kcal/mol")
    logger.info("")
    
    # Load force field
    logger.info("Loading force field (amber99sb + GB-OBC2)...")
    forcefield = ForceField('amber99sb.xml', 'implicit/obc2.xml')
    
    # Collect all structures to process
    structures = [('WT', wt_pdb)]
    for i in range(1, 24):
        mutant_name = f"Ala{i:02d}"
        mutant_path = mutant_dir / f"model_01_Ala{i:02d}.pdb"
        if mutant_path.exists():
            structures.append((mutant_name, mutant_path))
        else:
            logger.warning(f"Missing mutant file: {mutant_path}")
    
    logger.info(f"Found {len(structures)} structures to process")
    logger.info("")
    
    # Process all structures
    logger.info("=" * 70)
    logger.info("PROCESSING STRUCTURES")
    logger.info("=" * 70)
    
    all_results = []
    energy_log = []
    
    for name, pdb_path in structures:
        result = process_structure(name, pdb_path, forcefield, output_dir, logger)
        all_results.append(result)
        energy_log.extend(result['replicates'])
    
    # Save raw energy log
    logger.info("")
    logger.info("Saving raw energy data...")
    energy_csv = output_dir / "energy_log.csv"
    with open(energy_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['name', 'replicate', 'E_complex', 
                                                'E_receptor', 'E_peptide', 'dG_bind'])
        writer.writeheader()
        writer.writerows(energy_log)
    logger.info(f"  Saved: {energy_csv}")
    
    # Calculate scaling factor from WT
    logger.info("")
    logger.info("=" * 70)
    logger.info("CALIBRATION")
    logger.info("=" * 70)
    
    wt_result = next(r for r in all_results if r['name'] == 'WT')
    wt_dG_avg = wt_result['dG_avg']
    
    logger.info(f"  WT ΔG_bind (avg): {wt_dG_avg:.2f} kcal/mol")
    logger.info(f"  Experimental ΔG: {DG_EXP:.2f} kcal/mol")
    
    scale_factor = DG_EXP / wt_dG_avg
    logger.info(f"  Scale factor: {scale_factor:.4f}")
    
    # Apply scaling and calculate ΔΔG
    logger.info("")
    logger.info("=" * 70)
    logger.info("ΔΔG CALCULATION")
    logger.info("=" * 70)
    
    wt_dG_scaled = wt_dG_avg * scale_factor
    
    final_results = []
    for result in all_results:
        name = result['name']
        dG_raw = result['dG_avg']
        dG_scaled = dG_raw * scale_factor
        
        if name == 'WT':
            ddG = 0.0
            position = 0
            wt_res = '-'
        else:
            ddG = dG_scaled - wt_dG_scaled
            position = int(name.replace('Ala', ''))
            wt_res = HSN1_SEQUENCE[position - 1]
        
        # Classify
        if ddG > HOTSPOT_THRESHOLD:
            classification = 'Hotspot'
        elif ddG > MODERATE_THRESHOLD:
            classification = 'Moderate'
        else:
            classification = 'Scaffold'
        
        final_results.append({
            'position': position,
            'WT_residue': wt_res,
            'mutant_residue': 'A',
            'dG_raw': dG_raw,
            'dG_raw_std': result['dG_std'],
            'dG_scaled': dG_scaled,
            'dG_WT_scaled': wt_dG_scaled,
            'ddG': ddG,
            'classification': classification
        })
        
        if name != 'WT':
            logger.info(f"  Position {position:2d} ({wt_res}→A): ΔΔG = {ddG:+.2f} kcal/mol [{classification}]")
    
    # Sort by position
    final_results.sort(key=lambda x: x['position'])
    
    # Save final results CSV
    logger.info("")
    logger.info("Saving results...")
    
    results_csv = output_dir / "alanine_scan_results.csv"
    with open(results_csv, 'w', newline='') as f:
        fieldnames = ['position', 'WT_residue', 'mutant_residue', 'dG_raw', 
                      'dG_raw_std', 'dG_scaled', 'dG_WT_scaled', 'ddG', 'classification']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(final_results)
    logger.info(f"  Saved: {results_csv}")
    
    # Generate heatmap
    logger.info("Generating heatmap...")
    generate_heatmap(final_results, output_dir, logger)
    
    # Summary statistics
    logger.info("")
    logger.info("=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    
    mutant_results = [r for r in final_results if r['position'] > 0]
    hotspots = [r for r in mutant_results if r['classification'] == 'Hotspot']
    moderate = [r for r in mutant_results if r['classification'] == 'Moderate']
    scaffold = [r for r in mutant_results if r['classification'] == 'Scaffold']
    
    logger.info(f"  Total positions analyzed: {len(mutant_results)}")
    logger.info(f"  Hotspots (ΔΔG > {HOTSPOT_THRESHOLD}): {len(hotspots)}")
    if hotspots:
        hotspot_str = ', '.join([f"{r['position']}{r['WT_residue']}" for r in hotspots])
        logger.info(f"    Positions: {hotspot_str}")
    
    logger.info(f"  Moderate ({MODERATE_THRESHOLD}-{HOTSPOT_THRESHOLD}): {len(moderate)}")
    if moderate:
        moderate_str = ', '.join([f"{r['position']}{r['WT_residue']}" for r in moderate])
        logger.info(f"    Positions: {moderate_str}")
    
    logger.info(f"  Scaffold (< {MODERATE_THRESHOLD}): {len(scaffold)}")
    
    logger.info("")
    logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 70)


def generate_heatmap(results: List[Dict], output_dir: Path, logger: logging.Logger):
    """
    Generate heatmap visualization for Figure 1.
    
    Creates a bar plot showing ΔΔG for each position, colored by classification.
    """
    # Filter out WT (position 0)
    mutant_results = [r for r in results if r['position'] > 0]
    
    positions = [r['position'] for r in mutant_results]
    ddG_values = [r['ddG'] for r in mutant_results]
    labels = [f"{r['position']}{r['WT_residue']}" for r in mutant_results]
    
    # Colors by classification
    colors = []
    for r in mutant_results:
        if r['classification'] == 'Hotspot':
            colors.append('#d62728')  # Red
        elif r['classification'] == 'Moderate':
            colors.append('#ff7f0e')  # Orange
        else:
            colors.append('#2ca02c')  # Green
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 6))
    
    bars = ax.bar(positions, ddG_values, color=colors, edgecolor='black', linewidth=0.5)
    
    # Add threshold lines
    ax.axhline(y=HOTSPOT_THRESHOLD, color='#d62728', linestyle='--', linewidth=1.5, 
               label=f'Hotspot threshold ({HOTSPOT_THRESHOLD} kcal/mol)')
    ax.axhline(y=MODERATE_THRESHOLD, color='#ff7f0e', linestyle='--', linewidth=1.5,
               label=f'Moderate threshold ({MODERATE_THRESHOLD} kcal/mol)')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    # Labels
    ax.set_xlabel('H-SN1 Position', fontsize=12, fontweight='bold')
    ax.set_ylabel('ΔΔG (kcal/mol)', fontsize=12, fontweight='bold')
    ax.set_title('Alanine Scanning of H-SN1 Binding to TNFR1\nIdentification of Binding Hotspots',
                 fontsize=14, fontweight='bold')
    
    # X-axis
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#d62728', edgecolor='black', label='Hotspot (>1.5 kcal/mol)'),
        mpatches.Patch(facecolor='#ff7f0e', edgecolor='black', label='Moderate (0.8-1.5 kcal/mol)'),
        mpatches.Patch(facecolor='#2ca02c', edgecolor='black', label='Scaffold (<0.8 kcal/mol)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    # Grid
    ax.yaxis.grid(True, linestyle=':', alpha=0.7)
    ax.set_axisbelow(True)
    
    # Tight layout
    plt.tight_layout()
    
    # Save
    heatmap_path = output_dir / "alanine_scan_heatmap.png"
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  Saved: {heatmap_path}")


if __name__ == "__main__":
    main()

