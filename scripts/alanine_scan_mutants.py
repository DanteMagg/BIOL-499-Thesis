#!/usr/bin/env python3
"""
Alanine Scanning Mutant Generator
=================================
Creates 23 single-point alanine mutants of chain A (H-SN1).
Chain R (TNFR1) is left untouched.

Input:  model_01.pdb (FlexPepDock refined complex)
Output: model_01_Ala01.pdb ... model_01_Ala23.pdb
"""

from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.Data.IUPACData import protein_letters_3to1
import copy
import warnings

# Suppress Biopython PDB warnings
warnings.filterwarnings('ignore')


class ChainSelector(Select):
    """Select specific chains for output."""
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids


def mutate_to_ala(residue):
    """
    Mutate a residue to alanine by:
    1. Keeping backbone atoms (N, CA, C, O, H, HA)
    2. Replacing sidechain with CB only
    3. Updating residue name to ALA
    """
    backbone_atoms = {'N', 'CA', 'C', 'O', 'H', 'HA', 'H1', 'H2', 'H3', '1H', '2H', '3H'}
    
    # Get atoms to remove (sidechain except CB)
    atoms_to_remove = []
    has_cb = False
    
    for atom in residue.get_atoms():
        atom_name = atom.get_name().strip()
        if atom_name == 'CB':
            has_cb = True
        elif atom_name not in backbone_atoms:
            atoms_to_remove.append(atom.get_id())
    
    # Remove sidechain atoms (except CB)
    for atom_id in atoms_to_remove:
        residue.detach_child(atom_id)
    
    # Update residue name to ALA
    residue.resname = 'ALA'
    
    return residue


def generate_ala_mutants(input_pdb: Path, output_dir: Path):
    """Generate 23 alanine mutants for chain A positions 1-23."""
    
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    
    # Load structure
    print(f"Loading: {input_pdb}")
    structure = parser.get_structure('complex', input_pdb)
    
    # Get chain A info
    model = structure[0]
    chain_a = model['A']
    residues = list(chain_a.get_residues())
    
    print(f"Chain A: {len(residues)} residues")
    print(f"Chain R: {len(list(model['R'].get_residues()))} residues (unchanged)")
    print("-" * 50)
    
    # Generate mutants
    mutant_files = []
    
    for i, res in enumerate(residues, start=1):
        res_name = res.get_resname()
        res_id = res.get_id()[1]
        
        # Skip if already ALA or GLY (no sidechain to remove)
        if res_name == 'ALA':
            print(f"  Position {i:2d}: {res_name} → ALA (already ALA, copying as-is)")
        elif res_name == 'GLY':
            print(f"  Position {i:2d}: {res_name} → ALA (GLY→ALA adds CB)")
        else:
            print(f"  Position {i:2d}: {res_name} → ALA")
        
        # Deep copy the structure
        mutant_structure = copy.deepcopy(structure)
        mutant_model = mutant_structure[0]
        mutant_chain_a = mutant_model['A']
        mutant_residues = list(mutant_chain_a.get_residues())
        
        # Mutate the target residue
        target_res = mutant_residues[i - 1]
        
        if res_name != 'ALA':
            mutate_to_ala(target_res)
        
        # Write output
        output_file = output_dir / f"model_01_Ala{i:02d}.pdb"
        io.set_structure(mutant_structure)
        io.save(str(output_file))
        mutant_files.append(output_file)
    
    return mutant_files


def main():
    # Paths
    project = Path(__file__).parent.parent
    input_pdb = project / "rosie-51675.flex_input_A_R_refinement/output/model_01.pdb"
    output_dir = project / "results/alanine_scanning/mutants"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 50)
    print("ALANINE SCANNING MUTANT GENERATOR")
    print("=" * 50)
    
    # Generate mutants
    mutant_files = generate_ala_mutants(input_pdb, output_dir)
    
    print("-" * 50)
    print(f"\nGenerated {len(mutant_files)} mutant PDBs in:")
    print(f"  {output_dir}/")
    print("\nFiles:")
    for f in mutant_files:
        print(f"  {f.name}")
    
    # Write summary
    summary_file = output_dir.parent / "mutant_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("Alanine Scanning Mutant Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Input: {input_pdb.name}\n")
        f.write(f"Mutants generated: {len(mutant_files)}\n\n")
        f.write("Position → ALA mutations:\n")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('complex', input_pdb)
        chain_a = structure[0]['A']
        
        for i, res in enumerate(chain_a.get_residues(), start=1):
            res_name = res.get_resname()
            one_letter = protein_letters_3to1.get(res_name.capitalize(), 'X')
            f.write(f"  {i:2d}. {res_name} ({one_letter}) → ALA (A)\n")
    
    print(f"\nSummary written to: {summary_file}")
    print("=" * 50)


if __name__ == "__main__":
    main()

