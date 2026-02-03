#!/usr/bin/env python3
"""
Prepare model_1.pdb for FlexPepDock submission

This script takes the full-atom model_1.pdb from CABS-dock and prepares it
for FlexPepDock refinement by:
1. Adding proper header information
2. Preserving secondary structure records
3. Renumbering atoms sequentially
4. Ensuring proper chain order (R before A)
"""

# Update paths as needed - script assumes CABSdock directory is in same location
input_file = "CABSdock_ea91db3082d4878/model_1.pdb"
output_file = "H-SN1_TNFR1_flexpepdock.pdb"

# Read and filter the PDB
with open(input_file, 'r') as f:
    lines = f.readlines()

with open(output_file, 'w') as f:
    # Write clean header
    f.write("HEADER    PEPTIDE-RECEPTOR COMPLEX                   22-JAN-26   XXXX\n")
    f.write("TITLE     H-SN1 PEPTIDE BOUND TO TNFR1 FOR FLEXPEPDOCK REFINEMENT\n")
    f.write("COMPND    MOL_ID: 1;\n")
    f.write("COMPND   2 MOLECULE: TUMOR NECROSIS FACTOR RECEPTOR 1;\n")
    f.write("COMPND   3 CHAIN: R;\n")
    f.write("COMPND   4 MOL_ID: 2;\n")
    f.write("COMPND   5 MOLECULE: H-SN1 ANTI-INFLAMMATORY PEPTIDE;\n")
    f.write("COMPND   6 CHAIN: A\n")
    f.write("REMARK   1 CABS-DOCK MODEL 1 - PREPARED FOR FLEXPEPDOCK\n")
    f.write("REMARK   2 CHAIN R: TNFR1 (RESIDUES 15-153)\n")
    f.write("REMARK   3 CHAIN A: H-SN1 PEPTIDE (23 AA)\n")
    
    # Write HELIX/SHEET records from original
    for line in lines:
        if line.startswith('HELIX') or line.startswith('SHEET'):
            f.write(line)
    
    # Renumber atoms and write ATOM records
    atom_num = 1
    for line in lines:
        if line.startswith('ATOM'):
            # Renumber atoms
            new_line = f"ATOM  {atom_num:5d}{line[11:]}"
            f.write(new_line)
            atom_num += 1
    
    f.write("TER\n")
    f.write("END\n")

print(f"Created: {output_file}")
print(f"Total atoms: {atom_num - 1}")

# Verify chains
with open(output_file, 'r') as f:
    lines = f.readlines()

chain_a = sum(1 for l in lines if l.startswith('ATOM') and l[21] == 'A')
chain_r = sum(1 for l in lines if l.startswith('ATOM') and l[21] == 'R')
print(f"Chain A (peptide): {chain_a} atoms")
print(f"Chain R (receptor): {chain_r} atoms")

