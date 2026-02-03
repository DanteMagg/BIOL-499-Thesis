# RMSD Validation Script
# Compares FlexPepDock model_01 to 1TNR crystal structure
# Validates binding at CRD2/3 interface

reinitialize

# Load structures
load rosie-51675.flex_input_A_R_refinement/output/model_01.pdb, flexpepdock
load data/structures/1TNR_clean.pdb, zheng

# Align peptide+CRD2/3 region
align flexpepdock and (chain A or chain R), zheng and resi 80-140

# Measure RMSD (backbone CRD2/3)
# Note: zheng only has chain R, so we compare receptor regions
rms_cur flexpepdock and chain R and resi 80-140 and name N+CA+C, zheng and chain R and resi 80-140 and name N+CA+C

# Print results
print "="*60
print "RMSD Validation: FlexPepDock vs 1TNR Crystal Structure"
print "="*60

# Visualize
hide everything

# Show FlexPepDock model
show cartoon, flexpepdock and chain A
color cyan, flexpepdock and chain A
show surface, flexpepdock and chain R
color wheat, flexpepdock and chain R
set transparency, 0.3, flexpepdock and chain R

# Show reference structure
show cartoon, zheng and chain R
color gray70, zheng and chain R

# Highlight CRD2/3 region
select crd23, resi 80-140
show sticks, crd23 and flexpepdock
util.cbag crd23 and flexpepdock

# Labels
label flexpepdock and chain A and name CA and resi 1, "H-SN1 N-term"

# Set view
bg_color white
set ray_shadows, 0
zoom flexpepdock

# Save image
ray 1200, 800
png rmsd_validation_flexpepdock_vs_1TNR.png, dpi=300

print "="*60
print "Image saved: rmsd_validation_flexpepdock_vs_1TNR.png"
print "="*60

