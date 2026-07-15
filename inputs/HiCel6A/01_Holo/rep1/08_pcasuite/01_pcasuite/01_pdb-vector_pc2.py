# Rename object to remove dots
set_name protein.traj.proj2.renumbered, protein_traj_proj2_renumbered

# Split the desired states
split_states protein_traj_proj2_renumbered, 6, 6
split_states protein_traj_proj2_renumbered, 16, 16

# Rename the newly created objects
set_name protein_traj_proj2_renumbered_0006, protein_traj_proj2_0006
set_name protein_traj_proj2_renumbered_0016, protein_traj_proj2_0016

# Load modevector script
run modevector.py

# Show protein as cartoon
as cartoon, protein_traj_proj2_0006 or protein_traj_proj2_0016

# Define protein color from hex
set_color protein_col, [8, 104, 172]

# Apply color to both states
color protein_col, protein_traj_proj2_0006
color protein_col, protein_traj_proj2_0016

# Generate vectors (single call)
modevectors protein_traj_proj2_0006, protein_traj_proj2_0016, \
    headrgb=(1.0,0.0,0.0), tailrgb=(0.0,1.0,0.0)

# Show vectors
show cgo, modevectors

# Render
set ray_opaque_background, off
ray 2400,2400
png proj2.pymol, dpi=600


#This controls the colour of the cone/arrow head using RGB values (Fig.6)
#modevectors 6, 16, headrgb=(1.0,0.0,0.0)
#This controls the colour of the cylinder/arrow tails using RGB values (Fig.7)
#modevectors 6, 16, tailrgb=(1.0,0.0,0.0)
