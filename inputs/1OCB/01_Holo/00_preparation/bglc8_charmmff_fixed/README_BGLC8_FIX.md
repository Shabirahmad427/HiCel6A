# BGLC8 CHARMM36 carbohydrate repair

This folder keeps the original CHARMM-GUI water, ions, box, atom count, and
coordinates. Only the 171 substrate atoms were renamed/retyped and reordered
internally so each BGLC residue is contiguous for CHAMBER.

Changes:
- split the original one-residue CGenFF `BGL`/`BGLC` ligand into 8 `BGLC`
  residues in segment `SUB`;
- remapped the CGenFF atom labels to CHARMM36 carbohydrate atom labels;
- assigned CHARMM36 carbohydrate atom types and charges;
- applied the beta-1,4 linkage charge/type changes corresponding to `14bb`.

Generated files are in `amber/`.
