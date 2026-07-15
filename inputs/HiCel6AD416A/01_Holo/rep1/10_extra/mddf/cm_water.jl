using PDBTools, ComplexMixtures, DelimitedFiles

# Load PDB file of the system
atoms = readPDB("../../00_preparation/step3_input.pdb")
protein = select(atoms,"protein")
water = select(atoms,"resname WAT")

solute = ComplexMixtures.AtomSelection(protein,nmols=1)
solvent = ComplexMixtures.AtomSelection(water,natomspermol=3)

#READING TRAJECTORY
trajectory = ComplexMixtures.Trajectory("../../04_traj/2us_wat.nc",solute,solvent)

#GETTING MDDFs
options = ComplexMixtures.Options(bulk_range=(8.0, 10.0), firstframe=1, lastframe=5000)
results = ComplexMixtures.mddf(trajectory,options)

#SAVING RESULTS
ComplexMixtures.save(results,"./results-water.json")
