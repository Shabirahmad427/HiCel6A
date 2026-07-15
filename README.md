# HiCel6a Reproducibility Inputs

This directory is a GitHub-ready archive for the HiCel6a simulation input files.
It is intended to become the contents of a private GitHub repository first, and
later a public release connected to Zenodo for journal submission.

## Scope

Included here:

- AMBER topology and coordinate inputs (`step3_input.*`, `parm7`, `rst7`, `pdb`, `psf`)
- MD control files (`step4.0_minimization.mdin`, `step4.1_equilibration.mdin`, `step5_production.mdin`)
- System description files (`sysinfo.dat`, `dihe.restraint`, selected `README` files)
- Metadynamics inputs (`plumed.in`)
- Analysis input scripts commonly requested for reproducibility (`cpptraj*.in`, `mmpbsa.in`)
- Original preparation notes for the BGLC8 carbohydrate fix

Not included here:

- Long production trajectories
- Restart/output logs from completed runs
- Large derived analysis outputs and figures

Those larger files are usually deposited separately in Zenodo, OSF, or an
institutional repository if the journal requires raw trajectory access.

## Structure

- `inputs/`: copied input files preserving their original relative paths
- `metadata/`: archive manifest and submission checklist
- `scripts/collect_inputs.sh`: rebuilds `inputs/` from the main HiCel6a workspace

## Recommended GitHub -> Zenodo workflow

1. Create a new **private** GitHub repository.
2. Put the contents of this directory into that repository.
3. Review the archive and remove anything you do not want public later.
4. When the manuscript is ready, make the repository public or publish a tagged release.
5. Connect the GitHub repository to Zenodo.
6. Create a GitHub release such as `v1.0.0`.
7. Zenodo will mint a DOI for that release.
8. Cite the Zenodo DOI in the paper.

## What people usually upload for journal submission

For MD input deposits, the common minimum is:

- starting structures
- topology/parameter files
- all run control files
- exact analysis scripts used to make the paper figures
- a short README explaining directory meaning and software versions

For stronger reproducibility, people also add:

- a manifest of files
- software/version notes
- one DOI-backed frozen release in Zenodo
- separate large-data archive links for trajectories if GitHub is too small

## Before public release

Check the following carefully:

- no personal paths that should be hidden
- no unnecessary large binaries
- no confidential reviewer notes
- no temporary scripts that would confuse readers
- README text matches the final manuscript terminology

## Notes

This archive was assembled from:

- `1OC5`
- `1OCB`
- `1OCJ`

If you want, the next step can be either:

- initialize this folder as a local git repository, or
- refine the archive scope before pushing to GitHub
