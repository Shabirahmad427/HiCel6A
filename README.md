# HiCel6A Reproducibility Inputs

This directory is a GitHub-ready archive for the HiCel6A simulation input files.
It is intended to become the contents of a private GitHub repository first, and
later a public release connected to Zenodo for journal submission.

## Scope

Included here:

- AMBER topology and coordinate inputs (`step3_input.*`, `parm7`, `rst7`, `pdb`, `psf`)
- Production control files and related run inputs (`step5_production.mdin`, topology, coordinates)
- System description files (`sysinfo.dat`, `dihe.restraint`, selected `README` files)
- Analysis inputs and key reproducibility files for RMSD, PCA/PCASuite, correlation, MMGBSA, and MDDF

Not included here:

- Long production trajectories
- Restart/output logs from completed runs
- Large derived analysis outputs and figures

Those larger files are usually deposited separately in Zenodo, OSF, or an
institutional repository if the journal requires raw trajectory access.

## Structure

- `inputs/`: curated input files organized under paper-facing system names
- `metadata/`: archive manifest and submission checklist
- `scripts/collect_inputs.sh`: rebuilds `inputs/` from the main HiCel6A workspace

Current top-level system folders in `inputs/`:

- `HiCel6A`
- `HiCel6AD405N`
- `HiCel6AD416A`

Current curation notes:

- only holo systems used in the paper are included
- system names in `inputs/` follow paper-facing labels, not raw source folder names
- `HiCel6A` includes a curated `01_Holo/00_preparation` directory
- `Filesph7`, `bglc8_charmmff_fixed`, metadynamics, apo systems, and `input_run*` preparation subfolders are excluded

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
- the run control files actually used for the reported production work
- exact analysis scripts and key matrix/result files used to make the paper figures
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

This archive was assembled from these paper-facing system labels:

- `HiCel6A` from source folder `1OCB`
- `HiCel6AD405N` from source folder `1OC5`
- `HiCel6AD416A` from source folder `1OCJ`

Only the holo systems used in the paper are included. Apo-only materials were
excluded from this archive.

If you want, the next step can be either:

- review the current archive before pushing to GitHub, or
- prepare a clean Git commit matching the final curated scope
