# HiCel6A

Input files and analysis scripts for the molecular dynamics study

**Loop Dynamics, Substrate Binding, and Active Site Hydration in HiCel6A and its Catalytically Impaired D416A and D405N Mutants**

Shabir Ahmad, Clauber H. S. da Costa, Ander F. Pereira, Alberto M. dos Santos, Munir S. Skaf
Institute of Chemistry, University of Campinas (UNICAMP), Brazil

## Systems

| name in this repository | variant | PDB entry | substrate |
|---|---|---|---|
| `HiCel6A` | wild type | 1OCB | cellooctaose (C8) |
| `HiCel6AD416A` | D416A | 1OCJ | cellooctaose (C8) |
| `HiCel6AD405N` | D405N | 1OC5 | cellooctaose (C8) |

Each system was simulated in three independent replicas of 2 µs (`rep1`, `rep2`, `rep3`).

## Repository layout

```
HiCel6A/
├── inputs/                      simulation and analysis inputs, one folder per system
│   └── <system>/01_Holo/
│       ├── 00_preparation/      topology, coordinates and production input (Amber, CHARMM-GUI)
│       └── rep{1,2,3}/          per replica:
│           ├── 04_traj/         trajectory processing (cpptraj)
│           ├── 05_RMSD/         RMSD, RMSF and 2D RMSD
│           ├── 06_distance/     catalytic-residue to substrate distances
│           ├── 07_PCA/          principal components and free-energy landscapes
│           ├── 08_pcasuite/     collective-motion projections
│           ├── 09_MMGBSA/       MM-GBSA binding free energy, 200 ns windows
│           └── 10_extra/        MDDF hydration analysis and correlation inputs
├── scripts/
│   ├── collect_inputs.sh        rebuilds `inputs/` from the simulation folders
│   └── water_C1/                water near the anomeric carbon (Figure S12), see its README
└── metadata/                    file manifest and release checklist
```

Trajectories are not included because of their size (see Data availability).

## Analyses in this repository

| analysis | where |
|---|---|
| RMSD, RMSF, 2D RMSD | `inputs/<system>/01_Holo/rep*/05_RMSD/` |
| Substrate–residue distances | `inputs/<system>/01_Holo/rep*/06_distance/` |
| PCA, free-energy landscapes, collective motions | `inputs/<system>/01_Holo/rep*/07_PCA/`, `08_pcasuite/` |
| MM-GBSA binding free energy | `inputs/<system>/01_Holo/rep*/09_MMGBSA/` |
| MDDF hydration analysis | `inputs/<system>/01_Holo/rep*/10_extra/mddf/` |
| Catalytic-position water and water wire (Figure S12) | `scripts/water_C1/` |

## Data availability

- Input files and scripts: this repository.
- Archived release: Zenodo, DOI to be added at release.
- Trajectories: available from the corresponding author on request, or deposited separately and linked here.

## License

All files in this repository are released under the Creative Commons Attribution 4.0 International license (CC BY 4.0), see `LICENSE`.

## Citation

Please cite the article and the Zenodo record (see `CITATION.cff`).
