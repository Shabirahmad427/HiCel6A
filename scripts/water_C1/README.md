# Water in the catalytic position near the anomeric carbon (C1)

Scripts that produce **Figure S12** and the numbers of the corresponding Results paragraph
(HiCel6A, HiCel6A-D416A and HiCel6A-D405N, each with the cellooctaose substrate).

## Data analysed

- Three independent 2 µs replicas per system, read every 10th frame (200 ps) and concatenated: **30 000 frames per system**.
- Only the 80 waters closest to C1 in every frame were kept for the analysis.
- Trajectories are not included in this repository.

## Definitions

| item | definition |
|---|---|
| C1 | C30, anomeric carbon of the −1 glucose |
| O_glyc | scissile glycosidic oxygen (O19) |
| Hydrogen bond | donor–acceptor distance ≤ 3.5 Å and D–H···A angle ≥ 135° |
| Catalytic-position water | closest water oxygen within 4.0 Å of C1 with the angle Ow···C1–O_glyc ≥ 135° |
| Attack-ready geometry | catalytic-position water with d(Ow–C1) ≤ 3.5 Å and angle ≥ 150° |
| Water wire (≤ 2 waters) | H-bonds from a water within 4.0 Å of C1, through at most one further water, to the carboxylate of Asp180 |
| Strict wire | the same, with the wire required to start at the catalytic-position water |

## Scripts, in run order

| step | files | what it does |
|---|---|---|
| 1 | `run_merge.slurm`, `make_cpptraj_merged.sh`, `site_residues.sh` | merge the three replicas (stride 10) and write the reduced trajectories |
| 2 | `run_analysis_merged.slurm`, `analyse_water_C1.py` | per-frame occupancy, catalytic-position water (distance, angle), H-bond partners and water wires |
| 3 | `wire_candidates.py`, `wire_eval.py` | strict water wire from the catalytic-position water to Asp180 |
| 4 | `boxplots_system.py` | Figure S12 (one box per system) and the data of each panel |

Figure S12: whiskers are the 5th and 95th percentiles. Panels A and D give the percentage of frames in consecutive blocks of 50 frames
(10 ns), 600 values per system. Panels B and C use only the frames with a catalytic-position water
(810, 2231 and 2488 frames for HiCel6A, D416A and D405N).

## Requirements

Python 3 (NumPy, pandas, netCDF4, matplotlib), Amber/cpptraj, SLURM. Paths inside the scripts point to the cluster where the analysis was run and must be adapted.
