#!/bin/bash
# Merge rep1+rep2+rep3 of one system into ONE trajectory (stride 10 = 200 ps -> 3 x 10 000 = 30 000 frames, order rep1,rep2,rep3).
# usage: make_cpptraj_merged.sh SYSTEM
#   <sys>_merged_site.in : stride-10 frames of the three replicas, autoimage, Ca fit to step3_input, pooled water-O density grid
#                          (centred on C30, 0.5 A, relative to bulk), 80 waters closest to C30 per frame, protein stripped to residues
#                          within 12 A of C30 -> 02_reduced/<sys>_merged_site.nc (30 000 frames) + <sys>_merged_closest.dat
#   <sys>_merged_snap.in : every 100th ORIGINAL frame (2 ns) of the three replicas, whole protein + 80 closest waters
#                          -> 02_reduced/<sys>_merged_snap.nc (3 x 1000 = 3000 snapshots, order rep1,rep2,rep3)
set -e
s=$1
B=/path/to/data M=$B/manuscrito; W=$B/analyses/water_C1_occupancy_s10
case $s in HiCel6A) L=362;; HiCel6AD405N) L=365;; HiCel6AD416A) L=361;; *) echo bad system; exit 1;; esac
# grid centre = C30 position in step3_input.rst7 (the frame all trajectories are fitted to)
printf "parm $M/$s/step3_input.parm7\ntrajin $M/$s/step3_input.rst7\nvector c30 center :$L@C30 out $W/01_cpptraj/${s}_C30_rst7.dat\nrun\n" > $W/01_cpptraj/${s}_C30_rst7.in
cpptraj -i $W/01_cpptraj/${s}_C30_rst7.in > /dev/null
read x y z < <(tail -1 $W/01_cpptraj/${s}_C30_rst7.dat | awk "{print \$2,\$3,\$4}")
SITE=$(cat $W/01_cpptraj/${s}_site_residues.txt)
TR10=""; TR100=""
for r in rep1 rep2 rep3; do TR10="${TR10}trajin $M/$s/$r/2us_wat.nc 1 last 10\n"; TR100="${TR100}trajin $M/$s/$r/2us_wat.nc 1 last 100\n"; done
printf "parm $M/$s/step3_input.parm7\nreference $M/$s/step3_input.rst7 [ref]\n${TR10}autoimage anchor :1-$L\nrms fit @CA ref [ref]\ngrid $W/03_grids/${s}_merged_waterO.dx 64 0.5 64 0.5 64 0.5 gridcenter $x $y $z :WAT@O normframe\nclosest 80 :$L@C30 oxygen closestout $W/02_reduced/${s}_merged_closest.dat\nstrip !(:$SITE,$L|:WAT) outprefix $W/02_reduced/${s}_site\ntrajout $W/02_reduced/${s}_merged_site.nc netcdf\nrun\n" > $W/01_cpptraj/${s}_merged_site.in
printf "parm $M/$s/step3_input.parm7\nreference $M/$s/step3_input.rst7 [ref]\n${TR100}autoimage anchor :1-$L\nrms fit @CA ref [ref]\nclosest 80 :$L@C30 oxygen outprefix $W/02_reduced/${s}_full\ntrajout $W/02_reduced/${s}_merged_snap.nc netcdf\nrun\n" > $W/01_cpptraj/${s}_merged_snap.in
