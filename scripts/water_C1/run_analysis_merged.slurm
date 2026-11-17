#!/bin/bash
#SBATCH -p par48-i
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -t 12:00:00
#SBATCH -J wC1s10_ana
#SBATCH -o /path/to/water_C1_occupancy/logs/%x_%j.out
# Step 2: C1 analysis (3-4 A shell, and 6 A shell) on the MERGED 30 000-frame trajectory of each system.
# WC1_SEG = frames per replica: residence episodes / catalytic-water runs are cut at the rep1|rep2 and rep2|rep3 junctions.
source /path/to/activate_amber.sh
W=/path/to/water_C1_occupancy
PY=python3
export WC1_SEG=10000
# stop unless every merged trajectory has exactly 3 x 10 000 frames
$PY - <<PYEOF || exit 1
import sys, netCDF4
bad = 0
for s in ("HiCel6A", "HiCel6AD405N", "HiCel6AD416A"):
    n = netCDF4.Dataset("$W/02_reduced/%s_merged_site.nc" % s).variables["coordinates"].shape[0]
    print(s, n, "frames", flush=True); bad += n != 30000
sys.exit(bad)
PYEOF
for s in HiCel6A HiCel6AD405N HiCel6AD416A; do
  $PY $W/00_scripts/analyse_water_C1.py    $s merged > $W/logs/${s}_merged_analyse.log    2>&1 && echo "$s 3-4A done $(date)" &
done
wait
