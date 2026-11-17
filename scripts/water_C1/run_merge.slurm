#!/bin/bash
#SBATCH -p par48-i
#SBATCH -N 1
#SBATCH -n 3
#SBATCH -t 04:00:00
#SBATCH -J wC1s10_merge
#SBATCH -o /path/to/water_C1_occupancy/logs/%x_%j.out
# Step 1: merge rep1+rep2+rep3 (stride 10) into one 30 000-frame trajectory per system; no analysis here.
source /path/to/activate_amber.sh
W=/path/to/water_C1_occupancy
PY=python3
one() {
  s=$1
  $W/00_scripts/make_cpptraj_merged.sh $s
  cpptraj -i $W/01_cpptraj/${s}_merged_site.in > $W/logs/${s}_merged_site.log 2>&1
  gzip -f $W/02_reduced/${s}_merged_closest.dat
  cpptraj -i $W/01_cpptraj/${s}_merged_snap.in > $W/logs/${s}_merged_snap.log 2>&1
  $PY - <<PYEOF
import netCDF4
for k, n in (("site", 30000), ("snap", 3000)):
    f = netCDF4.Dataset("$W/02_reduced/${s}_merged_%s.nc" % k); nf = f.variables["coordinates"].shape[0]
    print("$s", k, nf, "frames", "OK" if nf == n else "WRONG (expected %d)" % n, flush=True)
PYEOF
  echo "$s merged $(date)"
}
for s in HiCel6A HiCel6AD405N HiCel6AD416A; do one $s & done
wait
