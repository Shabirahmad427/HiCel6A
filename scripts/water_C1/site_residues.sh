#!/bin/bash
# Residues (protein) with any atom within 12 A of ligand C30 in step3_input -> 01_cpptraj/<sys>_site_residues.txt
B=/path/to/data M=$B/manuscrito; W=$B/analyses/water_C1_occupancy_s10
for s in HiCel6A:362 HiCel6AD405N:365 HiCel6AD416A:361; do n=${s%%:*}; L=${s##*:}
 printf "parm $M/$n/step3_input.parm7\nreference $M/$n/step3_input.rst7\nselect (:$L@C30<:12.0)&!(:WAT,SOD,CLA,$L)\n" > $W/01_cpptraj/${n}_site_select.in
 cpptraj -i $W/01_cpptraj/${n}_site_select.in | grep "^Selected=" | tr ' ' '\n' | grep -E '^[0-9]+$' > $W/01_cpptraj/${n}_site_atoms.tmp
 printf "parm $M/$n/step3_input.parm7\natominfo @$(paste -sd, $W/01_cpptraj/${n}_site_atoms.tmp)\n" > $W/01_cpptraj/${n}_tmp.in
 cpptraj -i $W/01_cpptraj/${n}_tmp.in | awk '$1~/^[0-9]+$/ && NF>8{print $3}' | sort -nu | paste -sd, > $W/01_cpptraj/${n}_site_residues.txt
 rm $W/01_cpptraj/${n}_site_atoms.tmp $W/01_cpptraj/${n}_tmp.in
 echo "$n: $(tr ',' '\n' < $W/01_cpptraj/${n}_site_residues.txt | wc -l) residues"
done
