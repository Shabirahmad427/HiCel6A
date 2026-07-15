#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
src_root="$(cd "${repo_root}/.." && pwd)"
dest_root="${repo_root}/inputs"
manifest="${repo_root}/metadata/file_manifest.txt"

mkdir -p "${dest_root}" "$(dirname "${manifest}")"

copy_matches() {
    while IFS= read -r relpath; do
        [ -n "${relpath}" ] || continue
        mkdir -p "${dest_root}/$(dirname "${relpath}")"
        cp -p "${src_root}/${relpath}" "${dest_root}/${relpath}"
    done
}

find \
    1OC5/00_Apo \
    1OCB/00_Apo \
    1OCJ/00_Apo \
    -type f \
    \( \
        -name '*.pqr' -o \
        -path '*/amber/README' -o \
        -name 'step3_input.*' -o \
        -name 'step4.0_minimization.mdin' -o \
        -name 'step4.1_equilibration.mdin' -o \
        -name 'step5_production.mdin' -o \
        -name 'sysinfo.dat' \
    \) | sort | copy_matches

find \
    1OC5/01_Holo/00_preparation \
    1OCB/01_Holo/00_preparation \
    1OCJ/01_Holo/00_preparation \
    -type f \
    \( \
        -name 'README' -o \
        -name 'README_BGLC8_FIX.md' -o \
        -name 'step3_input.*' -o \
        -name 'step4.0_minimization.mdin' -o \
        -name 'step4.1_equilibration.mdin' -o \
        -name 'step5_production.mdin' -o \
        -name 'sysinfo.dat' -o \
        -name 'dihe.restraint' \
    \) | sort | copy_matches

find \
    1OC5/01_Holo/rep1 1OC5/01_Holo/rep2 1OC5/01_Holo/rep3 \
    1OCB/01_Holo/rep1 1OCB/01_Holo/rep2 1OCB/01_Holo/rep3 \
    1OCJ/01_Holo/rep1 1OCJ/01_Holo/rep2 1OCJ/01_Holo/rep3 \
    -type f \
    \( \
        -name 'README' -o \
        -name 'step3_input.*' -o \
        -name 'step4.0_minimization.mdin' -o \
        -name 'step4.1_equilibration.mdin' -o \
        -name 'step5_production.mdin' -o \
        -name 'sysinfo.dat' -o \
        -name 'dihe.restraint' -o \
        -name 'plumed.in' -o \
        -name 'mmpbsa.in' -o \
        -name 'cpptraj*.in' \
    \) | sort | copy_matches

(
    cd "${repo_root}"
    find inputs -type f | sort > "${manifest}"
)

printf 'Collected %s files into %s\n' "$(wc -l < "${manifest}")" "${dest_root}"
