#!/usr/bin/env python
"""Step 1 of the 'water wire to Asp180' figure: shortlist frames of the merged run in which a water wire (<= 2 waters) links Asp180 to a water
within 4 A of C1 AND a catalytic-position water exists, then extract those frames (ALL waters kept) from the original 2-us trajectories.
merged frame f (0-based) <-> replica f // 10000 + 1, original frame (f % 10000) * 10 + 1 (1-based, cpptraj)."""
import os, subprocess, numpy as np
W = '/path/to/water_C1_occupancy'; M = '/path/to/simulations'
LIG = {'HiCel6A': 362, 'HiCel6AD405N': 365, 'HiCel6AD416A': 361}; NSHORT = 40; MINSEP = 100          # frames (= 20 ns) between shortlisted frames
os.makedirs(f'{W}/16_wire', exist_ok=True)
for s in LIG:
    D = np.load(f'{W}/04_results/per_frame_{s}_merged.npz'); cd, ca = D['cat_d'], D['cat_ang']
    cat = ~np.isnan(cd); wire = D['wire2_Asp180'] & cat
    md, ma = np.nanmedian(cd), np.nanmedian(ca)
    print(f'{s}: frames {len(cd)}, catalytic water present {cat.sum()} ({100 * cat.mean():.1f}%), + wire to Asp180 (<=2 waters) {wire.sum()} ({100 * wire.mean():.2f}%), '
          f'by replica {[int(wire[i * 10000:(i + 1) * 10000].sum()) for i in range(3)]}')
    idx = np.where(wire)[0]; sc = np.abs(cd[idx] - md) / 0.2 + np.abs(ca[idx] - ma) / 10
    pick = []
    for j in np.argsort(sc):
        if all(abs(idx[j] - p) >= MINSEP for p in pick): pick.append(int(idx[j]))
        if len(pick) == NSHORT: break
    pick.sort(); np.save(f'{W}/16_wire/{s}_cand_frames.npy', np.array(pick))
    SITE = open(f'{W}/01_cpptraj/{s}_site_residues.txt').read().strip(); L = LIG[s]
    tr = ''.join(f'trajin {M}/{s}/rep{f // 10000 + 1}/2us_wat.nc {(f % 10000) * 10 + 1} {(f % 10000) * 10 + 1}\n' for f in pick)
    open(f'{W}/01_cpptraj/{s}_wire_cand.in', 'w').write(f"""parm {M}/{s}/step3_input.parm7
reference {M}/{s}/step3_input.rst7 [ref]
{tr}autoimage anchor :1-{L}
rms fit @CA ref [ref]
strip !(:{SITE},{L}|:WAT) outprefix {W}/16_wire/{s}_cand
trajout {W}/16_wire/{s}_cand.nc netcdf
run
""")
    subprocess.run(['cpptraj', '-i', f'{W}/01_cpptraj/{s}_wire_cand.in'], check=True, capture_output=True)
    print(f'   shortlisted {len(pick)} frames, e.g. merged frames {pick[:6]} -> {W}/16_wire/{s}_cand.nc')
