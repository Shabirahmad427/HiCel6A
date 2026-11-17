#!/usr/bin/env python
"""Step 2: for every shortlisted frame (16_wire/<sys>_cand.nc, all waters) find the real H-bond water wire catalytic-water -> (<=1 intermediate water) -> Asp180 carboxylate,
and the waters H-bonded to the side chain of residue 405.  H-bond = D-A <= 3.5 A and D-H...A >= 135 deg (as in the main analysis).
Writes 16_wire/<sys>_cand_eval.json (all frames) and prints the ranking."""
import os, sys, json, re, itertools, numpy as np, netCDF4
sys.path.insert(0, '/path/to/water_C1_occupancy/00_scripts')
os.environ['WC1_DIR'] = '/path/to/water_C1_occupancy'
import analyse_water_C1 as A
W = A.W; ang = A.angle
OFF = A.OFFSET; LIG = A.LIG
HB_D, HB_A = 3.5, 135.0

def hb_wp(Ow, H1, H2, P, dH):   # water (O,H1,H2 arrays n x 3) to polar atoms P (m x 3), dH = H atoms of the solute
    d = np.linalg.norm(Ow[:, None] - P[None], axis=-1); ok = d <= HB_D; r = np.zeros_like(ok)
    if not ok.any(): return r
    r |= ok & ((ang(Ow[:, None], H1[:, None], P[None]) >= HB_A) | (ang(Ow[:, None], H2[:, None], P[None]) >= HB_A))
    return r

def analyse(s):
    names, rl, ptr, resid, resname = A.parm_sections(f'{W}/16_wire/{s}_cand.step3_input.parm7')
    site = [int(x) for x in open(f'{W}/01_cpptraj/{s}_site_residues.txt').read().strip().split(',')]
    orig = site + [LIG[s]]; nprot = len(orig); pdbn = {i: orig[i] + OFF[s] for i in range(nprot - 1)}
    lig_i = nprot - 1; at = lambda ri, nm: [a for a in range(ptr[ri], ptr[ri + 1]) if names[a] == nm][0]
    C30, O19 = at(lig_i, 'C30'), at(lig_i, 'O19')
    ri = {p: i for i, p in pdbn.items()}
    acc180 = [a for a in range(ptr[ri[180]], ptr[ri[180] + 1]) if names[a] in ('OD1', 'OD2')]
    sc405 = [a for a in range(ptr[ri[405]], ptr[ri[405] + 1]) if names[a] in ('OD1', 'OD2', 'ND2')]
    bb405 = [at(ri[405], 'O')]
    wat = [i for i, r in enumerate(rl) if r == 'WAT']; wO = np.array([ptr[i] for i in wat])
    hyd = np.array([a for a in range(ptr[0], ptr[nprot]) if names[a][0] == 'H'])
    X = netCDF4.Dataset(f'{W}/16_wire/{s}_cand.nc').variables['coordinates']; fr = np.load(f'{W}/16_wire/{s}_cand_frames.npy')
    D = np.load(f'{W}/04_results/per_frame_{s}_merged.npz'); out = []
    for k in range(X.shape[0]):
        x = np.asarray(X[k], float); c = x[C30]; O = x[wO]; H1 = x[wO + 1]; H2 = x[wO + 2]
        near = np.where(np.linalg.norm(O - c, axis=1) < 14)[0]         # waters near the site (wire and 405 region are within ~12 A of C1)
        On, H1n, H2n = O[near], H1[near], H2[near]
        d = np.linalg.norm(On - c, axis=1); a_ = ang(On, c[None], x[O19][None])
        cand = np.where((d <= 4.0) & (a_ >= 135))[0]
        rec = dict(k=k + 1, merged_frame=int(fr[k]), strict=False)
        if len(cand):
            j = cand[np.argmin(d[cand])]; rec.update(cat=int(near[j]), cat_d=float(d[j]), cat_ang=float(a_[j]))
            # water-water H-bond adjacency among the near waters
            dd = np.linalg.norm(On[:, None] - On[None], axis=-1); cl = (dd <= HB_D) & (dd > 0.1); adj = np.zeros_like(cl)
            for Hs in (H1n, H2n): adj |= cl & (ang(On[:, None], Hs[:, None], On[None]) >= HB_A)
            adj = adj | adj.T
            h180 = hb_wp(On, H1n, H2n, x[acc180], None).any(1)           # waters donating to the Asp180 carboxylate
            path = None
            if h180[j]: path = [j]
            else:
                for m in np.where(adj[j] & h180)[0]: path = [j, int(m)]; break
            if path is not None:
                rec['strict'] = True; rec['wire'] = [int(near[p]) for p in path]
                rec['wire_d_C1'] = [float(d[p]) for p in path]
                rec['acc180'] = [names[a] for a in acc180 if hb_wp(On[[path[-1]]], H1n[[path[-1]]], H2n[[path[-1]]], x[[a]], None).any()]
            # waters at residue 405: side chain H-bond (water donates, or side-chain H donates) and backbone O
            Hs_all = x[hyd]
            def solute_hb(sel):
                res = np.zeros(len(On), bool)
                P = x[sel]; res |= hb_wp(On, H1n, H2n, P, None).any(1)
                for pj in range(len(sel)):
                    dh = np.linalg.norm(Hs_all - P[pj], axis=1); hs = Hs_all[dh < 1.2]
                    for h in hs: res |= (np.linalg.norm(On - P[pj], axis=1) <= HB_D) & (ang(P[pj][None], h[None], On) >= HB_A)
                return res
            w405 = np.where(solute_hb(sc405))[0]; w405b = np.where(solute_hb(bb405))[0]
            rec['w405_side'] = [int(near[p]) for p in w405]; rec['w405_side_dC1'] = [float(d[p]) for p in w405]
            rec['w405_bb'] = [int(near[p]) for p in w405b]
        out.append(rec)
    return out, names

for s in LIG:
    recs, _ = analyse(s); json.dump(recs, open(f'{W}/16_wire/{s}_cand_eval.json', 'w'))
    good = [r for r in recs if r['strict']]; both = [r for r in good if r['w405_side']]
    two = [r for r in good if len(r['wire']) == 2]
    print(f"{s}: {len(recs)} shortlisted; catalytic-water wire to Asp180 carboxylate confirmed in {len(good)} ({len(two)} with 2 waters); of these also a water H-bonded to the 405 side chain: {len(both)}")
    for r in sorted(both, key=lambda r: (-len(r['wire']), abs(r['cat_d'] - 3.8)))[:5]:
        print(f"    cand {r['k']:2d} merged frame {r['merged_frame']:5d}: cat d={r['cat_d']:.2f} A ang={r['cat_ang']:.0f}; wire {len(r['wire'])} water(s) -> Asp180 {r['acc180']}; wire d(C1)={[round(v,1) for v in r['wire_d_C1']]}; 405 side-chain waters {len(r['w405_side'])} (d(C1)={[round(v,1) for v in r['w405_side_dC1']]}); 405 backbone waters {len(r['w405_bb'])}")
