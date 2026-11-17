#!/usr/bin/env python
"""Water occupancy / catalytic-water analysis around the anomeric carbon (C30, subsite -1).

Input  : 02_reduced/<sys>_<rep>_site.nc (+ <sys>_site.step3_input.parm7, <sys>_<rep>_closest.dat.gz)
         every 10th frame of the 2-us trajectory (200 ps), 80 waters closest to C30, sorted by distance.
Output : 04_results/per_frame_<sys>_<rep>.npz, 04_results/summary_per_replica.csv,
         04_results/summary_per_system.csv, 04_results/hbond_partners.csv, 04_results/residence.csv

Definitions (ligand atom names of the cellooctaose residue in the topology):
  C30 = anomeric carbon of the -1 glucose, O19 = scissile glycosidic O (-1/+1), O25 = ring O of -1.
  n(r)            : number of water O within r = 3.0, 3.5, 4.0 A of C30.
  Catalytic water : water O within 4.0 A of C30 on the face opposite the leaving group,
                    angle Ow...C30-O19 >= 135 deg (closest such water per frame).
  Closest water   : d and angle of the water O closest to C30 (no angle filter) -> unbiased geometry.
  Attack-ready    : catalytic water with d(Ow-C30) <= 3.5 A and angle >= 150 deg.
  H-bond          : donor-acceptor <= 3.5 A and D-H...A angle >= 135 deg (water H's explicit;
                    for protonated Asp226 its HD2 is used as donor).
  Bridge X-W-C1   : a water within 4.0 A of C30 that is H-bonded to residue X.
  Residence       : continuous presence of the same water molecule within 3.5 A of C30
                    (gaps of one frame = 200 ps tolerated).
"""
import gzip, re, sys, os
import numpy as np, pandas as pd, netCDF4

W = os.environ.get('WC1_DIR', '/path/to/water_C1_occupancy')
SYSTEMS = ['HiCel6A', 'HiCel6AD405N', 'HiCel6AD416A']
REPS = ['rep1', 'rep2', 'rep3']
DT_PS = 200.0
SEG = int(os.environ.get('WC1_SEG', 10**9))   # frames per replica in a merged trajectory: residence episodes never span a junction
LIG = {'HiCel6A': 362, 'HiCel6AD405N': 365, 'HiCel6AD416A': 361}
OFFSET = {'HiCel6A': 89, 'HiCel6AD405N': 86, 'HiCel6AD416A': 90}   # PDB number = topology + offset (checked for D180/D226/D405/D416)
KEY = {180: 'Asp180', 226: 'Asp226', 405: 'Asp405', 416: 'Asp416'}  # PDB numbering


def parm_sections(path):
    txt = open(path).read()
    def sec(name):
        m = re.search(r'%FLAG ' + name + r'\s*\n%FORMAT\((\d+)([aIE])(\d+)[^)]*\)\s*\n(.*?)(?=%FLAG|\Z)', txt, re.S)
        w = int(m.group(3)); s = m.group(4).replace('\n', '')
        return [s[i:i + w].strip() for i in range(0, len(s), w)]
    names = sec('ATOM_NAME'); rl = sec('RESIDUE_LABEL')
    ptr = [int(p) - 1 for p in sec('RESIDUE_POINTER')] + [len(names)]
    resid = np.zeros(len(names), int); resname = [''] * len(names)
    for i in range(len(rl)):
        resid[ptr[i]:ptr[i + 1]] = i;
        for a in range(ptr[i], ptr[i + 1]): resname[a] = rl[i]
    return names, rl, ptr, resid, resname


def angle(a, b, c):
    """angle a-b-c in degrees; arrays (...,3)."""
    u = a - b; v = c - b
    cs = (u * v).sum(-1) / (np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1))
    return np.degrees(np.arccos(np.clip(cs, -1, 1)))


def setup(sys_):
    names, rl, ptr, resid, resname = parm_sections(f'{W}/02_reduced/{sys_}_site.step3_input.parm7')
    # original residue numbers are lost after strip -> recover from the stripped residue order
    site = [int(x) for x in open(f'{W}/01_cpptraj/{sys_}_site_residues.txt').read().strip().split(',')]
    orig = site + [LIG[sys_]]                       # strip keeps order: site residues, ligand, waters
    nprot = len(orig)
    pdbnum = {i: orig[i] + OFFSET[sys_] for i in range(nprot - 1)}
    lig_i = nprot - 1
    assert rl[lig_i] == 'BGLC', rl[lig_i]
    def atom(ri, nm):
        for a in range(ptr[ri], ptr[ri + 1]):
            if names[a] == nm: return a
        raise KeyError((ri, nm))
    L = {nm: atom(lig_i, nm) for nm in ('C30', 'O19', 'O25')}
    wat = [i for i, r in enumerate(rl) if r == 'WAT']
    wO = np.array([ptr[i] for i in wat]); wH1 = wO + 1; wH2 = wO + 2
    assert names[wO[0]] == 'O' and names[wH1[0]] == 'H1'
    # polar protein/ligand atoms for H-bond partner census
    polar = [a for a in range(ptr[0], ptr[nprot]) if names[a][0] in 'ON']
    # hydrogens bonded to polar atoms: assume H atoms immediately following within residue named H*
    # (used only for protein donors -> take all H within 1.2 A at run time)
    hyd = [a for a in range(ptr[0], ptr[nprot]) if names[a][0] == 'H']
    lab = []
    for a in polar:
        ri = resid[a]
        lab.append(f'{rl[ri]}{pdbnum[ri]}:{names[a]}' if ri != lig_i else f'Lig:{names[a]}')
    key = {}
    for ri, pn in pdbnum.items():
        if pn in KEY:
            at = [a for a in range(ptr[ri], ptr[ri + 1]) if names[a] in ('OD1', 'OD2', 'ND2')]
            if at: key[KEY[pn]] = at          # Ala416 (D416A) has no polar side chain -> skipped
    return dict(names=names, L=L, wO=wO, wH1=wH1, wH2=wH2, polar=np.array(polar), hyd=np.array(hyd),
                lab=np.array(lab), key=key, pdbnum=pdbnum, resid=resid, lig_i=lig_i)


def hbond_matrix(X, info, wsel):
    """H-bonds between selected waters (indices into water list) and polar solute atoms in one frame.
    returns bool (len(wsel), npolar)."""
    O = X[info['wO'][wsel]]; H1 = X[info['wH1'][wsel]]; H2 = X[info['wH2'][wsel]]
    P = X[info['polar']]
    d = np.linalg.norm(O[:, None] - P[None], axis=-1)
    hb = np.zeros_like(d, bool)
    close = d <= 3.5
    if not close.any(): return hb
    # water as donor
    a1 = angle(O[:, None], H1[:, None], P[None]); a2 = angle(O[:, None], H2[:, None], P[None])
    hb |= close & ((a1 >= 135) | (a2 >= 135))
    # solute as donor: any H within 1.2 A of polar atom
    Hs = X[info['hyd']]
    for j in np.where(close.any(0))[0]:
        dh = np.linalg.norm(Hs - P[j], axis=1); hs = Hs[dh < 1.2]
        for h in hs:
            ang = angle(P[j][None], h[None], O)
            hb[:, j] |= close[:, j] & (ang >= 135)
    return hb


def water_water_hb(X, info, wsel):
    """H-bonds between selected waters and all waters (O-O <= 3.5 A, O-H...O >= 135 deg, either direction)."""
    O = X[info['wO']]; H1 = X[info['wH1']]; H2 = X[info['wH2']]
    d = np.linalg.norm(O[wsel][:, None] - O[None], axis=-1)
    close = (d <= 3.5) & (d > 0.1)
    hb = np.zeros_like(close)
    if not close.any(): return hb
    a = O[wsel][:, None]; b = O[None]
    for Hs, don in ((H1[wsel][:, None], a), (H2[wsel][:, None], a)):
        hb |= close & (angle(don, Hs, b) >= 135)
    for Hs in (H1[None], H2[None]):
        hb |= close & (angle(b, Hs, a) >= 135)
    hb[np.arange(len(wsel)), wsel] = False
    return hb


def residence(ids_in, max_gap=1):
    """ids_in: list of sets (water ids within cutoff per frame). returns list of episode lengths (frames)."""
    open_ = {}; last = {}; ep = []
    for f, ids in enumerate(ids_in):
        for w in ids:
            if w in last and f - last[w] <= max_gap + 1:
                last[w] = f
            else:
                if w in last: ep.append(last[w] - open_[w] + 1)
                open_[w] = f; last[w] = f
        # close episodes that exceeded the gap
        for w in [w for w in last if f - last[w] > max_gap]:
            ep.append(last[w] - open_[w] + 1); del last[w]; del open_[w]
    ep += [last[w] - open_[w] + 1 for w in last]
    return np.array(ep)


def analyse(sys_, rep):
    info = setup(sys_)
    nc = netCDF4.Dataset(f'{W}/02_reduced/{sys_}_{rep}_site.nc'); X = nc.variables['coordinates']
    nf = X.shape[0]
    cl = pd.read_csv(f'{W}/02_reduced/{sys_}_{rep}_closest.dat.gz', sep=r'\s+', comment='#', header=None,
                     usecols=[0, 1], names=['f', 'mol'])
    mol = cl['mol'].values.reshape(nf, -1)
    assert mol.shape[1] == len(info['wO'])
    L = info['L']; keys = [k for k in KEY.values() if k in info['key']]
    out = {k: np.full(nf, np.nan) for k in ['n30', 'n35', 'n40', 'cat_d', 'cat_ang', 'cat_mol', 'near_d', 'near_ang']}
    for k in keys:
        out[f'dmin_{k}'] = np.full(nf, np.nan)           # direct min distance residue polar O/N - C30
        out[f'bridge_{k}'] = np.zeros(nf, bool)           # any water within 4 A of C30 H-bonded to residue
        out[f'cathb_{k}'] = np.zeros(nf, bool)            # catalytic water H-bonded to residue
        out[f'wire2_{k}'] = np.zeros(nf, bool)            # C1-W(<=4 A)-W-residue or C1-W-residue (<=2 waters)
    out['cat_nhb'] = np.zeros(nf)
    partner = np.zeros(len(info['polar']))
    ncat = 0
    ids35 = []
    CH = 2000
    for f0 in range(0, nf, CH):
        blk = np.asarray(X[f0:f0 + CH], dtype=np.float64)
        for k_ in range(blk.shape[0]):
            f = f0 + k_; x = blk[k_]
            c = x[L['C30']]; O = x[info['wO']]
            d = np.linalg.norm(O - c, axis=1)
            out['n30'][f] = (d <= 3.0).sum(); out['n35'][f] = (d <= 3.5).sum(); out['n40'][f] = (d <= 4.0).sum()
            ids35.append(set(mol[f][d <= 3.5]))
            ang = angle(O, c[None], x[L['O19']][None])
            j0 = int(np.argmin(d)); out['near_d'][f] = d[j0]; out['near_ang'][f] = ang[j0]   # closest water, any angle
            for k in keys:
                out[f'dmin_{k}'][f] = np.linalg.norm(x[info['key'][k]] - c, axis=1).min()
            near = np.where(d <= 4.0)[0]
            if len(near):
                hb = hbond_matrix(x, info, near)
                ww = water_water_hb(x, info, near)                 # (len(near), nwat)
                allw = np.where(ww.any(0))[0]
                hb2 = hbond_matrix(x, info, allw) if len(allw) else None
                for k in keys:
                    cols = np.isin(info['polar'], info['key'][k])
                    out[f'bridge_{k}'][f] = hb[:, cols].any()
                    w2 = hb2[:, cols].any(1) if hb2 is not None else np.zeros(0, bool)
                    out[f'wire2_{k}'][f] = out[f'bridge_{k}'][f] or (len(allw) > 0 and w2.any())
                cand = near[ang[near] >= 135]
                if len(cand):
                    j = cand[np.argmin(d[cand])]; jj = int(np.where(near == j)[0][0])
                    out['cat_d'][f] = d[j]; out['cat_ang'][f] = ang[j]; out['cat_mol'][f] = mol[f][j]
                    out['cat_nhb'][f] = hb[jj].sum(); partner += hb[jj]; ncat += 1
                    for k in keys:
                        cols = np.isin(info['polar'], info['key'][k])
                        out[f'cathb_{k}'][f] = hb[jj, cols].any()
        print(sys_, rep, f0 + blk.shape[0], '/', nf, flush=True)
    segs = list(range(0, nf, SEG))
    eps = [residence(ids35[a:a + SEG]) for a in segs]
    ep = np.concatenate(eps) * DT_PS; res_rep = np.concatenate([np.full(len(e), i) for i, e in enumerate(eps)]).astype(int)
    # residence of catalytic-water identities (same molecule staying catalytic), broken at every replica junction
    cm = out['cat_mol']; cat_ep = []; cat_rep = []
    run = 0; prev = None
    for f, v in enumerate(cm):
        if f % SEG == 0:
            if run: cat_ep.append(run); cat_rep.append((f - 1) // SEG)
            run = 0; prev = None
        if not np.isnan(v) and v == prev: run += 1
        else:
            if run: cat_ep.append(run); cat_rep.append((f - 1) // SEG)
            run = 1 if not np.isnan(v) else 0
        prev = v
    if run: cat_ep.append(run); cat_rep.append((nf - 1) // SEG)
    cat_ep = np.array(cat_ep) * DT_PS; cat_rep = np.array(cat_rep, int)
    out['rep_id'] = (np.arange(nf) // SEG).astype(np.int8)
    np.savez_compressed(f'{W}/04_results/per_frame_{sys_}_{rep}.npz', **out, res_ep_ps=ep, cat_ep_ps=cat_ep, res_ep_rep=res_rep, cat_ep_rep=cat_rep)
    pf = pd.DataFrame({'label': info['lab'], 'freq_given_catwater': partner / max(ncat, 1)})
    pf['system'] = sys_; pf['rep'] = rep
    pf.to_csv(f'{W}/04_results/hbond_partners_{sys_}_{rep}.csv', index=False)


if __name__ == '__main__':
    args = sys.argv[1:]
    todo = [(s, r) for s in SYSTEMS for r in REPS] if not args else [tuple(args)]
    for s, r in todo:
        analyse(s, r)
