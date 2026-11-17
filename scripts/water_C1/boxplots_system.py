#!/usr/bin/env python3
"""Figure S12: box plots of the water near the anomeric carbon C1, ONE BOX PER SYSTEM (HiCel6A, D416A, D405N), all 30 000 merged frames (200 ps per frame).
  panels: water within 4 A of C1; distance and angle of the catalytic-position water; water wire to Asp180.
Usage: set WC1_DIR to the analysis folder (contains 04_results/) and optionally WC1_OUT for the output folder.
  continuous quantities : every frame
  yes/no quantities    : percentage of frames in consecutive blocks of 50 frames (10 ns) -> box = distribution of the block percentages
Whiskers at the 5th/95th percentile, no outliers."""
import os, re, csv
import pandas as pd
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

plt.rcParams.update({'font.family': 'sans-serif', 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.bfit': 'Arial:italic:bold', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'pdf.fonttype': 42, 'ps.fonttype': 42})
W = os.environ.get('WC1_DIR', '/path/to/water_C1_occupancy'); O = os.environ.get('WC1_OUT', W)
OUT_SYS = f'{O}/10_boxplots_system'; OUT_DATA = f'{O}/10_boxplot_data'; OUT = OUT_SYS
for k in ('merged',):
    os.makedirs(f'{OUT_SYS}/{k}', exist_ok=True); os.makedirs(f'{OUT_DATA}/{k}', exist_ok=True)
SYSCOL = ['#4E79A7', '#59A14F', '#F28E5B']     # HiCel6A, D416A, D405N
SYSTEMS = [('HiCel6A', r'$\mathbfit{Hi}\mathbf{Cel6A}$'), ('HiCel6AD416A', r'$\mathbfit{Hi}\mathbf{Cel6A}^{\mathbf{D416A}}$'),
           ('HiCel6AD405N', r'$\mathbfit{Hi}\mathbf{Cel6A}^{\mathbf{D405N}}$')]
REPS = ['merged']; COLORS = ['#4E79A7', '#F28E5B', '#59A14F']; WIN = 500

N1 = lambda s, r: np.load(f'{W}/04_results/per_frame_{s}_merged.npz')
cache = {}


def arr(kind, s, r, key):
    k = (kind, s, r)
    if k not in cache: cache[k] = N1(s, r)
    return np.asarray(cache[k][key], float)


def abbrev_label(s):
    """Use consistent manuscript atom notation in all plotted labels."""
    for old, new in (
        ("Angle Ow···C1–O19", "Angle O$_w$···C$_{anom}$–O$_{glyc}$"),
        ("Ow···C1–O19", "O$_w$···C$_{anom}$–O$_{glyc}$"),
        ("d(Ow···C1)", "d(O$_w$···C$_{anom}$)"),
        ("Ow···C1", "O$_w$···C$_{anom}$"),
        ("O19", "O$_{glyc}$"), ("C1", "C$_{anom}$"),
        ("Catalytic-position water", "Cat. W"), ("catalytic-position water", "cat. W"),
        ("catalytic water", "cat. W"), ("Waters", "W"), ("waters", "W"),
        ("Water", "W"), ("Residence time of water", "Residence, W"),
        ("Nearest water", "Nearest W"), ("H-bonded water identities", "H-bonded W identities"),
        ("(% of frames)", "(% frames)"),
        (", 6 Å shell", ", 6 Å"), ("Waiting", ""),
    ):
        s = s.replace(old, new)
    return s


def auto_ylim(dicts):
    """y range from the whiskers: lowest 5th and highest 95th percentile of the systems, plus 8 % padding"""
    vs = [v[~np.isnan(v)] for d in dicts for v in d.values()]; vs = [v for v in vs if len(v)]
    lo = min(np.percentile(v, 5) for v in vs); hi = max(np.percentile(v, 95) for v in vs); pad = 0.08 * (hi - lo if hi > lo else 1.0)
    return (max(lo - pad, 0) if lo >= 0 else lo - pad, hi + pad)


CHECKS = []


def _check(fig, name, series, ylim, log):
    """automatic test of one figure: y label and legend inside the canvas, 5th-95th percentile whiskers inside the y range, nothing empty"""
    fig.canvas.draw(); rd = fig.canvas.get_renderer(); fb = fig.bbox
    inside = lambda bb: bb.x0 >= fb.x0 - 1 and bb.y0 >= fb.y0 - 1 and bb.x1 <= fb.x1 + 1 and bb.y1 <= fb.y1 + 1
    lab_ok = inside(fig._supylabel.get_window_extent(rd)) if getattr(fig, '_supylabel', None) is not None else False
    leg_ok = all(inside(l.get_window_extent(rd)) for l in fig.legends)
    if getattr(fig, '_supylabel', None) is not None:
        sb = fig._supylabel.get_window_extent(rd)
        lab_ok = lab_ok and not any(sb.overlaps(tk.get_window_extent(rd)) for tk in fig.axes[0].get_yticklabels() if tk.get_text())    # label must not touch the tick numbers
    lo = min(np.percentile(v, 5) for v in series if len(v)); hi = max(np.percentile(v, 95) for v in series if len(v))
    ylo = ylim[0]; yhi = ylim[1]
    whisk_ok = lo >= ylo - 1e-9 and hi <= yhi + 1e-9
    CHECKS.append(dict(figure=name, ylabel_inside=lab_ok, legend_inside=leg_ok, whiskers_inside_ylim=whisk_ok, lowest_whisker=f'{lo:.5g}', highest_whisker=f'{hi:.5g}', ylim_low=f'{ylo:.5g}', ylim_high=f'{yhi:.5g}',
                       n_series=len(series), empty_series=sum(1 for v in series if len(v) == 0)))


def write_checks(path):
    with open(path, 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.DictWriter(fh, fieldnames=list(CHECKS[0].keys())); w.writeheader(); w.writerows(CHECKS)
    bad = [c for c in CHECKS if not (c['ylabel_inside'] and c['legend_inside'] and c['whiskers_inside_ylim']) or c['empty_series']]
    print('checked', len(CHECKS), 'figures;', len(bad), 'with a problem')
    for c in bad: print('  PROBLEM', c['figure'], {k: c[k] for k in ('ylabel_inside', 'legend_inside', 'whiskers_inside_ylim', 'lowest_whisker', 'highest_whisker', 'ylim_low', 'ylim_high', 'empty_series')})


def one_system(data, stem, tag, ylabel, ylim, log, mean_marker):
    """one box per system all frames, blocks or episodes of the system."""
    fig, ax = plt.subplots(figsize=(7.6, 5.15))
    pooled = [np.concatenate([data[(s, r)][~np.isnan(data[(s, r)])] for r in REPS]) for s, _ in SYSTEMS]
    bp = ax.boxplot(pooled, positions=[1, 2, 3], widths=.34, patch_artist=True, showfliers=False, whis=(5, 95),
                    medianprops={'color': 'black', 'linewidth': 2.0}, whiskerprops={'color': '#333333', 'linewidth': 1.5}, capprops={'color': '#333333', 'linewidth': 1.5})
    for patch, c in zip(bp['boxes'], SYSCOL): patch.set_facecolor(c); patch.set_alpha(.78); patch.set_edgecolor('#333333')
    if mean_marker:
        ax.scatter([1, 2, 3], [v.max() for v in pooled], marker='D', color='#C00000', s=45, zorder=7)
        ax.scatter([0.035], [0.965], transform=ax.transAxes, marker='D', color='#C00000', s=38, zorder=9, clip_on=True)
        ax.text(0.085, 0.965, 'longest stay', transform=ax.transAxes, ha='left', va='center',
                fontsize=12, color='#222222', zorder=9,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.82, pad=1.5))
    ax.set_xticks([1, 2, 3], ['', '', '']); ax.set_xlim(.4, 3.6)
    for i, (_, tt) in enumerate(SYSTEMS):      # system names on one common baseline (labels with a superscript would otherwise sit lower)
        ax.text(i + 1, -0.085, tt, transform=ax.get_xaxis_transform(), ha='center', va='baseline', fontsize=17)
    if log:
        ax.set_yscale('log'); ax.set_ylim(*ylim)
        cand = [0.2, 1, 10, 100, 2000]
        ax.yaxis.set_major_locator(FixedLocator([t for t in cand if ylim[0] <= t <= ylim[1]])); ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    else:
        ax.set_ylim(*ylim)
    ax.grid(axis='y', color='#d9d9d9', linewidth=.8); ax.set_axisbelow(True)
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_linewidth(1.4); sp.set_color('#333333')
    ax.tick_params(axis='both', which='major', width=1.1, length=4, color='#333333', labelsize=15, pad=4)
    for tk in ax.get_yticklabels(): tk.set_fontsize(17); tk.set_fontweight('bold')
    lab = abbrev_label(ylabel); fs = min(17, 345 / (0.60 * len(lab)))
    fig.supylabel(lab, fontsize=fs, fontweight='bold', x=.02, y=.60, va='center')
    fig.subplots_adjust(left=.17, right=.99, top=.96, bottom=.29)
    _check(fig, f'Figure_{stem}', pooled, ylim, log)
    for ext, kw in (('png', {'dpi': 300}), ('pdf', {})): fig.savefig(f'{OUT_SYS}/{tag}/Figure_{stem}.{ext}', facecolor='white', **kw)
    plt.close(fig)


def _kind(data, log):
    n = max(len(v) for v in data.values())
    return 'episode' if log else ('block' if n <= 600 else 'frame')


def data_csv(data, stem, tag, ylabel, log, path):
    """ONE table per figure: the data of every system (frame by frame, block by block, or one row per residence episode), then an empty column and the summary
    (n, mean, SD, percentiles) per system."""
    kind = _kind(data, log); cols = {}; raw = {}
    for s, _ in SYSTEMS:
        v = pd.Series(np.asarray(data[(s, REPS[0])], float)); cols[s] = v; raw[s] = v.values
    df = pd.DataFrame(cols); df.insert(0, f'index_{kind}', np.arange(1, len(df) + 1))
    def st(v):
        v = v[~np.isnan(v)]
        if not len(v): return [0] + [np.nan] * 9
        q = np.percentile(v, [5, 25, 50, 75, 95]); return [len(v), v.mean(), v.std(ddof=1) if len(v) > 1 else np.nan, *q, v.min(), v.max()]
    sd = {'statistic': ['n_values', 'mean', 'sd', 'p5', 'p25', 'median', 'p75', 'p95', 'min', 'max']}
    for s, _ in SYSTEMS: sd[f'summary_{s}'] = st(raw[s])
    sdf = pd.DataFrame(sd); sdf['y_axis_label'] = [ylabel] + [''] * 9; sdf.insert(0, '', '')
    pd.concat([df, sdf], axis=1).to_csv(path, index=False, float_format='%.6g', encoding='utf-8-sig')


def draw(data_full, data_win, stem, ylabel, log=False, ylim=None):
    """ylim: fixed range for the type of analysis (percent: 0-100; residence: 0.015 ns to the full 2 us; pairs of the same quantity share one range)"""
    if ylim is None: ylim = (0.1, 2.5 * max(v.max() for v in data_full.values() if len(v))) if log else auto_ylim([data_full])   # residence: 0.2 ns (one frame) up to 2.5 x the longest stay
    for tag, data in (('merged', data_full),):
        one_system(data, stem, tag, ylabel, ylim, log, mean_marker=log)
        data_csv(data, stem, tag, ylabel, log, f'{OUT_DATA}/{tag}/{stem}.csv')
    print('saved', stem)


def per_frame(kind, key, fn=None):
    full, win = {}, {}
    for s, _ in SYSTEMS:
        for r in REPS:
            v = fn(lambda k: arr(kind, s, r, k)) if fn else arr(kind, s, r, key); full[(s, r)] = v; win[(s, r)] = v[:WIN]
    return full, win


def fractions(kind, fn):
    full, win = {}, {}
    for s, _ in SYSTEMS:
        for r in REPS:
            b = fn(lambda k: arr(kind, s, r, k)).astype(float)
            full[(s, r)] = 100 * b.reshape(-1, 50).mean(1); win[(s, r)] = 100 * b[:WIN].reshape(-1, 10).mean(1)
    return full, win


def main():
    for stem, key, ylabel in (('C1_catalytic_water_distance', 'cat_d', 'Catalytic-position water d(Ow···C1) (Å)'),
                              ('C1_catalytic_water_angle', 'cat_ang', 'Angle Ow···C1–O19 of catalytic water (°)')):
        f, w = per_frame('c1', key); draw(f, w, stem, ylabel)
    for stem, ylabel, fn in (('C1_occupancy_within_4A', 'Water within 4 Å of C1 (% of frames)', lambda g: g('n40') >= 1),
                             ('C1_water_wire_to_Asp180', 'Water wire (≤ 2 waters) C1–Asp180 (% of frames)', lambda g: g('wire2_Asp180') > 0)):
        f, w = fractions('c1', fn); draw(f, w, stem, ylabel, ylim=(-2, 102))     # percent-of-frames figures: 0-100


if __name__ == '__main__':
    main()
    write_checks(f'{OUT}/figure_checks.csv')
