#!/usr/bin/env python3
"""
Generate v2 filled Figure_2 SVGs with measurement error bars (n=10) AND
95% CI band of the linear regression fit on top of the existing data points.

Inputs:
  - Single-measurement totals: Figure1/python/Downsampling_Results/<ds>/<ds>_metrics_v8_{10,11}.json
  - n=10 (or n=97) repetition stats: Figure3/python/Outputs/size_sweep/n10_<ds>/n100_<ds>.json
    plus Figure3/python/Outputs/n100_downsampling/n100_Jikji6.json
  Both repetition sources are at ratio = 0.188 (downsampling_factor 5.33).

Method:
  - Per-point time/position taken from v8_10/11 (full FASTQ).
  - Measurement std at full FASTQ is estimated by CV scaling:
        CV = std_down / mean_down
        std_full ≈ CV × time_full
    where time_full comes from v8_10/11. The n=10 reps measure the *relative*
    noise of the pipeline, which is approximately scale-invariant.
  - 95% CI band of the regression line uses the standard OLS formula:
        SE(y_hat | x) = s_e * sqrt(1/n + (x - x_bar)^2 / SS_xx)
        CI = y_hat ± t(0.025, n-2) * SE(y_hat | x)
    over n=13 data points.

Outputs:
  Manuscript_Drafts_v2/Figure_2/Figure_2a_Time_Regression_filled_v2.svg
"""
import json
import math
import pathlib
from typing import List, Dict, Optional, Tuple

ROOT = pathlib.Path('/Users/young/Desktop/DNA ink3')
RES_DIR = ROOT / 'Figure1' / 'python' / 'Downsampling_Results'
SWEEP_DIR = ROOT / 'Figure3' / 'python' / 'Outputs' / 'size_sweep'
N100_DOWNSAMPLING = ROOT / 'Figure3' / 'python' / 'Outputs' / 'n100_downsampling' / 'n100_Jikji6.json'
OUT_DIR = ROOT / 'Manuscript_Drafts_v2' / 'Figure_2'

DATASET_ORDER = [
    'Jikji5', 'Jikji1', 'Jikji2', 'Jikji6', 'Jikji3', 'Jikji4', 'Jikji7',
    'seokga', 'goldengate', 'pisa', 'liberty', 'venus', 'tripitaka',
]

# Datasets with repetition data + source pointer
REP_SOURCES = {
    'Jikji6':    N100_DOWNSAMPLING,
    'seokga':    SWEEP_DIR / 'n10_seokga'    / 'n100_seokga.json',
    'pisa':      SWEEP_DIR / 'n10_pisa'      / 'n100_pisa.json',
    'liberty':   SWEEP_DIR / 'n10_liberty'   / 'n100_liberty.json',
    'venus':     SWEEP_DIR / 'n10_venus'     / 'n100_venus.json',
    'tripitaka': SWEEP_DIR / 'n10_tripitaka' / 'n100_tripitaka.json',
}

# Student-t critical value t_{0.025, df} (two-sided 95% CI)
T_CRIT = {
    1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
    6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228,
    11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131,
    20: 2.086, 30: 2.042, 60: 2.000, 120: 1.980,
}


def t_critical(df: int) -> float:
    if df in T_CRIT: return T_CRIT[df]
    keys = sorted(T_CRIT.keys())
    # linear interp
    for k1, k2 in zip(keys, keys[1:]):
        if k1 <= df <= k2:
            f = (df - k1) / (k2 - k1)
            return T_CRIT[k1] + f * (T_CRIT[k2] - T_CRIT[k1])
    return 1.96  # fall back to z


def load_full_measurements() -> List[Dict]:
    """Load full-FASTQ single-measurement totals (n=1 each)."""
    rows = []
    for ds in DATASET_ORDER:
        v11 = RES_DIR / ds / f'{ds}_metrics_v8_11.json'
        v10 = RES_DIR / ds / f'{ds}_metrics_v8_10.json'
        path = v11 if v11.exists() else (v10 if v10.exists() else None)
        if path is None: continue
        with open(path) as f: m = json.load(f)
        tc = m.get('timecell', {}) or {}
        total_s = sum(v for v in tc.values() if isinstance(v, (int, float)))
        rows.append({
            'name': ds,
            'fastq_GB': m.get('fastq_size_bytes', 0) / 1e9,
            'total_s': total_s,
            'total_min': total_s / 60.0,
        })
    return rows


def load_repetition_stats() -> Dict[str, Dict]:
    """Load per-dataset n=10/97 mean and std (in *downsampled* seconds)."""
    out = {}
    for ds, path in REP_SOURCES.items():
        if not path.exists(): continue
        with open(path) as f: m = json.load(f)
        s = m['summary']['decode_time_total']
        meta = m['meta']
        out[ds] = {
            'n_trials': s['n'],
            'down_mean_s': s['mean'],
            'down_std_s': s['std'],
            'down_p5_s': s['p5'],
            'down_p95_s': s['p95'],
            'cv': s['std'] / s['mean'] if s['mean'] else 0.0,
            'downsampling_factor': meta['downsampling_factor'],
        }
    return out


def linreg(xs: List[float], ys: List[float]) -> Tuple[float, float, float, float, float]:
    """Return (slope, intercept, R^2, s_e, SS_xx).
    s_e is the residual standard error (sqrt of MSE)."""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    ss_xy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    ss_xx = sum((x - mx) ** 2 for x in xs)
    slope = ss_xy / ss_xx
    intercept = my - slope * mx
    yhat = [slope * x + intercept for x in xs]
    ss_res = sum((y - yh) ** 2 for y, yh in zip(ys, yhat))
    ss_tot = sum((y - my) ** 2 for y in ys)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    df_res = max(n - 2, 1)
    s_e = math.sqrt(ss_res / df_res)
    return slope, intercept, r2, s_e, ss_xx


def ci_halfwidth(x: float, x_bar: float, ss_xx: float, n: int, s_e: float) -> float:
    """Half-width of the 95% CI of the mean response at x."""
    se_yhat = s_e * math.sqrt(1.0 / n + (x - x_bar) ** 2 / ss_xx)
    return t_critical(n - 2) * se_yhat


SMALL_CLUSTER = {'Jikji1','Jikji2','Jikji3','Jikji4','Jikji5','Jikji6','Jikji7','seokga'}


def render_2a(rows: List[Dict], reps: Dict[str, Dict]) -> str:
    # Plot bounds
    PX, PR = 80, 600
    PY_TOP, PY_BOT = 45, 370
    GB_MAX, MIN_MAX = 24, 110

    def sx(gb): return PX + gb / GB_MAX * (PR - PX)
    def sy(m):  return PY_BOT - m / MIN_MAX * (PY_BOT - PY_TOP)

    xs_gb = [r['fastq_GB'] for r in rows]
    ys_min = [r['total_min'] for r in rows]
    slope, intercept, r2, s_e, ss_xx = linreg(xs_gb, ys_min)
    n = len(rows)
    x_bar = sum(xs_gb) / n

    # ----- Per-point full-FASTQ std (estimated via CV scaling) -----
    for r in rows:
        rep = reps.get(r['name'])
        if rep is None:
            r['std_min'] = None
            r['n_trials'] = None
        else:
            r['std_min'] = rep['cv'] * r['total_min']
            r['n_trials'] = rep['n_trials']

    # ----- CI band points along regression line (50 points, clipped) -----
    band_n = 50
    xs_band = [GB_MAX * i / (band_n - 1) for i in range(band_n)]
    band_lo, band_hi = [], []
    for x in xs_band:
        yhat = slope * x + intercept
        hw = ci_halfwidth(x, x_bar, ss_xx, n, s_e)
        band_lo.append(max(0.0, min(MIN_MAX, yhat - hw)))
        band_hi.append(max(0.0, min(MIN_MAX, yhat + hw)))

    poly_pts = []
    for x, y in zip(xs_band, band_hi):
        poly_pts.append(f'{sx(x):.1f},{sy(y):.1f}')
    for x, y in zip(reversed(xs_band), reversed(band_lo)):
        poly_pts.append(f'{sx(x):.1f},{sy(y):.1f}')
    band_polygon = ' '.join(poly_pts)

    y_ticks = [0, 20, 40, 60, 80, 100]
    x_ticks = [0, 5, 10, 15, 20]

    parts = []
    parts.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 660 430" '
                 'font-family="Arial, Helvetica, sans-serif">')
    parts.append('<defs><clipPath id="plotArea2a">'
                 f'<rect x="{PX}" y="{PY_TOP}" width="{PR-PX}" height="{PY_BOT-PY_TOP}"/>'
                 '</clipPath></defs>')
    parts.append('<rect width="660" height="430" fill="white"/>')
    parts.append('<text x="340" y="22" font-size="16.5" fill="#1A1A1A" '
                 'text-anchor="middle" font-weight="600">'
                 'Total Decoding Time vs. FASTQ File Size</text>')
    parts.append('<text x="340" y="38" font-size="11.5" fill="#4B5563" '
                 'text-anchor="middle">'
                 'Python serial baseline — n=13 datasets; error bars from n≥10 trials; shaded band = 95% CI of fit'
                 '</text>')

    # Axes
    parts.append(f'<line x1="{PX}" y1="{PY_TOP}" x2="{PX}" y2="{PY_BOT}" '
                 'stroke="#1A1A1A" stroke-width="2.03"/>')
    parts.append(f'<line x1="{PX}" y1="{PY_BOT}" x2="{PR}" y2="{PY_BOT}" '
                 'stroke="#1A1A1A" stroke-width="2.03"/>')
    parts.append('<text x="28" y="210" font-size="17.3" fill="#1A1A1A" '
                 'text-anchor="middle" transform="rotate(-90, 28, 210)">'
                 'Total Decoding Time (min)</text>')
    parts.append('<text x="340" y="405" font-size="17.3" fill="#1A1A1A" '
                 'text-anchor="middle">FASTQ File Size (GB)</text>')

    # Y ticks
    for v in y_ticks:
        y = sy(v)
        if v > 0:
            parts.append(f'<line x1="{PX}" y1="{y:.1f}" x2="{PR}" y2="{y:.1f}" '
                         'stroke="#EBEBEB" stroke-width="1.01"/>')
        parts.append(f'<line x1="{PX-5}" y1="{y:.1f}" x2="{PX}" y2="{y:.1f}" '
                     'stroke="#1A1A1A" stroke-width="1.69"/>')
        parts.append(f'<text x="{PX-8}" y="{y+4:.1f}" font-size="14.7" '
                     f'fill="#4B5563" text-anchor="end">{v}</text>')

    # X ticks
    for v in x_ticks:
        x = sx(v)
        parts.append(f'<line x1="{x:.1f}" y1="{PY_BOT}" x2="{x:.1f}" '
                     f'y2="{PY_BOT+5}" stroke="#1A1A1A" stroke-width="1.69"/>')
        parts.append(f'<text x="{x:.1f}" y="{PY_BOT+16}" font-size="14.7" '
                     f'fill="#4B5563" text-anchor="middle">{v}</text>')
        if v > 0:
            parts.append(f'<line x1="{x:.1f}" y1="{PY_TOP}" x2="{x:.1f}" '
                         f'y2="{PY_BOT}" stroke="#EBEBEB" stroke-width="1.01"/>')

    # --- 95% CI band BEHIND regression line (translucent pastel blue) ---
    parts.append(f'<polygon points="{band_polygon}" fill="#1E40AF" '
                 'fill-opacity="0.10" stroke="#1E40AF" stroke-opacity="0.25" '
                 'stroke-width="0.75" clip-path="url(#plotArea2a)"/>')

    # Regression dashed line
    rx0, ry0 = 0.0, intercept
    rx1, ry1 = GB_MAX, slope * GB_MAX + intercept
    parts.append(f'<line x1="{sx(rx0):.1f}" y1="{sy(ry0):.1f}" '
                 f'x2="{sx(rx1):.1f}" y2="{sy(ry1):.1f}" '
                 'stroke="#1A3590" stroke-width="2.5" stroke-dasharray="6,3"/>')

    # --- Error bars (vertical) for datasets with rep data ---
    for r in rows:
        if r['std_min'] is None: continue
        x = sx(r['fastq_GB'])
        y_mid = sy(r['total_min'])
        y_lo  = sy(r['total_min'] - r['std_min'])
        y_hi  = sy(r['total_min'] + r['std_min'])
        # vertical line
        parts.append(f'<line x1="{x:.1f}" y1="{y_lo:.1f}" x2="{x:.1f}" '
                     f'y2="{y_hi:.1f}" stroke="#1A3590" stroke-width="1.4" '
                     'stroke-opacity="0.85"/>')
        # cap top/bot
        for ye in (y_lo, y_hi):
            parts.append(f'<line x1="{x-4:.1f}" y1="{ye:.1f}" x2="{x+4:.1f}" '
                         f'y2="{ye:.1f}" stroke="#1A3590" stroke-width="1.4" '
                         'stroke-opacity="0.85"/>')

    # --- Data points: only draw labels for non-cluster points here ---
    for r in rows:
        x, y = sx(r['fastq_GB']), sy(r['total_min'])
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6.5" '
                     'fill="#1E40AF" stroke="white" stroke-width="2"/>')
        if r['name'] in SMALL_CLUSTER:
            continue  # labels go into callout box below
        suffix = f' (n={r["n_trials"]})' if r['n_trials'] else ''
        # Custom label placement per-dataset to avoid overlap
        LABEL_POS = {
            # name: (dx, dy, anchor)
            'pisa':       (-9, -8, 'end'),
            'liberty':    ( 9, -8, 'start'),
            'goldengate': ( 9,  4, 'start'),
            'venus':      (-9, 16, 'end'),
            'tripitaka':  (-9, -2, 'end'),
        }
        dx, dy, anchor = LABEL_POS.get(r['name'], (9, -3, 'start'))
        parts.append(f'<text x="{x+dx:.1f}" y="{y+dy:.1f}" font-size="11.5" '
                     f'fill="#1A1A1A" font-weight="600" text-anchor="{anchor}">'
                     f'{r["name"]}</text>')
        parts.append(f'<text x="{x+dx:.1f}" y="{y+dy+13:.1f}" font-size="10" '
                     f'fill="#4B5563" text-anchor="{anchor}">'
                     f'{r["total_min"]:.1f} min{suffix}</text>')

    # --- Callout box for the small-FASTQ cluster (Jikji1-7 + seokga) ---
    cluster_x = sx(1.6)  # near centroid of cluster
    cluster_y = sy(2.2)
    parts.append(f'<path d="M {cluster_x:.1f} {cluster_y:.1f} Q 145 240 165 200" '
                 'fill="none" stroke="#9CA3AF" stroke-width="1.2" stroke-dasharray="4,3"/>')
    parts.append(f'<circle cx="{cluster_x:.1f}" cy="{cluster_y:.1f}" r="2" fill="#9CA3AF"/>')
    parts.append('<rect x="135" y="55" width="250" height="155" rx="5" fill="white" '
                 'stroke="#1A3590" stroke-width="1.35" fill-opacity="0.96"/>')
    parts.append('<text x="148" y="74" font-size="13" fill="#1A3590" '
                 'font-weight="bold">Small-file cluster (&lt; 3 GB)</text>')
    parts.append('<text x="148" y="89" font-size="10.5" fill="#4B5563">'
                 '8 samples · time in min · (n) = rep count</text>')
    parts.append('<line x1="148" y1="95" x2="372" y2="95" stroke="#E5E7EB" stroke-width="1"/>')
    cluster_rows = [r for r in rows if r['name'] in SMALL_CLUSTER]
    cluster_rows.sort(key=lambda r: r['name'])
    yi = 112
    for i, r in enumerate(cluster_rows):
        col = i % 2
        cx_lbl = 148 if col == 0 else 268
        cy_lbl = yi + (i // 2) * 18
        suffix = f' (n={r["n_trials"]})' if r['n_trials'] else ''
        parts.append(f'<text x="{cx_lbl}" y="{cy_lbl}" font-size="11.5" fill="#1A1A1A">'
                     f'<tspan font-weight="600">{r["name"]}</tspan>  '
                     f'{r["total_min"]:.2f}{suffix}</text>')

    # Regression annotation box — mid-plot area below tripitaka and above venus
    parts.append('<rect x="408" y="170" width="180" height="80" rx="4" '
                 'fill="white" stroke="#1A3590" stroke-width="1.35" '
                 'fill-opacity="0.94"/>')
    parts.append('<text x="418" y="187" font-size="13.5" fill="#1A3590" '
                 'font-weight="bold">Linear Regression</text>')
    parts.append(f'<text x="418" y="203" font-size="12" fill="#1A1A1A">'
                 f'R² = {r2:.4f}</text>')
    parts.append(f'<text x="418" y="218" font-size="12" fill="#1A1A1A">'
                 f'Slope = {slope:.3f} min/GB</text>')
    parts.append(f'<text x="418" y="233" font-size="10.5" fill="#4B5563">'
                 f's_e = {s_e:.2f} min (n={n})</text>')
    parts.append(f'<text x="418" y="246" font-size="10.5" fill="#4B5563">'
                 f't<tspan dy="2" font-size="8">.025,{n-2}</tspan>'
                 f'<tspan dy="-2"> = {t_critical(n-2):.3f}</tspan></text>')

    parts.append('</svg>')
    return '\n'.join(parts)


# =============================================================================
# Per-step n=10/97 repetition data (downsampled, ratio = 0.188, df = 5.33)
# =============================================================================
PER_STEP_KEYS = ['t_Evert','t_PrimerSort','t_KeyFind','t_Consensus',
                 't_RS','t_NullFind','t_XOR']
# Mapping between per-step n=10 JSON keys and metric-JSON timecell keys
STEP_TO_TIMECELL = {
    't_Evert':      'Evert',
    't_PrimerSort': 'Primer Sorting',
    't_KeyFind':    'Key Finding',
    't_Consensus':  'Index Sorting',
    't_RS':         'RS Decoding',
    't_NullFind':   'Null Finding',
    't_XOR':        'XOR Recovery',
}


def load_per_step_stats() -> Dict[str, Dict[str, Dict]]:
    """Per-step CV (relative std). Returns {ds: {timecell_key: cv}}."""
    out = {}
    for ds, path in REP_SOURCES.items():
        if not path.exists(): continue
        with open(path) as f: m = json.load(f)
        s = m['summary']
        d = {}
        for k in PER_STEP_KEYS:
            if k in s and s[k].get('mean', 0) > 0:
                d[STEP_TO_TIMECELL[k]] = {
                    'mean_down_s': s[k]['mean'],
                    'std_down_s':  s[k]['std'],
                    'cv':          s[k]['std'] / s[k]['mean'],
                    'n':           s[k]['n'],
                }
        out[ds] = d
    return out


def load_full_per_step() -> Dict[str, Dict[str, float]]:
    """Per-step full-FASTQ timings (seconds, n=1 from v8_10/11)."""
    out = {}
    for ds in DATASET_ORDER:
        v11 = RES_DIR / ds / f'{ds}_metrics_v8_11.json'
        v10 = RES_DIR / ds / f'{ds}_metrics_v8_10.json'
        path = v11 if v11.exists() else (v10 if v10.exists() else None)
        if path is None: continue
        with open(path) as f: m = json.load(f)
        tc = m.get('timecell', {}) or {}
        kf = m.get('keyfind_subtiming', {}) or {}
        out[ds] = {
            'fastq_GB': m.get('fastq_size_bytes', 0) / 1e9,
            'total_s': sum(v for v in tc.values() if isinstance(v,(int,float))),
            'Evert':         tc.get('Evert', 0) or 0,
            'Primer Sorting':tc.get('Primer Sorting', 0) or 0,
            'Key Finding':   tc.get('Key Finding', 0) or 0,
            'Index Sorting': tc.get('Index Sorting', 0) or 0,
            'RS Decoding':   tc.get('RS Decoding', 0) or 0,
            'Null Finding':  tc.get('Null Finding', 0) or 0,
            'XOR Recovery':  tc.get('XOR Recovery', 0) or 0,
            'STL Decode':    tc.get('STL Decode', 0) or 0,
            'kis_mask_s':    kf.get('kis_mask_s', 0) or 0,
            'oligo_len':     m.get('oligo_len'),
        }
    return out


FIVE_REPS = ['seokga', 'Jikji6', 'goldengate', 'liberty', 'tripitaka']
PALETTE = {
    'Primer Sorting': '#1E40AF',
    'Index Sorting':  '#E85520',
    'RS+XOR':         '#20B8A6',
    'KIS Detection':  '#E53935',
    'Other':          '#6B7280',
}


# =============================================================================
# Figure 2b v2: stacked composition + step-level error caps
# =============================================================================
def render_2b(per_step: Dict, per_step_reps: Dict) -> str:
    rep = [{'name': n, **per_step[n]} for n in FIVE_REPS if n in per_step]
    rep.sort(key=lambda r: r['oligo_len'] or 0)

    PX, PR = 90, 600
    PY_TOP, PY_BOT = 60, 355
    H = PY_BOT - PY_TOP

    bar_w = 70
    total_w = (PR - PX) - 30
    n_bars = len(rep)
    spacing = (total_w - bar_w * n_bars) / (n_bars + 1)
    start_x = PX + spacing

    parts = []
    parts.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 660 430" '
                 'font-family="Arial, Helvetica, sans-serif">')
    parts.append('<rect width="660" height="430" fill="white"/>')
    parts.append('<text x="340" y="22" font-size="16.5" fill="#1A1A1A" '
                 'text-anchor="middle" font-weight="600">'
                 'Processing Time Composition by Decoding Step</text>')
    parts.append('<text x="340" y="38" font-size="11.5" fill="#4B5563" '
                 'text-anchor="middle">'
                 'Python serial baseline · 5 representative datasets · '
                 'error caps = ± step SD (n≥10 trials, CV-scaled)</text>')

    # Axes
    parts.append(f'<line x1="{PX}" y1="{PY_TOP}" x2="{PX}" y2="{PY_BOT}" '
                 'stroke="#1A1A1A" stroke-width="1.69"/>')
    parts.append(f'<line x1="{PX}" y1="{PY_BOT}" x2="{PR}" y2="{PY_BOT}" '
                 'stroke="#1A1A1A" stroke-width="1.69"/>')
    parts.append('<text x="30" y="200" font-size="17.3" fill="#1A1A1A" '
                 'text-anchor="middle" transform="rotate(-90, 30, 200)">'
                 'Proportion of Total Time (%)</text>')

    for v in (0, 20, 40, 60, 80, 100):
        y = PY_BOT - v / 100 * H
        if v > 0:
            parts.append(f'<line x1="{PX}" y1="{y:.1f}" x2="{PR}" y2="{y:.1f}" '
                         'stroke="#EBEBEB" stroke-width="0.7"/>')
        parts.append(f'<text x="{PX-10}" y="{y+4:.1f}" font-size="14.7" '
                     f'fill="#4B5563" text-anchor="end">{v}</text>')

    # Stacked bars + error caps
    for i, r in enumerate(rep):
        bx = start_x + i * (bar_w + spacing)
        total = r['total_s'] or 1.0
        cats = [
            ('Primer Sorting', r['Primer Sorting'], PALETTE['Primer Sorting']),
            ('KIS Detection',  r['Key Finding'],    PALETTE['KIS Detection']),
            ('Index Sorting',  r['Index Sorting'],  PALETTE['Index Sorting']),
            ('RS+XOR',         r['RS Decoding'] + r['XOR Recovery'], PALETTE['RS+XOR']),
            ('Other',          r['Evert'] + r['Null Finding'] + r['STL Decode'],
                               PALETTE['Other']),
        ]
        cum = 0.0
        rep_data = per_step_reps.get(r['name'], {})
        for label, sec, color in cats:
            pct = sec / total * 100
            seg_h = pct / 100 * H
            seg_y = PY_BOT - cum - seg_h
            parts.append(f'<rect x="{bx:.1f}" y="{seg_y:.1f}" '
                         f'width="{bar_w}" height="{seg_h:.1f}" fill="{color}"/>')
            if pct >= 6:
                ly = seg_y + seg_h / 2 + 4
                parts.append(f'<text x="{bx + bar_w/2:.1f}" y="{ly:.1f}" '
                             f'font-size="11" fill="white" text-anchor="middle" '
                             f'font-weight="600">{pct:.0f}%</text>')

            # Error cap on segment TOP (only for steps with rep data)
            cv = None
            if label == 'Primer Sorting': cv = rep_data.get('Primer Sorting',{}).get('cv')
            elif label == 'KIS Detection': cv = rep_data.get('Key Finding',{}).get('cv')
            elif label == 'Index Sorting': cv = rep_data.get('Index Sorting',{}).get('cv')
            elif label == 'RS+XOR':
                # combined std: SD = sqrt(SD_RS^2 + SD_XOR^2) (independent)
                rs = rep_data.get('RS Decoding', {})
                xr = rep_data.get('XOR Recovery', {})
                if rs and xr:
                    rs_std_full = rs['cv'] * r['RS Decoding']
                    xr_std_full = xr['cv'] * r['XOR Recovery']
                    sd_full = math.sqrt(rs_std_full**2 + xr_std_full**2)
                    if sec > 0: cv = sd_full / sec
            if cv:
                std_pct = cv * pct
                err_h = std_pct / 100 * H
                cap_y_lo = seg_y - err_h
                cap_y_hi = seg_y + err_h
                cap_x = bx + bar_w / 2
                # Error bar line + caps at segment TOP
                parts.append(f'<line x1="{cap_x}" y1="{cap_y_lo:.1f}" '
                             f'x2="{cap_x}" y2="{cap_y_hi:.1f}" '
                             'stroke="#1A1A1A" stroke-width="1.1" '
                             'stroke-opacity="0.65"/>')
                parts.append(f'<line x1="{cap_x-4}" y1="{cap_y_lo:.1f}" '
                             f'x2="{cap_x+4}" y2="{cap_y_lo:.1f}" '
                             'stroke="#1A1A1A" stroke-width="1.1" '
                             'stroke-opacity="0.65"/>')
                parts.append(f'<line x1="{cap_x-4}" y1="{cap_y_hi:.1f}" '
                             f'x2="{cap_x+4}" y2="{cap_y_hi:.1f}" '
                             'stroke="#1A1A1A" stroke-width="1.1" '
                             'stroke-opacity="0.65"/>')
            cum += seg_h

        # Labels under bar
        has_rep = r['name'] in per_step_reps and per_step_reps[r['name']]
        n_lbl = '(n=10)' if has_rep and r['name'] != 'Jikji6' else ('(n=97)' if r['name']=='Jikji6' else '(n=1)')
        parts.append(f'<text x="{bx + bar_w/2:.1f}" y="{PY_BOT+18:.1f}" '
                     f'font-size="13" fill="#1A1A1A" text-anchor="middle" '
                     f'font-weight="600">{r["name"]}</text>')
        parts.append(f'<text x="{bx + bar_w/2:.1f}" y="{PY_BOT+32:.1f}" '
                     f'font-size="11" fill="#4B5563" text-anchor="middle">'
                     f'{r["fastq_GB"]:.2f} GB</text>')
        parts.append(f'<text x="{bx + bar_w/2:.1f}" y="{PY_BOT+45:.1f}" '
                     f'font-size="11" fill="#4B5563" text-anchor="middle">'
                     f'{r["total_s"]/60:.1f} min</text>')
        parts.append(f'<text x="{bx + bar_w/2:.1f}" y="{PY_BOT+58:.1f}" '
                     f'font-size="9.5" fill="#9CA3AF" text-anchor="middle">{n_lbl}</text>')

    # Legend — two rows to avoid overflow
    legend_items = [
        ('Primer Sorting',  PALETTE['Primer Sorting']),
        ('KIS (Key Find)',  PALETTE['KIS Detection']),
        ('Index (Consensus)', PALETTE['Index Sorting']),
        ('RS + XOR',        PALETTE['RS+XOR']),
        ('Other',           PALETTE['Other']),
    ]
    lx_start = 90
    for i, (lbl, c) in enumerate(legend_items):
        ly = 422
        offset = i * 110
        parts.append(f'<rect x="{lx_start + offset}" y="{ly-9}" width="10" height="10" '
                     f'rx="1" fill="{c}"/>')
        parts.append(f'<text x="{lx_start + offset + 14}" y="{ly}" font-size="10.5" '
                     f'fill="#1A1A1A">{lbl}</text>')

    parts.append('</svg>')
    return '\n'.join(parts)


# =============================================================================
# Figure 2c v2: per-step time vs FASTQ size with error bars
# =============================================================================
def render_2c(per_step: Dict, per_step_reps: Dict) -> str:
    PX, PR = 85, 595
    PY_TOP, PY_BOT = 60, 360
    GB_MAX = 24
    SEC_MAX = 5500

    def sx(gb): return PX + gb / GB_MAX * (PR - PX)
    def sy(s):  return PY_BOT - s / SEC_MAX * (PY_BOT - PY_TOP)

    rows = sorted(per_step.values(), key=lambda r: r['fastq_GB'])
    name_by_gb = sorted(per_step.keys(), key=lambda n: per_step[n]['fastq_GB'])

    parts = []
    parts.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 660 430" '
                 'font-family="Arial, Helvetica, sans-serif">')
    parts.append('<rect width="660" height="430" fill="white"/>')
    parts.append('<text x="340" y="22" font-size="16.5" fill="#1A1A1A" '
                 'text-anchor="middle" font-weight="600">'
                 'Per-Step Processing Time vs. FASTQ File Size</text>')
    parts.append('<text x="340" y="38" font-size="11.5" fill="#4B5563" '
                 'text-anchor="middle">'
                 'Python serial · n=13 datasets · error bars = ± step SD on 6 with n≥10 trials</text>')
    parts.append(f'<line x1="{PX}" y1="{PY_TOP}" x2="{PX}" y2="{PY_BOT}" '
                 'stroke="#1A1A1A" stroke-width="1.69"/>')
    parts.append(f'<line x1="{PX}" y1="{PY_BOT}" x2="{PR}" y2="{PY_BOT}" '
                 'stroke="#1A1A1A" stroke-width="1.69"/>')
    parts.append('<text x="28" y="210" font-size="17.3" fill="#1A1A1A" '
                 'text-anchor="middle" transform="rotate(-90, 28, 210)">'
                 'Absolute Processing Time (sec)</text>')
    parts.append('<text x="340" y="400" font-size="17.3" fill="#1A1A1A" '
                 'text-anchor="middle">FASTQ File Size (GB)</text>')

    # Y axis
    for v in (0, 1000, 2000, 3000, 4000, 5000):
        y = sy(v)
        if v > 0:
            parts.append(f'<line x1="{PX}" y1="{y:.1f}" x2="{PR}" y2="{y:.1f}" '
                         'stroke="#EBEBEB" stroke-width="0.7"/>')
        parts.append(f'<text x="{PX-8}" y="{y+4:.1f}" font-size="14.7" '
                     f'fill="#4B5563" text-anchor="end">{v}</text>')

    # X axis
    for v in (0, 5, 10, 15, 20):
        x = sx(v)
        parts.append(f'<text x="{x:.1f}" y="{PY_BOT+17}" font-size="14.7" '
                     f'fill="#4B5563" text-anchor="middle">{v}</text>')
        if v > 0:
            parts.append(f'<line x1="{x:.1f}" y1="{PY_TOP}" x2="{x:.1f}" '
                         f'y2="{PY_BOT}" stroke="#EBEBEB" stroke-width="0.7"/>')

    series = [
        ('Primer Sort',  PALETTE['Primer Sorting'], 'Primer Sorting', 'Primer Sorting'),
        ('Index Cluster',PALETTE['Index Sorting'],  'Index Sorting',  'Index Sorting'),
        ('RS+XOR',       PALETTE['RS+XOR'],          None,             None),
        ('KIS-mask',     PALETTE['KIS Detection'],   'kis_mask_s',     'Key Finding'),
    ]
    r2_values = {}
    for label, color, attr, rep_key in series:
        pts, ds_names = [], []
        for n in name_by_gb:
            r = per_step[n]
            if attr is None:  # RS+XOR
                y_val = r['RS Decoding'] + r['XOR Recovery']
            else:
                y_val = r[attr]
            pts.append((r['fastq_GB'], y_val))
            ds_names.append(n)
        # Polyline
        polyline = ' '.join(f'{sx(x):.1f},{sy(y):.1f}' for x, y in pts)
        parts.append(f'<polyline points="{polyline}" fill="none" '
                     f'stroke="{color}" stroke-width="2.5" stroke-opacity="0.9"/>')
        # Points + error bars
        for (x_, y_), ds in zip(pts, ds_names):
            # Error bar if rep data exists for this step
            rep_data = per_step_reps.get(ds, {})
            cv = None
            if attr == 'Primer Sorting': cv = rep_data.get('Primer Sorting',{}).get('cv')
            elif attr == 'Index Sorting': cv = rep_data.get('Index Sorting',{}).get('cv')
            elif attr is None:  # RS+XOR
                rs = rep_data.get('RS Decoding', {})
                xr = rep_data.get('XOR Recovery', {})
                r_full = per_step[ds]
                if rs and xr and y_ > 0:
                    sd = math.sqrt((rs['cv'] * r_full['RS Decoding'])**2
                                 + (xr['cv'] * r_full['XOR Recovery'])**2)
                    cv = sd / y_
            elif attr == 'kis_mask_s':
                # use Key Finding CV as proxy for KIS sub-step
                cv = rep_data.get('Key Finding',{}).get('cv')
            if cv:
                err = cv * y_
                parts.append(f'<line x1="{sx(x_):.1f}" y1="{sy(y_-err):.1f}" '
                             f'x2="{sx(x_):.1f}" y2="{sy(y_+err):.1f}" '
                             f'stroke="{color}" stroke-width="1.1" '
                             f'stroke-opacity="0.65"/>')
                parts.append(f'<line x1="{sx(x_)-3:.1f}" y1="{sy(y_-err):.1f}" '
                             f'x2="{sx(x_)+3:.1f}" y2="{sy(y_-err):.1f}" '
                             f'stroke="{color}" stroke-width="1.1" '
                             f'stroke-opacity="0.65"/>')
                parts.append(f'<line x1="{sx(x_)-3:.1f}" y1="{sy(y_+err):.1f}" '
                             f'x2="{sx(x_)+3:.1f}" y2="{sy(y_+err):.1f}" '
                             f'stroke="{color}" stroke-width="1.1" '
                             f'stroke-opacity="0.65"/>')
            parts.append(f'<circle cx="{sx(x_):.1f}" cy="{sy(y_):.1f}" r="4" '
                         f'fill="{color}" stroke="white" stroke-width="1.2"/>')
        # R²
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        _, _, r2, _, _ = linreg(xs, ys)
        r2_values[label] = r2
        # Endpoint label — anchor to end so labels stay inside the canvas
        parts.append(f'<text x="{sx(pts[-1][0])-6:.1f}" '
                     f'y="{sy(pts[-1][1])-6:.1f}" font-size="11" '
                     f'fill="{color}" font-weight="600" text-anchor="end">'
                     f'{label}</text>')

    # R² box
    parts.append('<rect x="100" y="70" width="200" height="80" rx="4" '
                 'fill="white" stroke="#4B5563" stroke-width="0.99" '
                 'fill-opacity="0.93"/>')
    parts.append('<text x="110" y="87" font-size="13" fill="#1A1A1A" '
                 'font-weight="bold">Per-step R² (linear fit)</text>')
    yi = 102
    for lbl, c, _, _ in series:
        parts.append(f'<rect x="115" y="{yi-7}" width="8" height="8" rx="1" '
                     f'fill="{c}"/>')
        parts.append(f'<text x="128" y="{yi}" font-size="12" fill="#1A1A1A">'
                     f'{lbl}: R² = {r2_values[lbl]:.3f}</text>')
        yi += 13

    parts.append('</svg>')
    return '\n'.join(parts)


def main():
    rows = load_full_measurements()
    reps = load_repetition_stats()
    per_step = load_full_per_step()
    per_step_reps = load_per_step_stats()

    print(f'[v2] Loaded {len(rows)} datasets; {len(reps)} with total-time reps; '
          f'{len(per_step_reps)} with per-step reps.')

    out_a = OUT_DIR / 'Figure_2a_Time_Regression_filled_v2.svg'
    out_b = OUT_DIR / 'Figure_2b_Time_Composition_filled_v2.svg'
    out_c = OUT_DIR / 'Figure_2c_PerStep_Scaling_filled_v2.svg'

    out_a.write_text(render_2a(rows, reps), encoding='utf-8')
    out_b.write_text(render_2b(per_step, per_step_reps), encoding='utf-8')
    out_c.write_text(render_2c(per_step, per_step_reps), encoding='utf-8')

    print(f'\nWrote:\n  {out_a.name}\n  {out_b.name}\n  {out_c.name}')


if __name__ == '__main__':
    main()
