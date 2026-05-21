#!/usr/bin/env python3
# =============================================================================
# Figure 3 — N=100 Random-Downsampling Statistical Replication (Jikji6)
# =============================================================================
# Builds the statistical evidence behind Figure 3a / 3b / 3c claims:
#   • ×5.33 downsampling (4.57M → ~857K reads) preserves 100% reconstruction.
#   • Decoding time falls from ~34.5 s (full) to ~2.47 s (downsampled).
#   • Per-trial RS-fail count + XOR recovery resilience.
#
# How it differs from run_downsampling_experiment.py:
#   • Single operating point (configurable RATIO, default 1/5.33).
#   • N independent trials with varied seeds.
#   • Master pool loaded ONCE (full FASTQ R1+R2) → reused across trials.
#   • Per-stage timings + raw coverage (reads/OligoLen) reported, matching the
#     manuscript's "×504 / 2.47 s" claim convention.
#   • Excel checkpointing after every trial — safe to interrupt.
#   • Per-trial Depth array dump for ONE representative trial → feeds Fig 3b.
#
# Usage:
#   python3 Figure3/python/run_downsampling_n100.py
#   python3 Figure3/python/run_downsampling_n100.py --n 100 --ratio 0.18762
#   python3 Figure3/python/run_downsampling_n100.py --n 5  # pilot
# =============================================================================

import argparse
import json
import os
import random
import sys
import time
from collections import Counter

import numpy as np

# Locate Final pipeline + reuse helpers from the single-shot script ----------
HERE = os.path.dirname(os.path.abspath(__file__))
FIG1_PYTHON = os.path.normpath(os.path.join(HERE, '..', '..', 'Figure1', 'python'))
sys.path.insert(0, FIG1_PYTHON)
sys.path.insert(0, HERE)

from run_downsampling_experiment import (  # noqa: E402
    R1_PATH, R2_PATH, FILENAME,
    load_pool, run_final_pipeline,
)


# =============================================================================
# DATASET PATH RESOLUTION  (mirrors dna_decode_baseline_260507_Final.main)
# =============================================================================

DATA_ROOT = '/Users/young/Desktop/DNA ink3/3D culture heritage seq data'

DATASET_PATHS = {
    # Jikji series — paired Jikji{N}_1.fastq / Jikji{N}_2.fastq
    **{f'Jikji{i}': (
        os.path.join(DATA_ROOT, 'Jikji', f'Jikji{i}_1.fastq'),
        os.path.join(DATA_ROOT, 'Jikji', f'Jikji{i}_2.fastq'),
    ) for i in range(1, 8)},
    'tripitaka': (
        os.path.join(DATA_ROOT, 'Tripitaka_extracted', 'tripitaka_1.fastq'),
        os.path.join(DATA_ROOT, 'Tripitaka_extracted', 'tripitaka_2.fastq'),
    ),
    'liberty': (
        os.path.join(DATA_ROOT, 'Statue of Liberty', 'liberty_1.fastq'),
        os.path.join(DATA_ROOT, 'Statue of Liberty', 'liberty_2.fastq'),
    ),
    'pisa': (
        os.path.join(DATA_ROOT, 'Tower of Pisa', 'pisa_1.fastq'),
        os.path.join(DATA_ROOT, 'Tower of Pisa', 'pisa_2.fastq'),
    ),
    'seokga': (
        os.path.join(DATA_ROOT, 'seokgatap', 'seokga_1.fastq'),
        os.path.join(DATA_ROOT, 'seokgatap', 'seokga_2.fastq'),
    ),
    'venus': (
        os.path.join(DATA_ROOT, 'Venus Statue', 'V_1.fastq'),
        os.path.join(DATA_ROOT, 'Venus Statue', 'V_2.fastq'),
    ),
    'goldengate': (
        os.path.join(DATA_ROOT, 'Golden gate', 'goldengate_1.fastq'),
        os.path.join(DATA_ROOT, 'Golden gate', 'goldengate_2.fastq'),
    ),
}


def resolve_paths(filename):
    """Return (filename, r1, r2). Accepts J{N} alias for Jikji{N}."""
    if filename in DATASET_PATHS:
        return filename, *DATASET_PATHS[filename]
    if filename in (f'J{i}' for i in range(1, 8)):
        canonical = f'Jikji{filename[1:]}'
        return canonical, *DATASET_PATHS[canonical]
    raise KeyError(f'Unknown dataset: {filename}. '
                   f'Known: {sorted(DATASET_PATHS)}')

# =============================================================================
# CONFIG
# =============================================================================

# 1/5.33 ≈ 0.18762 — the "×5.33 downsampling" operating point in the manuscript.
DEFAULT_RATIO = 1.0 / 5.33
DEFAULT_N_TRIALS = 100
DEFAULT_SEED_BASE = 1000  # trial t uses seed = SEED_BASE + t
DEFAULT_OUTPUT_DIR = os.path.join(HERE, 'Outputs', 'n100_downsampling')


# =============================================================================
# DRIVER
# =============================================================================

def trial_seed(base, idx):
    return int(base) + int(idx)


def deterministic_sample(pool, n_target, seed):
    """random.sample with a fixed seed; full pool if n_target >= len(pool)."""
    if n_target >= len(pool):
        return list(pool)
    rng = random.Random(seed)
    return rng.sample(pool, n_target)


def write_checkpoint(records, meta, output_dir, dump_depth_for_trial=None,
                     depth_array=None):
    """Write Excel + JSON checkpoint. Safe to call after every trial."""
    import pandas as pd  # local import — lighter cold-start

    os.makedirs(output_dir, exist_ok=True)
    fname = meta.get('filename', FILENAME)

    df = pd.DataFrame(records)
    excel_path = os.path.join(output_dir, f'n100_{fname}.xlsx')
    json_path = os.path.join(output_dir, f'n100_{fname}.json')

    # --- Trials sheet ---
    trial_cols = [
        'trial_id', 'seed', 'ratio',
        'n_reads_used', 'valid_reads', 'OligoLen',
        'raw_coverage', 'mean_depth', 'min_depth', 'max_depth', 'n_low_depth_le16',
        'rs_uncorrectable', 'rs_corrected',
        'recovery_stage', 'success', 'vertices', 'faces',
        't_Evert', 't_PrimerSort', 't_KeyFind', 't_Consensus',
        't_RS', 't_NullFind', 't_XOR', 't_total',
    ]
    # Some columns may not exist in the first few records; reindex tolerantly.
    df = df.reindex(columns=trial_cols)

    # --- Summary sheet ---
    summary = build_summary(records, meta)
    summary_df = pd.DataFrame(list(summary.items()), columns=['metric', 'value'])

    # --- Stage distribution sheet ---
    stage_counter = Counter([r.get('recovery_stage') or 'FAIL' for r in records])
    stage_df = pd.DataFrame(
        [(s, c, c / len(records)) for s, c in stage_counter.most_common()],
        columns=['recovery_stage', 'count', 'fraction'])

    # --- (V, F) tuples sheet — mesh identity audit ---
    vf_counter = Counter(
        [(r.get('vertices'), r.get('faces')) for r in records if r.get('success')])
    vf_df = pd.DataFrame(
        [(v, f, c) for (v, f), c in vf_counter.most_common()],
        columns=['vertices', 'faces', 'count'])

    with __import__('pandas').ExcelWriter(excel_path, engine='openpyxl') as w:
        df.to_excel(w, sheet_name='trials', index=False)
        summary_df.to_excel(w, sheet_name='summary', index=False)
        stage_df.to_excel(w, sheet_name='recovery_stages', index=False)
        vf_df.to_excel(w, sheet_name='mesh_VF_audit', index=False)

    # --- JSON: full granularity for downstream figure code ---
    payload = {
        'meta': meta,
        'summary': summary,
        'records': records,
    }
    if dump_depth_for_trial is not None and depth_array is not None:
        payload['representative_depth'] = {
            'trial_id': dump_depth_for_trial,
            'depth_per_oligo': [int(x) for x in depth_array.tolist()],
        }
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, default=str)


def build_summary(records, meta):
    """Aggregate stats matching what Figures 3a/3b/3c will consume."""
    if not records:
        return {'n_trials': 0}

    n = len(records)
    successes = [r for r in records if r.get('success')]

    def stat(field, source=None):
        src = source if source is not None else records
        vals = [r.get(field) for r in src if isinstance(r.get(field), (int, float))]
        if not vals:
            return None
        arr = np.asarray(vals, dtype=float)
        return {
            'mean': float(np.mean(arr)),
            'std': float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
            'median': float(np.median(arr)),
            'p5': float(np.percentile(arr, 5)),
            'p95': float(np.percentile(arr, 95)),
            'min': float(np.min(arr)),
            'max': float(np.max(arr)),
            'n': int(len(arr)),
        }

    summary = {
        'n_trials': n,
        'n_success': len(successes),
        'success_rate': len(successes) / n,
        'ratio': meta['ratio'],
        'pool_size': meta['pool_actual'],
        'total_reads_in_fastq': meta['total_reads_in_fastq'],

        'decode_time_total': stat('t_total', successes or records),
        'raw_coverage': stat('raw_coverage', successes or records),
        'mean_depth': stat('mean_depth', successes or records),
        'min_depth': stat('min_depth', successes or records),
        'n_low_depth_le16': stat('n_low_depth_le16', successes or records),
        'rs_uncorrectable': stat('rs_uncorrectable', successes or records),
        'rs_corrected': stat('rs_corrected', successes or records),
        'OligoLen': stat('OligoLen', successes or records),
        'valid_reads': stat('valid_reads', successes or records),

        # per-stage timings (success only — failed trials may have skewed timings)
        't_Evert': stat('t_Evert', successes or records),
        't_PrimerSort': stat('t_PrimerSort', successes or records),
        't_KeyFind': stat('t_KeyFind', successes or records),
        't_Consensus': stat('t_Consensus', successes or records),
        't_RS': stat('t_RS', successes or records),
        't_NullFind': stat('t_NullFind', successes or records),
        't_XOR': stat('t_XOR', successes or records),
    }
    return summary


def flatten_record(trial_id, seed, ratio, n_reads_used, ok, timings, results,
                   oligo_len_for_raw_cov):
    """Flatten run_final_pipeline output into a per-trial row."""
    rec = {
        'trial_id': trial_id,
        'seed': seed,
        'ratio': ratio,
        'n_reads_used': n_reads_used,
        'valid_reads': results.get('valid_reads'),
        'OligoLen': results.get('OligoLen'),
        'raw_coverage': (
            n_reads_used / oligo_len_for_raw_cov
            if oligo_len_for_raw_cov else None
        ),
        'mean_depth': results.get('mean_depth'),
        'min_depth': results.get('min_depth'),
        'max_depth': results.get('max_depth'),
        'rs_uncorrectable': results.get('rs_uncorrectable'),
        'rs_corrected': results.get('rs_corrected'),
        'recovery_stage': results.get('recovery_stage'),
        'success': bool(ok),
        'vertices': results.get('vertices'),
        'faces': results.get('faces'),
        # timings
        't_Evert': timings.get('Evert'),
        't_PrimerSort': timings.get('PrimerSort'),
        't_KeyFind': timings.get('KeyFind'),
        't_Consensus': timings.get('Consensus'),
        't_RS': timings.get('RS'),
        't_NullFind': timings.get('NullFind'),
        't_XOR': timings.get('XOR'),
        't_total': timings.get('total'),
    }
    return rec


def main():
    ap = argparse.ArgumentParser(description='N-trial random-downsample replication.')
    ap.add_argument('--n', type=int, default=DEFAULT_N_TRIALS,
                    help='Number of trials (default 100).')
    ap.add_argument('--ratio', type=float, default=DEFAULT_RATIO,
                    help='Subsample ratio. Default 1/5.33 ≈ 0.18762.')
    ap.add_argument('--seed-base', type=int, default=DEFAULT_SEED_BASE,
                    help='Trial t uses seed=SEED_BASE+t (deterministic & reproducible).')
    ap.add_argument('--pool', type=int, default=None,
                    help='Optional pool size cap (reservoir). Omit = full FASTQ.')
    ap.add_argument('--filename', type=str, default=FILENAME,
                    help='Dataset name (Jikji1..7, tripitaka, liberty, pisa, '
                         'seokga, venus, goldengate, or J{N} alias).')
    ap.add_argument('--output-dir', type=str, default=None,
                    help='Defaults to Outputs/n100_<filename>/.')
    ap.add_argument('--depth-trial', type=int, default=1,
                    help='Trial id whose Depth array is saved for Figure 3b. Default 1.')
    ap.add_argument('--resume', action='store_true',
                    help='Load existing JSON output and only run missing trial ids.')
    args = ap.parse_args()

    fname, r1, r2 = resolve_paths(args.filename)
    out_dir = args.output_dir or os.path.join(HERE, 'Outputs', f'n100_{fname}')

    print('=' * 78)
    print(f'Figure 3 — N={args.n} replication for {fname}')
    print(f'  R1          : {r1}')
    print(f'  R2          : {r2}')
    print(f'  ratio       : {args.ratio:.5f}  (×{1.0/args.ratio:.3f} downsampling)')
    print(f'  pool        : {"FULL" if args.pool is None else f"{args.pool:,} reservoir"}')
    print(f'  seed base   : {args.seed_base}')
    print(f'  output dir  : {out_dir}')
    print('=' * 78)

    pool, total_reads = load_pool(r1, r2, args.pool)
    if not pool:
        print('[error] empty pool. abort.')
        return 1
    n_target = max(1, int(round(total_reads * args.ratio)))
    if n_target > len(pool):
        print(f'[warn] target {n_target:,} > pool {len(pool):,} → using whole pool.')
        n_target = len(pool)
    print(f'\n[plan] each trial samples {n_target:,} reads from pool of {len(pool):,}')
    print(f'[plan] total {args.n} trials; checkpointing after each.\n')

    meta = {
        'filename': fname,
        'pipeline': 'dna_decode_baseline_260507_Final',
        'pool_actual': len(pool),
        'total_reads_in_fastq': total_reads,
        'ratio': args.ratio,
        'downsampling_factor': 1.0 / args.ratio,
        'n_target_per_trial': n_target,
        'seed_base': args.seed_base,
        'n_trials_planned': args.n,
        'started_at': time.strftime('%Y-%m-%d %H:%M:%S'),
    }

    records = []
    representative_depth = None
    representative_trial_id = None

    # --- Resume support: load already-completed trials and skip them. -------
    completed_trial_ids = set()
    if args.resume:
        existing_json = os.path.join(out_dir, f'n100_{fname}.json')
        if os.path.exists(existing_json):
            with open(existing_json, 'r', encoding='utf-8') as f:
                prev = json.load(f)
            for r in prev.get('records', []):
                tid = r.get('trial_id')
                if tid is not None and r.get('seed') == trial_seed(args.seed_base, tid):
                    records.append(r)
                    completed_trial_ids.add(int(tid))
            print(f'[resume] loaded {len(records)} prior trials from {existing_json}')
            print(f'[resume] will run remaining trials: '
                  f'{sorted(set(range(1, args.n+1)) - completed_trial_ids)[:10]}'
                  f'{" ..." if (args.n - len(completed_trial_ids)) > 10 else ""}')
        else:
            print(f'[resume] no prior file at {existing_json}; starting fresh.')

    t_run_start = time.time()
    for t in range(1, args.n + 1):
        if t in completed_trial_ids:
            continue
        seed = trial_seed(args.seed_base, t)

        sample = deterministic_sample(pool, n_target, seed)

        ok, timings, results = run_final_pipeline(sample, fname)

        # Capture Depth array of the chosen representative trial for Fig 3b.
        if t == args.depth_trial and ok and '_depth_array' in results:
            representative_depth = results['_depth_array']
            representative_trial_id = t

        rec = flatten_record(
            trial_id=t, seed=seed, ratio=args.ratio,
            n_reads_used=len(sample), ok=ok, timings=timings, results=results,
            oligo_len_for_raw_cov=results.get('OligoLen'))
        records.append(rec)

        elapsed = time.time() - t_run_start
        n_this_run = t - max(completed_trial_ids, default=0) if completed_trial_ids else t
        rate = elapsed / max(n_this_run, 1)
        eta = rate * (args.n - t)
        # All print fields below are None-safe — a recovery failure (all 3
        # stages exception out) leaves recovery_stage=None and depth fields
        # may also be missing.
        stage = results.get("recovery_stage") or "FAIL"
        v = results.get("vertices") if results.get("vertices") is not None else "?"
        f_ = results.get("faces") if results.get("faces") is not None else "?"
        ttot = timings.get("total")
        ttot_s = f"{ttot}s" if ttot is not None else "?"
        rc_val = rec.get("raw_coverage")
        rc_s = f"{rc_val:.1f}x" if isinstance(rc_val, (int, float)) else "?"
        md = results.get("mean_depth")
        md_s = f"{md}x" if md is not None else "?"
        print(f'  [{t:3d}/{args.n}] seed={seed} reads={len(sample):,} '
              f'ok={ok} stage={stage:>6} '
              f'V/F={v}/{f_} t={ttot_s} cov={rc_s}  depth={md_s}  '
              f'ETA {eta/60:.1f}m')

        # Checkpoint after every trial.
        write_checkpoint(records, meta, out_dir,
                         dump_depth_for_trial=representative_trial_id,
                         depth_array=representative_depth)

    meta['ended_at'] = time.strftime('%Y-%m-%d %H:%M:%S')
    meta['total_runtime_seconds'] = round(time.time() - t_run_start, 1)
    write_checkpoint(records, meta, args.output_dir,
                     dump_depth_for_trial=representative_trial_id,
                     depth_array=representative_depth)

    # Final report ------------------------------------------------------------
    summary = build_summary(records, meta)
    print('\n' + '=' * 78)
    print(f'COMPLETE — {len(records)} trials in {meta["total_runtime_seconds"]:.1f}s')
    print('-' * 78)
    print(f"success rate   : {summary['n_success']}/{summary['n_trials']} "
          f"= {summary['success_rate']*100:.1f}%")
    if summary.get('decode_time_total'):
        s = summary['decode_time_total']
        print(f"decode time    : {s['mean']:.3f} ± {s['std']:.3f} s "
              f"(median {s['median']:.3f}, p5–p95 {s['p5']:.3f}–{s['p95']:.3f})")
    if summary.get('raw_coverage'):
        s = summary['raw_coverage']
        print(f"raw coverage   : ×{s['mean']:.1f} ± {s['std']:.1f}")
    if summary.get('mean_depth'):
        s = summary['mean_depth']
        print(f"realized depth : ×{s['mean']:.1f} ± {s['std']:.1f}")
    if summary.get('n_low_depth_le16'):
        s = summary['n_low_depth_le16']
        print(f"low-depth (≤16): {s['mean']:.1f} ± {s['std']:.1f} oligos/trial")

    stage_dist = Counter([r.get('recovery_stage') or 'FAIL' for r in records])
    print(f"recovery stages: {dict(stage_dist)}")

    print('=' * 78)
    print(f"Excel: {os.path.join(out_dir, f'n100_{fname}.xlsx')}")
    print(f"JSON : {os.path.join(out_dir, f'n100_{fname}.json')}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
