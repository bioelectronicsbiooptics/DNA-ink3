#!/usr/bin/env python3
# =============================================================================
# Figure 3 — Random Downsampling Experiment (Final pipeline edition)
# =============================================================================
# Base: Figure1/python/dna_decode_baseline_260507_Final.py  (production validator)
# Target: Jikji6  (R1 + R2 paired FASTQ)
#
# WHAT THIS DOES
#   1. Load reads from R1+R2 — either ALL or a uniformly-random POOL
#      (reservoir sampling, single streaming pass per file → no front-bias).
#   2. For each target coverage, randomly subsample from the pool
#      (random.sample with fixed seed → reproducible).
#   3. Run the SAME pipeline as the Final validator on the downsampled set,
#      including the 3-stage Recovery cascade (Global → Depth-cut → Local).
#   4. Save per-coverage timings + decode metrics into a single JSON.
#
# OPEN KNOB
#   POOL_SIZE — None means "load everything"; an int caps memory by reservoir-
#   sampling that many reads from the combined R1+R2 stream. The user asked to
#   keep the batch/pool size open since loading all of Jikji6 is heavy and
#   loading from the front would bias the sample.
#
# Usage:
#   python3 Figure3/python/run_downsampling_experiment.py
#   python3 Figure3/python/run_downsampling_experiment.py --pool 200000
#   python3 Figure3/python/run_downsampling_experiment.py --targets 1421,500,100,20
# =============================================================================

import argparse
import json
import os
import random
import sys
import time

import numpy as np

# --- Locate the Final pipeline module ---------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
FIG1_PYTHON = os.path.normpath(os.path.join(HERE, '..', '..', 'Figure1', 'python'))
sys.path.insert(0, FIG1_PYTHON)

from dna_decode_baseline_260507_Final import (  # noqa: E402
    PrimerSelect, ReverseEvert, PrimerSort,
    KeyFind, Keysort, Consensus, Binarization, RSdecoding,
    NullFind, XORrun, nestedXOR, reRSdec, tridecode, trierrXOR,
    _make_triangulation, _dna_to_binary, _de2bi,
    read_fastq_file,
    KEYFIND_LAST_TIMING,
    HAS_TRIMESH,
)

# =============================================================================
# CONFIG
# =============================================================================

FILENAME = 'Jikji6'
DATA_ROOT = '/Users/young/Desktop/DNA ink3/3D culture heritage seq data'
R1_PATH = os.path.join(DATA_ROOT, 'Jikji', 'Jikji6_1.fastq')
R2_PATH = os.path.join(DATA_ROOT, 'Jikji', 'Jikji6_2.fastq')

# Empirical mean depth of the full Jikji6 paired FASTQ (matches reference scripts).
ORIGINAL_COVERAGE = 1421

# Coverage targets — 1421 == full data (no downsampling).
DEFAULT_TARGETS = [1421, 500, 200, 100, 50, 20, 10]

# OPEN KNOB. None = load all reads (peak memory). int = reservoir-sample to
# this many reads in one streaming pass. The largest coverage target you can
# reach is bounded by POOL_SIZE.
DEFAULT_POOL_SIZE = None

OUTPUT_DIR = os.path.join(HERE, 'Outputs', 'downsampling_results')
RNG_SEED = 42


# =============================================================================
# LOADING — full or reservoir
# =============================================================================

def _stream_fastq_seqs(path):
    """Yield sequence strings from a FASTQ file without loading it all.

    Uses the plain-text every-4-lines convention (matches read_fastq_file's
    fallback) so we don't depend on Biopython for the streaming pass.
    Supports .gz transparently.
    """
    import gzip
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        i = 0
        for line in fh:
            if i % 4 == 1:
                yield line.rstrip('\n')
            i += 1


def reservoir_sample_pair(r1, r2, k, seed=RNG_SEED):
    """Uniformly sample k reads across R1 + R2 in a single pass each.

    Returns:
        (sampled list[str], total_reads_seen).
    Uniform over the concatenation of both files — same distribution as
    `read_fastq_file(r1) + read_fastq_file(r2)` then random.sample().
    """
    rng = random.Random(seed)
    reservoir = []
    n = 0
    for path in (r1, r2):
        if not os.path.exists(path):
            print(f"[warn] file not found: {path}")
            continue
        for s in _stream_fastq_seqs(path):
            n += 1
            if len(reservoir) < k:
                reservoir.append(s)
            else:
                j = rng.randint(0, n - 1)
                if j < k:
                    reservoir[j] = s
    return reservoir, n


def load_pool(r1, r2, pool_size):
    """Build the master pool of reads.

    pool_size=None → load everything (uses Final's read_fastq_file).
    pool_size=int  → reservoir-sample to that size in one streaming pass.
    """
    t = time.time()
    if pool_size is None:
        seq = read_fastq_file(r1) + read_fastq_file(r2)
        total = len(seq)
        print(f"[load] full FASTQ R1+R2: {total:,} reads ({round(time.time()-t,1)}s)")
        return seq, total

    if pool_size <= 0:
        raise ValueError(f"POOL_SIZE must be positive or None, got {pool_size}")

    pool, total = reservoir_sample_pair(r1, r2, pool_size, seed=RNG_SEED)
    print(f"[load] reservoir-sampled {len(pool):,} of {total:,} reads "
          f"(pool target={pool_size:,}, {round(time.time()-t,1)}s)")
    return pool, total


# =============================================================================
# DOWNSAMPLING
# =============================================================================

def random_downsample(pool, n_target, seed=RNG_SEED):
    """Random.sample without replacement; deterministic with fixed seed."""
    if n_target >= len(pool):
        return list(pool)
    rng = random.Random(seed)
    return rng.sample(pool, n_target)


def compute_n_target(total_reads, target_cov, original_cov):
    """Reads needed to reach target_cov, scaled from the FULL FASTQ size.

    Note: scales against `total_reads` (the FASTQ count we observed during
    loading), NOT against pool length. That keeps target depth interpretable
    even when we use a reservoir-sampled pool.
    """
    ratio = target_cov / original_cov
    return max(1, int(total_reads * ratio))


# =============================================================================
# PIPELINE — Final validator, modularized so we can call it per-coverage
# =============================================================================

def run_final_pipeline(seq, filename):
    """Run the Final pipeline (260507) on a list of read strings.

    Returns (success: bool, timings: dict, results: dict). Mirrors main() in
    dna_decode_baseline_260507_Final.py — including the cascade where Depth
    and Local stages are attempted even after Global succeeds (later
    successes overwrite). The recovery_stage in `results` reflects the LAST
    stage that produced a usable mesh.
    """
    timings = {}
    results = {}

    try:
        F_primers, R_primers = PrimerSelect(filename)

        # --- Step 1: Evert ---------------------------------------------------
        t = time.time()
        forward_section = [s[:21] if len(s) >= 21 else s for s in seq]
        tevert_seq, reversed_log = ReverseEvert(seq, forward_section, R_primers, F_primers)
        rev_indices = np.where(reversed_log)[0]
        for i, idx in enumerate(rev_indices):
            if tevert_seq[i] is not None:
                seq[idx] = tevert_seq[i]
        timings['Evert'] = round(time.time() - t, 3)

        # --- Step 2: PrimerSort ----------------------------------------------
        t = time.time()
        psSeq = PrimerSort(seq, F_primers, R_primers)
        timings['PrimerSort'] = round(time.time() - t, 3)
        results['valid_reads'] = len(psSeq)

        if len(psSeq) == 0:
            results['error'] = 'no valid reads after PrimerSort'
            return False, timings, results

        # --- Step 3: KeyFind / Keysort ---------------------------------------
        t = time.time()
        k = 156
        n = k + 8
        m = len(_de2bi(n))
        KeyMap, KmatDNA = KeyFind(psSeq)
        OligoLen, idxlen, ncolpo, nbitpo, ncolco, nbitco = Keysort(KeyMap, KmatDNA, k, filename)
        timings['KeyFind'] = round(time.time() - t, 3)
        results.update({
            'OligoLen': int(OligoLen),
            'idxlen': int(idxlen),
            'ncolpo': int(ncolpo), 'nbitpo': int(nbitpo),
            'ncolco': int(ncolco), 'nbitco': int(nbitco),
        })

        # --- Step 4: Consensus (Index Sorting) -------------------------------
        t = time.time()
        MIndex = Consensus(psSeq, idxlen, OligoLen)
        Depth = np.array([MIndex[i].shape[0] if MIndex[i] is not None else 0
                          for i in range(len(MIndex))])
        timings['Consensus'] = round(time.time() - t, 3)
        results['mean_depth'] = round(float(np.mean(Depth)), 1)
        results['min_depth'] = int(np.min(Depth)) if len(Depth) else 0
        results['max_depth'] = int(np.max(Depth)) if len(Depth) else 0
        # Figure 3b needs per-oligo depth distribution. Cheap (≈ a few KB).
        results['_depth_array'] = Depth.copy()
        results['n_low_depth_le16'] = int(np.sum(Depth <= 16))

        # --- Step 5: RS decoding ---------------------------------------------
        t = time.time()
        SeqtoBin = Binarization(MIndex)
        rsECC, RSerrinf = RSdecoding(SeqtoBin, k, m, n)
        rsECC[RSerrinf < 0, :] = 0
        rsECC[RSerrinf > 4, :] = 0
        for a in range(min(OligoLen, len(rsECC))):
            if np.any(rsECC[a, :] > 1):
                rsECC[a, :] = 0
        payload = rsECC[:, :k - idxlen]
        timings['RS'] = round(time.time() - t, 3)
        results['rs_uncorrectable'] = int(np.sum(RSerrinf < 0))
        results['rs_corrected'] = int(np.sum((RSerrinf > 0) & (RSerrinf <= 4)))

        # --- Step 6: NullFind ------------------------------------------------
        t = time.time()
        term = ncolpo * 3 * nbitpo + ncolco * 3 * nbitco
        Nullkey_1 = term + 1 if term % 2 == 1 else term
        s_val = payload.shape[1]
        Nullkey_2 = Nullkey_1 // 2

        Null_count = 0
        nk1 = Nullkey_1
        while nk1 % s_val != 0:
            nk1 += 1
            Null_count += 1
        Null_xor_count = 0
        nk2 = Nullkey_2
        while nk2 % s_val != 0:
            nk2 += 1
            Null_xor_count += 1

        data_Nd, xor_Nd, dNc, xNc = NullFind(seq, Null_count, Null_xor_count, s_val, F_primers)
        bindataN = _dna_to_binary(data_Nd[20:])
        binxorN = _dna_to_binary(xor_Nd[20:])
        timings['NullFind'] = round(time.time() - t, 3)

        # update Depth with null/xor read counts (matches Final main)
        tri_row_idx = payload[:-1, :].size * 2 // (3 * s_val) - 1
        if tri_row_idx < len(Depth) and dNc:
            Depth[tri_row_idx] = sum(len(d) for d in dNc if d)
        if len(Depth) > 1 and xNc:
            Depth[-2] = sum(len(x) for x in xNc if x)

        # --- Step 7: XOR Recovery (3-stage cascade) --------------------------
        t = time.time()
        matsz = k - idxlen
        Temp_mat = payload[:-1, :].copy()

        data_insert_row = Temp_mat.size * 2 // (3 * s_val) - 1
        if data_insert_row < Temp_mat.shape[0]:
            il = min(len(bindataN), Temp_mat.shape[1])
            Temp_mat[data_insert_row, :il] = bindataN[:il]
        if Temp_mat.shape[0] > 0:
            il = min(len(binxorN), Temp_mat.shape[1])
            Temp_mat[-1, :il] = binxorN[:il]

        if len(Depth) > Temp_mat.shape[0]:
            Depth = Depth[:Temp_mat.shape[0]]

        trid = None
        recovery_stage = None
        err_indices = np.array([], dtype=int)

        # 7-1: Global RS
        try:
            ECC_part = Temp_mat[:Temp_mat.size * 2 // (3 * s_val), :]
            ECC_flat = ECC_part.flatten(order='C')
            if Null_count > 0:
                ECC_flat = ECC_flat[:-Null_count]
            tc, tp = tridecode(ECC_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
            trid = _make_triangulation(tc, tp)
            recovery_stage = 'Global'
        except Exception as e:
            results.setdefault('stage_errors', {})['Global'] = str(e)

        # 7-2: Depth-cut + XOR
        Depthcut = 16
        ErrD = np.where(Depth <= Depthcut)[0]
        try:
            err_mask = RSerrinf[ErrD] != 0 if len(ErrD) > 0 else np.array([], dtype=bool)
            err_indices = ErrD[err_mask] if len(err_mask) > 0 else np.array([], dtype=int)
            if len(err_indices) > 0:
                Temp_mat[err_indices, :] = 0
            xtp, _, _, _, XA, XB, XC = XORrun(
                Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)
            ECC_flat = xtp.flatten(order='C')
            if Null_count > 0:
                ECC_flat = ECC_flat[:-Null_count]
            tc, tp = tridecode(ECC_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
            trid = trierrXOR(Temp_mat, tc, tp, Null_count, Null_xor_count,
                             ncolpo, nbitpo, ncolco, nbitco)
            recovery_stage = 'Depth'
        except Exception as e:
            results.setdefault('stage_errors', {})['Depth'] = str(e)

        # 7-3: Local RS (nestedXOR + reRSdec)
        try:
            if len(err_indices) > 0:
                Temp_mat[err_indices, :] = 0
            xtp, _, tri_p, xor_p, XA, XB, XC = XORrun(
                Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)
            if len(ErrD) > 0:
                XORMIndex, Del_xor = nestedXOR(
                    MIndex, xtp, ErrD, RSerrinf, XA, XB, XC,
                    payload, tri_p, xor_p, matsz, Null_count)
                Temp_mat = reRSdec(XORMIndex, Temp_mat, k, m, n, Del_xor, matsz)
                _, missing_data, _, _, _, _, _ = XORrun(
                    Temp_mat, s_val, Null_count, Null_xor_count,
                    ncolpo, ncolco, nbitpo, nbitco)
                XORmData = [MIndex[i] for i in missing_data if i < len(MIndex)]
                Temp_mat = reRSdec(XORmData, Temp_mat, k, m, n, missing_data, matsz)
            xtp_final, _, _, _, _, _, _ = XORrun(
                Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)
            ECC_flat = xtp_final.flatten(order='C')
            if Null_count > 0:
                ECC_flat = ECC_flat[:-Null_count]
            tc, tp = tridecode(ECC_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
            trid = trierrXOR(Temp_mat, tc, tp, Null_count, Null_xor_count,
                             ncolpo, nbitpo, ncolco, nbitco)
            recovery_stage = 'Local'
        except Exception as e:
            results.setdefault('stage_errors', {})['Local'] = str(e)

        timings['XOR'] = round(time.time() - t, 3)

        success = trid is not None
        results['success'] = success
        results['recovery_stage'] = recovery_stage
        if trid is not None:
            if hasattr(trid, 'vertices'):
                results['vertices'] = int(len(trid.vertices))
                results['faces'] = int(len(trid.faces))
            elif isinstance(trid, dict):
                results['vertices'] = len(trid.get('vertices', []))
                results['faces'] = len(trid.get('faces', []))

        timings['total'] = round(sum(v for v in timings.values() if isinstance(v, (int, float))), 3)
        results['keyfind_subtiming'] = dict(KEYFIND_LAST_TIMING)
        return success, timings, results

    except Exception as e:
        results['error'] = f"{type(e).__name__}: {e}"
        results['success'] = False
        timings['total'] = round(sum(v for v in timings.values() if isinstance(v, (int, float))), 3)
        return False, timings, results


# =============================================================================
# DRIVER
# =============================================================================

def parse_targets(s):
    return [int(x) for x in s.split(',') if x.strip()]


def main():
    ap = argparse.ArgumentParser(description='Random downsampling experiment (Final pipeline).')
    ap.add_argument('--pool', type=int, default=DEFAULT_POOL_SIZE,
                    help='Reservoir pool size. Omit (or set <=0) to load ALL reads.')
    ap.add_argument('--targets', type=parse_targets, default=DEFAULT_TARGETS,
                    help='Comma-separated coverage targets (e.g. 1421,500,100,20).')
    ap.add_argument('--seed', type=int, default=RNG_SEED,
                    help='Random seed for reservoir + per-coverage sampling.')
    ap.add_argument('--output-dir', type=str, default=OUTPUT_DIR)
    args = ap.parse_args()

    pool_size = args.pool if (args.pool is None or args.pool > 0) else None
    os.makedirs(args.output_dir, exist_ok=True)

    print('=' * 72)
    print(f'Random Downsampling Experiment — {FILENAME}')
    print(f'  pipeline : Figure1/python/dna_decode_baseline_260507_Final.py')
    print(f'  R1       : {R1_PATH}')
    print(f'  R2       : {R2_PATH}')
    print(f'  pool     : {"ALL" if pool_size is None else f"{pool_size:,} (reservoir)"}')
    print(f'  targets  : {args.targets}')
    print(f'  baseline : {ORIGINAL_COVERAGE}x mean depth (full FASTQ)')
    print(f'  seed     : {args.seed}')
    print('=' * 72)

    # Build the master pool ----------------------------------------------------
    pool, total_reads = load_pool(R1_PATH, R2_PATH, pool_size)
    if not pool:
        print('[error] pool is empty. Aborting.')
        return 1

    all_results = []

    for target_cov in args.targets:
        n_full = compute_n_target(total_reads, target_cov, ORIGINAL_COVERAGE)
        n_use = min(n_full, len(pool))
        clipped = n_full > len(pool)

        print(f'\n{"-" * 60}')
        print(f'Coverage target: {target_cov}x  '
              f'(needs {n_full:,} reads of {total_reads:,} total)')
        if clipped:
            print(f'  [warn] pool has only {len(pool):,} reads — using all of them. '
                  f'Effective coverage will be ≤ '
                  f'{round(len(pool) / total_reads * ORIGINAL_COVERAGE, 1)}x.')
        ds_seq = random_downsample(pool, n_use, seed=args.seed)
        print(f'  sampled: {len(ds_seq):,} reads '
              f'({len(ds_seq) / total_reads * 100:.2f}% of FASTQ)')

        ok, timings, results = run_final_pipeline(ds_seq, FILENAME)

        record = {
            'target_coverage': target_cov,
            'n_reads_used': len(ds_seq),
            'n_reads_target_full': n_full,
            'pool_clipped': clipped,
            'pool_size': len(pool),
            'total_reads_in_fastq': total_reads,
            'timings': timings,
            **results,
        }
        all_results.append(record)

        # Console summary
        print(f'  success      : {ok}')
        if ok:
            print(f"  realized depth: {results.get('mean_depth', '?')}x "
                  f"(min={results.get('min_depth', '?')})")
            print(f"  RS uncorr.   : {results.get('rs_uncorrectable', '?')}")
            print(f"  recovery     : {results.get('recovery_stage', '?')}")
            print(f"  V/F          : {results.get('vertices', '?')}/{results.get('faces', '?')}")
        else:
            err = results.get('error') or results.get('stage_errors') or 'unknown'
            print(f"  error        : {err}")
        print(f"  total time   : {timings.get('total', '?')}s")

    # Save JSON ---------------------------------------------------------------
    out_json = os.path.join(args.output_dir, f'downsampling_experiment_{FILENAME}.json')
    payload = {
        'meta': {
            'filename': FILENAME,
            'pipeline': 'dna_decode_baseline_260507_Final',
            'pool_size': pool_size,
            'pool_actual': len(pool),
            'total_reads_in_fastq': total_reads,
            'original_coverage': ORIGINAL_COVERAGE,
            'targets': args.targets,
            'seed': args.seed,
            'has_trimesh': bool(HAS_TRIMESH),
        },
        'results': all_results,
    }
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2, default=str, ensure_ascii=False)
    print(f'\n[saved] {out_json}')

    # Summary table -----------------------------------------------------------
    print('\n' + '=' * 78)
    print(f"{'Target':>8} {'Reads':>12} {'Depth':>8} {'RSerr':>7} "
          f"{'Stage':>8} {'V/F':>14} {'Time(s)':>10}")
    print('-' * 78)
    for r in all_results:
        cov = r.get('target_coverage', '?')
        reads = r.get('n_reads_used', 0)
        depth = r.get('mean_depth', '?')
        rs_err = r.get('rs_uncorrectable', '?')
        stage = r.get('recovery_stage') or ('FAIL' if not r.get('success') else '?')
        v = r.get('vertices', '?')
        fc = r.get('faces', '?')
        tot = r.get('timings', {}).get('total', '?')
        print(f"{cov:>7}x {reads:>12,} {depth:>8} {rs_err:>7} "
              f"{str(stage):>8} {f'{v}/{fc}':>14} {tot:>10}")
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
