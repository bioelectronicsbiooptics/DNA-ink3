#!/usr/bin/env python3
# =============================================================================
# Figure 4 — Option 1: Front-batch sequential decoding with PRE-COMPUTED KIS
# Created: 2026-05-08
# =============================================================================
# IDEA
#   Sequential (NOT random) batched decoding. The full FASTQ is loaded once.
#   We pre-compute the KIS map (OligoLen, idxlen, ncolpo/nbitpo/ncolco/nbitco)
#   from the ENTIRE FASTQ first — this is the "trick": the decoder is given
#   ground-truth structural metadata, so partial-batch failures can ONLY be
#   blamed on insufficient depth, not on Keysort detection variance.
#
#   We then sweep batch sizes from the front of the FASTQ:
#       seq[:batch_size]  →  decode  →  success? grow batch.
#   The minimum batch where decoding succeeds is the answer.
#
# WHY pre-compute KIS:
#   In our N=100 J6 run, 3 failures all had OligoLen=938 instead of 1018 —
#   Keysort misdetection from a particular random subset. Here, for a CLEAN
#   test of "how few reads to decode?", we eliminate that nuisance variable.
#
# Pipeline reuse:
#   We import functions from the Final pipeline (Figure1/python/
#   dna_decode_baseline_260507_Final.py). We do NOT modify Final — we copy
#   only the orchestration logic from main() into a small `decode_with_kis()`
#   helper here, replacing the KeyFind+Keysort step with the pre-computed
#   values.
#
# Usage:
#   python3 Figure4/python/Option1_front_batch_decode_260508.py
#   python3 Figure4/python/Option1_front_batch_decode_260508.py --batches 50000,200000,500000,858921
#   python3 Figure4/python/Option1_front_batch_decode_260508.py --grow-from 50000 --grow-step 50000
# =============================================================================

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG1_PYTHON = os.path.normpath(os.path.join(HERE, '..', '..', 'Figure1', 'python'))
sys.path.insert(0, FIG1_PYTHON)

# Import correctness primitives from the Final pipeline. We do NOT modify
# them — we just reuse them as building blocks.
from dna_decode_baseline_260507_Final import (  # noqa: E402
    PrimerSelect, ReverseEvert, PrimerSort, KeyFind, Keysort,
    Consensus, Binarization, RSdecoding,
    NullFind, XORrun, nestedXOR, reRSdec, tridecode, trierrXOR,
    _make_triangulation, _dna_to_binary, _de2bi,
    read_fastq_file, HAS_TRIMESH,
)

# =============================================================================
# CONFIG
# =============================================================================

FILENAME = 'Jikji6'
DATA_ROOT = '/Users/young/Desktop/DNA ink3/3D culture heritage seq data'
R1_PATH = os.path.join(DATA_ROOT, 'Jikji', 'Jikji6_1.fastq')
R2_PATH = os.path.join(DATA_ROOT, 'Jikji', 'Jikji6_2.fastq')

OUTPUT_DIR = os.path.join(HERE, 'Outputs', 'option1_front_batch_260508')

# Default batch sweep — from very small to nearly full.
DEFAULT_BATCH_SIZES = [
    50_000, 100_000, 200_000, 400_000, 600_000, 858_921, 1_500_000, 4_578_048
]

# Grow mode: start small, add this many reads each step until success.
DEFAULT_GROW_FROM = 50_000
DEFAULT_GROW_STEP = 50_000
DEFAULT_GROW_MAX = 4_578_048


# =============================================================================
# KIS pre-computation (run ONCE on the full FASTQ)
# =============================================================================

def precompute_kis(seq, filename):
    """Run Evert + PrimerSort + KeyFind + Keysort on the FULL FASTQ once.

    Returns a dict of structural parameters that downstream batches reuse.
    Mirrors steps 1–3 of dna_decode_baseline_260507_Final.main().
    """
    t0 = time.time()
    print('[kis] pre-compute on FULL FASTQ ...')

    F_primers, R_primers = PrimerSelect(filename)

    # Step 1: Evert
    forward_section = [s[:21] if len(s) >= 21 else s for s in seq]
    tevert_seq, reversed_log = ReverseEvert(seq, forward_section, R_primers, F_primers)
    rev_indices = np.where(reversed_log)[0]
    seq_oriented = list(seq)  # don't mutate caller's list reference unexpectedly
    for i, idx in enumerate(rev_indices):
        if tevert_seq[i] is not None:
            seq_oriented[idx] = tevert_seq[i]
    t_evert = round(time.time() - t0, 2)

    # Step 2: PrimerSort (on the WHOLE FASTQ)
    t1 = time.time()
    ps_full = PrimerSort(seq_oriented, F_primers, R_primers)
    t_ps = round(time.time() - t1, 2)

    # Step 3: KeyFind / Keysort
    t2 = time.time()
    k = 156
    n = k + 8  # 164
    m = len(_de2bi(n))
    KeyMap, KmatDNA = KeyFind(ps_full)
    OligoLen, idxlen, ncolpo, nbitpo, ncolco, nbitco = Keysort(KeyMap, KmatDNA, k, filename)
    t_kf = round(time.time() - t2, 2)

    print(f'[kis] reads={len(seq):,}  oriented_evert={t_evert}s  primer_sort={t_ps}s  '
          f'key_find={t_kf}s  → OligoLen={OligoLen}, idxlen={idxlen}')

    return {
        'F_primers': F_primers,
        'R_primers': R_primers,
        'OligoLen': int(OligoLen),
        'idxlen': int(idxlen),
        'ncolpo': int(ncolpo),
        'nbitpo': int(nbitpo),
        'ncolco': int(ncolco),
        'nbitco': int(nbitco),
        'k': k, 'n': n, 'm': m,
        'time_kis_total': round(time.time() - t0, 2),
        'seq_oriented': seq_oriented,  # carry oriented reads for batch use
    }


# =============================================================================
# Decode a batch using the pre-computed KIS
# =============================================================================

def decode_batch_with_kis(batch_seq, kis, filename):
    """Run the rest of the pipeline (PrimerSort → Consensus → RS → NullFind →
    XOR cascade) on a contiguous prefix of reads, using the KIS values that
    were already computed on the FULL data.

    The PrimerSort step is repeated on the batch (cheap, scales linearly).
    KeyFind/Keysort are SKIPPED — that's the whole point of Option 1.

    Returns (success: bool, timings: dict, results: dict).
    """
    F_primers = kis['F_primers']
    R_primers = kis['R_primers']
    OligoLen = kis['OligoLen']
    idxlen = kis['idxlen']
    ncolpo = kis['ncolpo']; nbitpo = kis['nbitpo']
    ncolco = kis['ncolco']; nbitco = kis['nbitco']
    k = kis['k']; n = kis['n']; m = kis['m']

    timings = {}
    results = {'OligoLen': OligoLen, 'idxlen': idxlen}

    try:
        # PrimerSort on the BATCH only.
        t = time.time()
        ps = PrimerSort(batch_seq, F_primers, R_primers)
        timings['PrimerSort_batch'] = round(time.time() - t, 3)
        results['valid_reads'] = len(ps)
        if len(ps) == 0:
            results['error'] = 'no valid reads after PrimerSort'
            return False, timings, results

        # Consensus using KNOWN OligoLen/idxlen.
        t = time.time()
        MIndex = Consensus(ps, idxlen, OligoLen)
        Depth = np.array(
            [MIndex[i].shape[0] if MIndex[i] is not None else 0
             for i in range(len(MIndex))])
        timings['Consensus'] = round(time.time() - t, 3)
        results['mean_depth'] = round(float(np.mean(Depth)), 1)
        results['min_depth'] = int(np.min(Depth))
        results['max_depth'] = int(np.max(Depth))
        results['n_low_depth_le16'] = int(np.sum(Depth <= 16))

        # RS decoding
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

        # NullFind needs the original (non-PrimerSort'd) batch — same as Final.
        t = time.time()
        term = ncolpo * 3 * nbitpo + ncolco * 3 * nbitco
        Nullkey_1 = term + 1 if term % 2 == 1 else term
        s_val = payload.shape[1]
        Nullkey_2 = Nullkey_1 // 2

        Null_count = 0; nk1 = Nullkey_1
        while nk1 % s_val != 0:
            nk1 += 1
            Null_count += 1
        Null_xor_count = 0; nk2 = Nullkey_2
        while nk2 % s_val != 0:
            nk2 += 1
            Null_xor_count += 1

        data_Nd, xor_Nd, dNc, xNc = NullFind(
            batch_seq, Null_count, Null_xor_count, s_val, F_primers)
        bindataN = _dna_to_binary(data_Nd[20:])
        binxorN = _dna_to_binary(xor_Nd[20:])
        timings['NullFind'] = round(time.time() - t, 3)

        # update Depth from null/xor counts (matches Final main)
        tri_row_idx = payload[:-1, :].size * 2 // (3 * s_val) - 1
        if tri_row_idx < len(Depth) and dNc:
            Depth[tri_row_idx] = sum(len(d) for d in dNc if d)
        if len(Depth) > 1 and xNc:
            Depth[-2] = sum(len(x) for x in xNc if x)

        # 3-stage XOR Recovery cascade (Global → Depth-cut → Local)
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

        # Global RS
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

        # Depth-cut + XOR
        Depthcut = 16
        ErrD = np.where(Depth <= Depthcut)[0]
        try:
            err_mask = RSerrinf[ErrD] != 0 if len(ErrD) > 0 else np.array([], dtype=bool)
            err_indices = ErrD[err_mask] if len(err_mask) > 0 else np.array([], dtype=int)
            if len(err_indices) > 0:
                Temp_mat[err_indices, :] = 0
            xtp, _, _, _, XA, XB, XC = XORrun(
                Temp_mat, s_val, Null_count, Null_xor_count,
                ncolpo, ncolco, nbitpo, nbitco)
            ECC_flat = xtp.flatten(order='C')
            if Null_count > 0:
                ECC_flat = ECC_flat[:-Null_count]
            tc, tp = tridecode(ECC_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
            trid = trierrXOR(Temp_mat, tc, tp, Null_count, Null_xor_count,
                             ncolpo, nbitpo, ncolco, nbitco)
            recovery_stage = 'Depth'
        except Exception as e:
            results.setdefault('stage_errors', {})['Depth'] = str(e)

        # Local RS (nestedXOR + reRSdec)
        try:
            if len(err_indices) > 0:
                Temp_mat[err_indices, :] = 0
            xtp, _, tri_p, xor_p, XA, XB, XC = XORrun(
                Temp_mat, s_val, Null_count, Null_xor_count,
                ncolpo, ncolco, nbitpo, nbitco)
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
                Temp_mat, s_val, Null_count, Null_xor_count,
                ncolpo, ncolco, nbitpo, nbitco)
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

        timings['total_decode'] = round(
            sum(v for v in timings.values() if isinstance(v, (int, float))), 3)
        return success, timings, results

    except Exception as e:
        results['error'] = f'{type(e).__name__}: {e}'
        results['success'] = False
        timings['total_decode'] = round(
            sum(v for v in timings.values() if isinstance(v, (int, float))), 3)
        return False, timings, results


# =============================================================================
# Driver
# =============================================================================

def run_table(seq_oriented, kis, batch_sizes, filename, output_dir):
    """Decode each batch size; write a per-batch summary table."""
    rows = []
    for bs in batch_sizes:
        bs_clamp = min(bs, len(seq_oriented))
        batch = seq_oriented[:bs_clamp]
        print(f'\n[batch] size={bs_clamp:,} (front of FASTQ)')
        t0 = time.time()
        ok, timings, results = decode_batch_with_kis(batch, kis, filename)
        wall = round(time.time() - t0, 2)
        row = {
            'batch_size': bs_clamp,
            'wall_seconds': wall,
            'success': ok,
            'recovery_stage': results.get('recovery_stage'),
            'mean_depth': results.get('mean_depth'),
            'rs_uncorrectable': results.get('rs_uncorrectable'),
            'n_low_depth_le16': results.get('n_low_depth_le16'),
            'vertices': results.get('vertices'),
            'faces': results.get('faces'),
            'valid_reads': results.get('valid_reads'),
            't_PrimerSort_batch': timings.get('PrimerSort_batch'),
            't_Consensus': timings.get('Consensus'),
            't_RS': timings.get('RS'),
            't_NullFind': timings.get('NullFind'),
            't_XOR': timings.get('XOR'),
            't_total_decode': timings.get('total_decode'),
        }
        rows.append(row)
        print(f'  → ok={ok} stage={results.get("recovery_stage")}  '
              f'V/F={results.get("vertices", "?")}/{results.get("faces", "?")}  '
              f't_decode={timings.get("total_decode", "?")}s  wall={wall}s')

    return rows


def run_grow(seq_oriented, kis, grow_from, grow_step, grow_max, filename, output_dir):
    """Grow the batch from `grow_from` by `grow_step` until decoding succeeds."""
    rows = []
    bs = grow_from
    success = False
    while bs <= grow_max and bs <= len(seq_oriented):
        batch = seq_oriented[:bs]
        print(f'\n[grow] try batch_size={bs:,}')
        t0 = time.time()
        ok, timings, results = decode_batch_with_kis(batch, kis, filename)
        wall = round(time.time() - t0, 2)
        rows.append({
            'batch_size': bs, 'wall_seconds': wall, 'success': ok,
            'recovery_stage': results.get('recovery_stage'),
            'mean_depth': results.get('mean_depth'),
            'rs_uncorrectable': results.get('rs_uncorrectable'),
            'n_low_depth_le16': results.get('n_low_depth_le16'),
            'vertices': results.get('vertices'),
            'faces': results.get('faces'),
            'valid_reads': results.get('valid_reads'),
            't_total_decode': timings.get('total_decode'),
        })
        print(f'  → ok={ok} stage={results.get("recovery_stage")}  '
              f'V/F={results.get("vertices", "?")}/{results.get("faces", "?")}  '
              f't_decode={timings.get("total_decode", "?")}s  wall={wall}s')
        if ok:
            success = True
            break
        bs += grow_step
    return rows, success


def write_results(rows, kis, mode, args, output_dir):
    import pandas as pd
    os.makedirs(output_dir, exist_ok=True)
    df = pd.DataFrame(rows)
    xlsx = os.path.join(output_dir, f'option1_{FILENAME}_{mode}.xlsx')
    json_path = os.path.join(output_dir, f'option1_{FILENAME}_{mode}.json')
    with pd.ExcelWriter(xlsx, engine='openpyxl') as w:
        df.to_excel(w, sheet_name='batches', index=False)
        kis_view = {k: v for k, v in kis.items() if k != 'seq_oriented'
                    and not k.startswith('F_') and not k.startswith('R_')}
        pd.DataFrame(list(kis_view.items()), columns=['key', 'value']).to_excel(
            w, sheet_name='kis_constants', index=False)
    with open(json_path, 'w', encoding='utf-8') as f:
        kis_safe = {k: v for k, v in kis.items() if k != 'seq_oriented'
                    and not k.startswith('F_') and not k.startswith('R_')}
        json.dump({
            'mode': mode,
            'kis': kis_safe,
            'rows': rows,
            'meta': vars(args),
        }, f, indent=2, default=str)
    print(f'\n[wrote] {xlsx}')
    print(f'[wrote] {json_path}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['table', 'grow'], default='table',
                    help='table = sweep specific batch sizes; '
                         'grow = increase from --grow-from until success.')
    ap.add_argument('--batches', type=str, default=None,
                    help=f'Comma-separated batch sizes (default {DEFAULT_BATCH_SIZES}).')
    ap.add_argument('--grow-from', type=int, default=DEFAULT_GROW_FROM)
    ap.add_argument('--grow-step', type=int, default=DEFAULT_GROW_STEP)
    ap.add_argument('--grow-max', type=int, default=DEFAULT_GROW_MAX)
    ap.add_argument('--output-dir', type=str, default=OUTPUT_DIR)
    args = ap.parse_args()

    print('=' * 78)
    print(f'Figure 4 Option 1 — Front-batch decoding with PRE-COMPUTED KIS')
    print(f'  dataset    : {FILENAME}')
    print(f'  mode       : {args.mode}')
    print(f'  output dir : {args.output_dir}')
    print('=' * 78)

    # 1. Load FULL FASTQ (sequential, NOT random).
    t0 = time.time()
    seq = read_fastq_file(R1_PATH) + read_fastq_file(R2_PATH)
    print(f'[load] FASTQ R1+R2 = {len(seq):,} reads ({round(time.time()-t0, 1)}s)')

    # 2. Pre-compute KIS on the WHOLE FASTQ (so partial-batch failures are
    #    only attributable to depth, not Keysort variance).
    kis = precompute_kis(seq, FILENAME)
    seq_oriented = kis['seq_oriented']  # use the eversion-corrected reads

    # 3. Run the chosen mode.
    if args.mode == 'table':
        if args.batches:
            batches = [int(x) for x in args.batches.split(',') if x.strip()]
        else:
            batches = DEFAULT_BATCH_SIZES
        rows = run_table(seq_oriented, kis, batches, FILENAME, args.output_dir)
    else:
        rows, _ = run_grow(seq_oriented, kis, args.grow_from, args.grow_step,
                           args.grow_max, FILENAME, args.output_dir)

    write_results(rows, kis, args.mode, args, args.output_dir)

    # 4. Console summary table (Figure-2 style).
    print('\n' + '=' * 88)
    print(f"{'batch':>10} {'reads_OK':>10} {'depth':>8} {'rs_unc':>7} "
          f"{'stage':>7} {'V/F':>14} {'decode_s':>10} {'wall_s':>8}")
    print('-' * 88)
    for r in rows:
        bs = r['batch_size']
        vr = r.get('valid_reads', '?')
        d = r.get('mean_depth', '?')
        ru = r.get('rs_uncorrectable', '?')
        st = r.get('recovery_stage') or 'FAIL'
        vf = f"{r.get('vertices', '?')}/{r.get('faces', '?')}"
        td = r.get('t_total_decode', '?')
        w = r.get('wall_seconds', '?')
        print(f"{bs:>10,} {vr:>10} {d:>8} {ru:>7} {st:>7} {vf:>14} {td:>10} {w:>8}")
    print('=' * 88)
    print(f"\nKIS constants used (pre-computed on FULL FASTQ):")
    print(f"  OligoLen = {kis['OligoLen']}, idxlen = {kis['idxlen']}")
    print(f"  ncolpo/nbitpo = {kis['ncolpo']}/{kis['nbitpo']}, "
          f"ncolco/nbitco = {kis['ncolco']}/{kis['nbitco']}")
    print(f"  KIS pre-compute time: {kis['time_kis_total']}s "
          f"(amortized over all batch attempts)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
