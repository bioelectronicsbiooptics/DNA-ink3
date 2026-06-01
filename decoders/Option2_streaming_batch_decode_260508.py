#!/usr/bin/env python3
# =============================================================================
# Figure 4 — Option 2: FASTQ line-streamed batched decoding
# Created: 2026-05-08
# =============================================================================
# IDEA
#   Genuinely streaming. The FASTQ is read 4 lines at a time (one record per
#   block), and reads are accumulated into a growing buffer. Every BATCH_READS
#   new reads (default 2000 reads = 8000 lines), we PAUSE and try to decode
#   the accumulated buffer. As soon as decoding succeeds we STOP — we do
#   not load the rest of the FASTQ.
#
#   Contrast with Option 1: Option 1 loads the full FASTQ first and pre-
#   computes KIS to isolate the depth question. Option 2 makes no such
#   assumption — each retry attempts the FULL pipeline (KeyFind included),
#   the way a pure streaming consumer would have to.
#
# FASTQ FORMAT NOTES
#   Each record is 4 lines: '@header', 'sequence', '+', 'quality'.
#   2000 reads ≈ 8000 lines IFF the file has no stray blank lines.
#   This implementation reads in record-units (always 4 lines) — robust
#   against blank lines, and the +/quality lines are correctly skipped.
#
# DATASET / SIZE
#   Same as before: J6 (Jikji6_1.fastq + Jikji6_2.fastq).
#   Default batch = 2000 reads per accumulation step (user-specified).
#
# Pipeline reuse:
#   Imports correctness primitives from Figure1/python/dna_decode_baseline_
#   260507_Final.py — does NOT modify Final. The orchestration logic from
#   Final.main() is copied locally into `decode_full()` so we can call it
#   on a partial buffer without changing Final.
#
# Usage:
#   python3 Figure4/python/Option2_streaming_batch_decode_260508.py
#   python3 Figure4/python/Option2_streaming_batch_decode_260508.py --batch-reads 2000
#   python3 Figure4/python/Option2_streaming_batch_decode_260508.py --max-batches 50
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

from dna_decode_baseline_260507_Final import (  # noqa: E402
    PrimerSelect, ReverseEvert, PrimerSort, KeyFind, Keysort,
    Consensus, Binarization, RSdecoding,
    NullFind, XORrun, nestedXOR, reRSdec, tridecode, trierrXOR,
    _make_triangulation, _dna_to_binary, _de2bi,
    HAS_TRIMESH,
)


# =============================================================================
# CONFIG
# =============================================================================

DATA_ROOT = '/Users/taeseokkang/Desktop/DNA ink3/3D culture heritage seq data'


def resolve_dataset(filename):
    """Map dataset short-name -> (R1_path, R2_path). Mirrors baseline logic."""
    if filename in ('Jikji1', 'Jikji2', 'Jikji3', 'Jikji4', 'Jikji5', 'Jikji6', 'Jikji7'):
        folder = 'Jikji'
        r1 = os.path.join(DATA_ROOT, folder, f'{filename}_1.fastq')
        r2 = os.path.join(DATA_ROOT, folder, f'{filename}_2.fastq')
    elif filename in ('J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7'):
        folder = 'Jikji'
        n = filename[1:]
        r1 = os.path.join(DATA_ROOT, folder, f'Jikji{n}_1.fastq')
        r2 = os.path.join(DATA_ROOT, folder, f'Jikji{n}_2.fastq')
    elif filename == 'tripitaka':
        folder = 'Tripitaka_extracted'
        r1 = os.path.join(DATA_ROOT, folder, 'tripitaka_1.fastq')
        r2 = os.path.join(DATA_ROOT, folder, 'tripitaka_2.fastq')
    elif filename == 'liberty':
        folder = 'Statue of Liberty'
        r1 = os.path.join(DATA_ROOT, folder, 'liberty_1.fastq')
        r2 = os.path.join(DATA_ROOT, folder, 'liberty_2.fastq')
    elif filename == 'pisa':
        folder = 'Tower of Pisa'
        r1 = os.path.join(DATA_ROOT, folder, 'pisa_1.fastq')
        r2 = os.path.join(DATA_ROOT, folder, 'pisa_2.fastq')
    elif filename == 'seokga':
        folder = 'seokgatap'
        r1 = os.path.join(DATA_ROOT, folder, 'seokga_1.fastq')
        r2 = os.path.join(DATA_ROOT, folder, 'seokga_2.fastq')
    elif filename == 'venus':
        folder = 'Venus Statue'
        r1 = os.path.join(DATA_ROOT, folder, 'V_1.fastq')
        r2 = os.path.join(DATA_ROOT, folder, 'V_2.fastq')
    elif filename == 'goldengate':
        folder = 'Golden gate'
        r1 = os.path.join(DATA_ROOT, folder, 'goldengate_1.fastq')
        r2 = os.path.join(DATA_ROOT, folder, 'goldengate_2.fastq')
    else:
        raise ValueError(f"Unknown dataset: {filename!r}")
    return r1, r2


FILENAME = 'Jikji6'
R1_PATH, R2_PATH = resolve_dataset(FILENAME)

OUTPUT_DIR = os.path.join(HERE, 'Outputs', 'option2_streaming_260508')

DEFAULT_BATCH_READS = 2000   # user-specified accumulation step
DEFAULT_MAX_BATCHES = 5000   # safety cap (10M reads max accumulator)


# =============================================================================
# Robust line-batched FASTQ stream
# =============================================================================

def stream_fastq_records(filepath):
    """Yield sequence strings from a FASTQ file in record units (4 lines).

    Robust to blank lines that may appear between records: if the next
    non-empty line does not start with '@' we skip it. Format guard:
    line 1 must start with '@', line 3 must start with '+'.
    """
    import gzip
    opener = gzip.open if filepath.endswith('.gz') else open
    with opener(filepath, 'rt') as f:
        while True:
            # Read 4 lines, tolerating leading blank lines.
            lines = []
            while len(lines) < 4:
                line = f.readline()
                if not line:
                    break
                if line.strip() == '' and not lines:
                    continue   # skip pre-record blank line
                lines.append(line.rstrip('\n'))
            if len(lines) < 4:
                break
            if not lines[0].startswith('@') or not lines[2].startswith('+'):
                # malformed record — skip and try to resync
                continue
            yield lines[1]  # the sequence line


def stream_fastq_pair(r1_path, r2_path, batch_reads):
    """Yield successive batches of ~batch_reads sequences from R1 then R2.

    Each yielded batch is a list[str]. The total sequence is reads from R1
    concatenated with reads from R2 (matches Final's `read_fastq_file(r1) +
    read_fastq_file(r2)` ordering).
    """
    buf = []
    for path in (r1_path, r2_path):
        for s in stream_fastq_records(path):
            buf.append(s)
            if len(buf) >= batch_reads:
                yield buf
                buf = []
    if buf:
        yield buf  # tail batch (may be smaller than batch_reads)


# =============================================================================
# Decode the accumulated buffer using the FULL Final pipeline
# =============================================================================

def decode_full(seq, filename):
    """Run Final's full main()-equivalent on the given sequence list.

    Returns (success, timings, results). Mirrors the steps of
    dna_decode_baseline_260507_Final.main() exactly — copied here so we can
    invoke it with a partial buffer without modifying Final.
    """
    timings = {}
    results = {}

    try:
        F_primers, R_primers = PrimerSelect(filename)

        # 1: Evert
        t = time.time()
        forward_section = [s[:21] if len(s) >= 21 else s for s in seq]
        tevert_seq, reversed_log = ReverseEvert(seq, forward_section, R_primers, F_primers)
        rev_indices = np.where(reversed_log)[0]
        for i, idx in enumerate(rev_indices):
            if tevert_seq[i] is not None:
                seq[idx] = tevert_seq[i]
        timings['Evert'] = round(time.time() - t, 3)

        # 2: PrimerSort
        t = time.time()
        ps = PrimerSort(seq, F_primers, R_primers)
        timings['PrimerSort'] = round(time.time() - t, 3)
        results['valid_reads'] = len(ps)
        if len(ps) == 0:
            results['error'] = 'no valid reads after PrimerSort'
            return False, timings, results

        # 3: KeyFind / Keysort  ←  this is the streaming-failure mode
        t = time.time()
        k = 156; n = k + 8; m = len(_de2bi(n))
        KeyMap, KmatDNA = KeyFind(ps)
        OligoLen, idxlen, ncolpo, nbitpo, ncolco, nbitco = Keysort(KeyMap, KmatDNA, k, filename)
        timings['KeyFind'] = round(time.time() - t, 3)
        results['OligoLen'] = int(OligoLen)
        results['idxlen'] = int(idxlen)

        # 4: Consensus
        t = time.time()
        MIndex = Consensus(ps, idxlen, OligoLen)
        Depth = np.array([MIndex[i].shape[0] if MIndex[i] is not None else 0
                          for i in range(len(MIndex))])
        timings['Consensus'] = round(time.time() - t, 3)
        results['mean_depth'] = round(float(np.mean(Depth)), 1)
        results['min_depth'] = int(np.min(Depth))
        results['max_depth'] = int(np.max(Depth))
        results['n_low_depth_le16'] = int(np.sum(Depth <= 16))

        # 5: RS
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

        # 6: NullFind
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

        data_Nd, xor_Nd, dNc, xNc = NullFind(seq, Null_count, Null_xor_count, s_val, F_primers)
        bindataN = _dna_to_binary(data_Nd[20:])
        binxorN = _dna_to_binary(xor_Nd[20:])
        timings['NullFind'] = round(time.time() - t, 3)

        tri_row_idx = payload[:-1, :].size * 2 // (3 * s_val) - 1
        if tri_row_idx < len(Depth) and dNc:
            Depth[tri_row_idx] = sum(len(d) for d in dNc if d)
        if len(Depth) > 1 and xNc:
            Depth[-2] = sum(len(x) for x in xNc if x)

        # 7: 3-stage cascade
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
# Driver — stream + accumulate + decode
# =============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dataset', nargs='?', default='Jikji6',
                    help='Dataset name (default: Jikji6). '
                         'Supported: Jikji1..7, J1..J7, tripitaka, liberty, '
                         'pisa, seokga, venus, goldengate.')
    ap.add_argument('--batch-reads', type=int, default=DEFAULT_BATCH_READS,
                    help='Reads per streaming batch (default 2000 ≈ 8000 lines).')
    ap.add_argument('--max-batches', type=int, default=DEFAULT_MAX_BATCHES,
                    help='Hard cap on batch count (safety; default 5000).')
    ap.add_argument('--decode-every', type=int, default=1,
                    help='Try decode every N batches. 1 = after every batch '
                         '(safest, slow). E.g. 50 = decode every 100k reads.')
    ap.add_argument('--stop-on-success', action='store_true', default=True)
    ap.add_argument('--no-stop-on-success', dest='stop_on_success',
                    action='store_false',
                    help='Continue accumulating beyond first success (collects '
                         'time scaling data).')
    ap.add_argument('--gold-vf', type=str, default='',
                    help='Comma-separated "V,F" the mesh must match before we '
                         'count a decode as truly successful. At very low read '
                         'counts the cascade can produce non-exception meshes '
                         'with wrong vertex counts. Empty string (default) '
                         'disables the check; the cascade-success boolean alone '
                         'is then trusted.')
    ap.add_argument('--output-dir', type=str, default=OUTPUT_DIR)
    args = ap.parse_args()

    # Resolve dataset paths
    global FILENAME, R1_PATH, R2_PATH
    FILENAME = args.dataset
    R1_PATH, R2_PATH = resolve_dataset(FILENAME)

    print('=' * 78)
    print('Figure 4 Option 2 — FASTQ line-streamed batched decoding')
    print(f'  dataset      : {FILENAME}')
    print(f'  R1           : {R1_PATH}')
    print(f'  R2           : {R2_PATH}')
    print(f'  batch_reads  : {args.batch_reads:,}  (lines per batch ≈ '
          f'{args.batch_reads * 4:,} for 4-line FASTQ records)')
    print(f'  max_batches  : {args.max_batches:,}')
    print(f'  decode_every : {args.decode_every}  '
          f'(decode attempt every {args.batch_reads * args.decode_every:,} reads)')
    print(f'  stop_on_OK   : {args.stop_on_success}')
    print(f'  output_dir   : {args.output_dir}')
    print('=' * 78)

    # Gold V/F gate — guards against "false success" where the recovery
    # cascade builds a non-exception triangulation from grossly insufficient
    # data (we observed V/F = 433/864 at 2k reads; gold is 434/864 on J6).
    gold_vf = None
    if args.gold_vf:
        try:
            gv, gf = (int(x) for x in args.gold_vf.split(','))
            gold_vf = (gv, gf)
        except Exception:
            print(f'[warn] could not parse --gold-vf "{args.gold_vf}"; disabling')

    accumulator = []
    rows = []
    batches_seen = 0
    success_at = None

    grand_t0 = time.time()
    for batch in stream_fastq_pair(R1_PATH, R2_PATH, args.batch_reads):
        batches_seen += 1
        t_stream = time.time()
        accumulator.extend(batch)
        n_total = len(accumulator)
        t_stream_done = round(time.time() - t_stream, 3)

        if batches_seen % args.decode_every != 0 and batches_seen != 1:
            continue

        # Try a decode on the current accumulator. Pass a COPY so the pipeline
        # can mutate orientation freely without affecting our growing list.
        seq_copy = list(accumulator)
        t_decode = time.time()
        ok, timings, results = decode_full(seq_copy, FILENAME)
        decode_wall = round(time.time() - t_decode, 2)

        v_got = results.get('vertices')
        f_got = results.get('faces')
        mesh_matches_gold = (
            ok and gold_vf is not None
            and v_got == gold_vf[0] and f_got == gold_vf[1]
        ) if gold_vf is not None else None
        # If gold_vf is set, only count "true" successes for stop_on_success.
        true_success = (mesh_matches_gold is True) if gold_vf else ok

        row = {
            'batch_index': batches_seen,
            'total_reads_streamed': n_total,
            'stream_seconds_this_step': t_stream_done,
            'decode_wall_seconds': decode_wall,
            'cumulative_wall_seconds': round(time.time() - grand_t0, 2),
            'success': ok,
            'mesh_matches_gold': mesh_matches_gold,
            'recovery_stage': results.get('recovery_stage'),
            'OligoLen_detected': results.get('OligoLen'),
            'mean_depth': results.get('mean_depth'),
            'rs_uncorrectable': results.get('rs_uncorrectable'),
            'n_low_depth_le16': results.get('n_low_depth_le16'),
            'vertices': results.get('vertices'),
            'faces': results.get('faces'),
            'valid_reads': results.get('valid_reads'),
            't_Evert': timings.get('Evert'),
            't_PrimerSort': timings.get('PrimerSort'),
            't_KeyFind': timings.get('KeyFind'),
            't_Consensus': timings.get('Consensus'),
            't_RS': timings.get('RS'),
            't_NullFind': timings.get('NullFind'),
            't_XOR': timings.get('XOR'),
            't_total_decode': timings.get('total_decode'),
        }
        rows.append(row)

        ol = results.get('OligoLen', '?')
        v = results.get('vertices', '?')
        f_ = results.get('faces', '?')
        st = results.get('recovery_stage') or 'FAIL'
        gate_tag = ''
        if gold_vf is not None:
            gate_tag = ' GOLD✓' if mesh_matches_gold else (' BAD-VF' if ok else '')
        print(f'  [batch {batches_seen:>4}] streamed {n_total:>9,} reads  '
              f'OL={ol}  decode={decode_wall:.2f}s  ok={ok} stage={st:>5} V/F={v}/{f_}{gate_tag}')

        if true_success and success_at is None:
            success_at = (batches_seen, n_total)
            if args.stop_on_success:
                gate = ' (gold V/F match)' if gold_vf else ''
                print(f'\n[stop] success at batch {batches_seen} ({n_total:,} reads){gate}')
                break

        if batches_seen >= args.max_batches:
            print(f'\n[stop] max-batches cap hit ({args.max_batches})')
            break

    # Save
    os.makedirs(args.output_dir, exist_ok=True)
    out_xlsx = os.path.join(args.output_dir, f'option2_{FILENAME}.xlsx')
    out_json = os.path.join(args.output_dir, f'option2_{FILENAME}.json')

    import pandas as pd
    df = pd.DataFrame(rows)
    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as w:
        df.to_excel(w, sheet_name='streaming', index=False)
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump({
            'filename': FILENAME,
            'batch_reads': args.batch_reads,
            'success_at': success_at,
            'rows': rows,
        }, f, indent=2, default=str)
    print(f'\n[wrote] {out_xlsx}')
    print(f'[wrote] {out_json}')

    # Console summary table
    print('\n' + '=' * 100)
    print(f"{'batch':>6} {'reads':>10} {'OL':>5} {'depth':>8} {'rs_unc':>7} "
          f"{'stage':>5} {'V/F':>14} {'decode_s':>10} {'wall_s':>8}")
    print('-' * 100)
    for r in rows:
        bi = r['batch_index']
        rd = r['total_reads_streamed']
        ol = r.get('OligoLen_detected', '?')
        d = r.get('mean_depth', '?')
        ru = r.get('rs_uncorrectable', '?')
        st = r.get('recovery_stage') or 'FAIL'
        vf = f"{r.get('vertices', '?')}/{r.get('faces', '?')}"
        td = r.get('t_total_decode', '?')
        w = r.get('decode_wall_seconds', '?')
        print(f"{bi:>6} {rd:>10,} {str(ol):>5} {d:>8} {ru:>7} {st:>5} {vf:>14} {td:>10} {w:>8}")
    print('=' * 100)
    if success_at:
        bi, rd = success_at
        print(f'\nFirst success: batch #{bi} at {rd:,} reads streamed.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
