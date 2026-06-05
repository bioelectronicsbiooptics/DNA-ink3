#!/usr/bin/env python3
"""
Patch-based parallel decoding experiment (v3 — BASELINE + patch).

For §4 verification (Patch-Based Decoding System).

Imports `dna_decode_baseline_260507_Final` — the iSCAN baseline pipeline
that is a direct Python port of the MATLAB original (no V14 4-layer
optimization, no Numba, no uint8 matrix).

This script measures:
  sequential (baseline)  vs  baseline + patch_{2,4,8,10}

so the manuscript's §4 claim ("23 min → 3 min, 7.6×") is checked against
the same code lineage it was projected from (the baseline iSCAN, not V14).

Difference vs v2:
  - Imports baseline_260507_Final (not 260328_0125)
  - Function signatures: no mat_uint8 / ps_mat_u8 / sorted_u8 / boundaries
  - Consensus returns MIndex only (not tuple)
  - No _seqs_to_matrix, no _init_rs_reedsolo
  - Special datasets (pisa, seokga, liberty) skip MapLocation pruning here
    (we measure the no-MapLocation path so all 13 datasets share one timing
    model; full-pipeline pisa/seokga/liberty results are in the existing
    baseline metrics files already).
"""
import sys
import os
import time
import math
import json
import argparse
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from dna_decode_baseline_260507_Final import (
    read_fastq_file, PrimerSelect, ReverseEvert,
    PrimerSort, KeyFind, Keysort, Consensus, Binarization, RSdecoding,
    _de2bi,
)

DATA_ROOT = '/Users/taeseokkang/Desktop/DNA ink3/3D culture heritage seq data'
OVERLAP = 100  # ±100 index overlap


def resolve_paired_fastq(filename):
    """Resolve dataset name -> (r1_path, r2_path). Mirrors baseline main()."""
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
        r1 = os.path.join(DATA_ROOT, f'{filename}_1.fastq')
        r2 = os.path.join(DATA_ROOT, f'{filename}_2.fastq')
    return r1, r2


def run_sequential(ps_seq, idxlen, OligoLen, k, m, n):
    """Sequential baseline decode."""
    t = time.time()
    MIndex = Consensus(ps_seq, idxlen, OligoLen)
    Depth = np.array([MIndex[i].shape[0] if MIndex[i] is not None else 0 for i in range(len(MIndex))])
    SeqtoBin = Binarization(MIndex)
    rsECC, RSerrinf = RSdecoding(SeqtoBin, k, m, n)
    rsECC[RSerrinf < 0, :] = 0
    rsECC[RSerrinf > 4, :] = 0
    for a in range(min(OligoLen, len(rsECC))):
        if np.any(rsECC[a, :] > 1):
            rsECC[a, :] = 0
    payload = rsECC[:, :k - idxlen]
    return payload, Depth, RSerrinf, time.time() - t


def decode_patch(args):
    """Single patch decode (baseline)."""
    patch_seqs, idxlen, oligo_start, oligo_end, k, m, n = args
    patch_oligo_len = oligo_end - oligo_start
    if len(patch_seqs) == 0 or patch_oligo_len <= 0:
        return oligo_start, oligo_end, None, None, None, 0
    t = time.time()
    try:
        MIndex = Consensus(patch_seqs, idxlen, patch_oligo_len)
        SeqtoBin = Binarization(MIndex)
        rsECC, RSerrinf = RSdecoding(SeqtoBin, k, m, n)
        rsECC[RSerrinf < 0, :] = 0
        rsECC[RSerrinf > 4, :] = 0
        payload = rsECC[:, :k - idxlen]
        Depth = np.array([MIndex[i].shape[0] if MIndex[i] is not None else 0 for i in range(len(MIndex))])
        return oligo_start, oligo_end, payload, Depth, RSerrinf, time.time() - t
    except Exception:
        return oligo_start, oligo_end, None, None, None, time.time() - t


def run_patch_parallel(ps_seq, idxlen, OligoLen, k, m, n, n_patches=4):
    """Patch-parallel baseline decode."""
    patch_size = math.ceil(OligoLen / n_patches)
    n_reads = len(ps_seq)
    reads_per_patch = math.ceil(n_reads / n_patches)
    patches = []
    for p in range(n_patches):
        rs = p * reads_per_patch
        re = min((p + 1) * reads_per_patch, n_reads)
        patch_seqs = ps_seq[rs:re]
        oligo_start = max(0, p * patch_size - OVERLAP)
        oligo_end = min(OligoLen, (p + 1) * patch_size + OVERLAP)
        patches.append((patch_seqs, idxlen, oligo_start, oligo_end, k, m, n))

    t_par = time.time()
    ctx = mp.get_context('fork')
    with ctx.Pool(min(n_patches, mp.cpu_count())) as pool:
        results = pool.map(decode_patch, patches)
    t_par_done = time.time()
    patch_times = [r[5] for r in results]

    total_payload = np.zeros((OligoLen, k - idxlen), dtype=np.uint8)
    total_depth = np.zeros(OligoLen, dtype=np.int64)
    total_rserr = np.full(OligoLen, -1, dtype=np.int64)
    for p, r in enumerate(results):
        os_, oe_, payload, Depth, RSerr, _ = r
        if payload is None: continue
        write_start = p * patch_size
        write_end = min(OligoLen, (p + 1) * patch_size)
        local_offset = write_start - os_
        rows = write_end - write_start
        if rows > 0 and payload.shape[0] >= local_offset + rows:
            total_payload[write_start:write_end] = payload[local_offset:local_offset + rows]
            if Depth is not None and len(Depth) >= local_offset + rows:
                total_depth[write_start:write_end] = Depth[local_offset:local_offset + rows]
            if RSerr is not None and len(RSerr) >= local_offset + rows:
                total_rserr[write_start:write_end] = RSerr[local_offset:local_offset + rows]

    total_time = t_par_done - t_par
    return total_payload, total_depth, total_rserr, total_time, t_par_done - t_par, patch_times


def measure_one(filename, n_patches_list, out_dir):
    """Run sequential + patch sweep on one dataset, save JSON, return summary."""
    print("=" * 72)
    print(f"BASELINE + Patch — dataset: {filename}  patches: {n_patches_list}")
    print("=" * 72)

    r1, r2 = resolve_paired_fastq(filename)
    print(f"R1: {r1}")
    print(f"R2: {r2}")

    t0 = time.time()
    seq = read_fastq_file(r1) + read_fastq_file(r2)
    load_t = time.time() - t0
    print(f"[Step 0] loaded reads: {len(seq):,}  ({load_t:.1f}s)")

    F_primers, R_primers = PrimerSelect(filename)

    print("[Step 1] Evert ...")
    t = time.time()
    forward_section = [s[:21] if len(s) >= 21 else s for s in seq]
    tevert_seq, rev_log = ReverseEvert(seq, forward_section, R_primers, F_primers)
    rev_indices = np.where(rev_log)[0]
    for i, idx in enumerate(rev_indices):
        if tevert_seq[i] is not None:
            seq[idx] = tevert_seq[i]
    evert_t = time.time() - t
    print(f"  done: {evert_t:.1f}s  (reverse: {len(rev_indices):,})")

    print("[Step 2] PrimerSort ...")
    t = time.time()
    ps_seq = PrimerSort(seq, F_primers, R_primers)
    primersort_t = time.time() - t
    print(f"  done: {primersort_t:.1f}s  valid: {len(ps_seq):,}")
    if len(ps_seq) == 0:
        print("[error] no valid reads after primer sort — abort")
        return None

    print("[Step 3] KeyFind + Keysort ...")
    t = time.time()
    k = 156
    n = k + 8
    m = len(_de2bi(n))
    KeyMap, KeymatDNA = KeyFind(ps_seq)
    OligoLen, idxlen, ncolpo, nbitpo, ncolco, nbitco = Keysort(KeyMap, KeymatDNA, k, filename)
    key_t = time.time() - t
    print(f"  done: {key_t:.1f}s  OligoLen={OligoLen}  idxlen={idxlen}")

    preproc_t = time.time() - t0
    print(f"\nPreprocessing total: {preproc_t:.1f}s")
    print(f"Valid reads after PrimerSort: {len(ps_seq):,}")

    all_results = []

    print(f"\n{'─'*60}\nSequential (baseline) ...")
    _, _, _, seq_time = run_sequential(ps_seq, idxlen, OligoLen, k, m, n)
    print(f"  decode time: {seq_time:.1f}s")
    all_results.append({'method': 'sequential', 'patches': 1, 'time': round(seq_time, 3)})

    for n_p in n_patches_list:
        print(f"\n{'─'*60}\nPatch {n_p} (±{OVERLAP} overlap) ...")
        _, _, _, total_t, par_t, patch_ts = run_patch_parallel(ps_seq, idxlen, OligoLen, k, m, n, n_patches=n_p)
        speedup = seq_time / total_t if total_t > 0 else 0
        print(f"  total: {total_t:.1f}s  speedup: {speedup:.2f}×")
        all_results.append({
            'method': f'patch_{n_p}',
            'patches': n_p,
            'time': round(total_t, 3),
            'parallel_time': round(par_t, 3),
            'patch_times': [round(t, 3) for t in patch_ts],
            'speedup': round(speedup, 2),
        })

    out_file = os.path.join(out_dir, f'patch_baseline_{filename}.json')
    payload = {
        'dataset': filename,
        'stack': 'baseline_Final-260507',
        'OligoLen': int(OligoLen),
        'idxlen': int(idxlen),
        'n_reads': int(len(seq)),
        'n_valid': int(len(ps_seq)),
        'preprocessing_seconds': {
            'load': round(load_t, 3),
            'evert': round(evert_t, 3),
            'primer_sort': round(primersort_t, 3),
            'key_find_sort': round(key_t, 3),
            'total': round(preproc_t, 3),
        },
        'results': all_results,
    }
    with open(out_file, 'w') as f:
        json.dump(payload, f, indent=2)

    print(f"\n{'='*72}")
    print(f"{'Method':>15} {'Patches':>8} {'Time':>10} {'Speedup':>9}")
    print(f"{'─'*45}")
    for r in all_results:
        print(f"{r['method']:>15} {r['patches']:>8} {r['time']:>9.3f}s {r.get('speedup', 1.0):>7.2f}×")
    print(f"\nSaved: {out_file}\n")
    return payload


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('datasets', nargs='+', help='Dataset names to measure (e.g. Jikji6 tripitaka)')
    ap.add_argument('--output-dir', default=None)
    ap.add_argument('--n-patches', default='2,4,8,10')
    args = ap.parse_args()

    n_patches_list = [int(x) for x in args.n_patches.split(',')]
    out_dir = args.output_dir or os.path.join(HERE, '..', 'patch_baseline_v3_out')
    os.makedirs(out_dir, exist_ok=True)

    print(f">>> v3 BASELINE + Patch sweep")
    print(f">>> datasets: {args.datasets}")
    print(f">>> n_patches: {n_patches_list}")
    print(f">>> output_dir: {out_dir}\n")

    all_summaries = []
    overall_t0 = time.time()
    for ds in args.datasets:
        t_ds = time.time()
        try:
            summary = measure_one(ds, n_patches_list, out_dir)
            if summary:
                all_summaries.append(summary)
                print(f"[{ds}] elapsed for this dataset: {time.time() - t_ds:.1f}s")
                print(f"[overall] elapsed so far: {(time.time() - overall_t0)/60:.1f} min\n")
        except Exception as e:
            print(f"[error] {ds}: {e}")
            import traceback; traceback.print_exc()

    # save combined sweep
    out_combined = os.path.join(out_dir, 'patch_baseline_sweep_all.json')
    with open(out_combined, 'w') as f:
        json.dump({
            'datasets': args.datasets,
            'n_patches': n_patches_list,
            'total_elapsed_seconds': round(time.time() - overall_t0, 1),
            'summaries': all_summaries,
        }, f, indent=2)
    print(f"\n>>> all done.  total elapsed: {(time.time() - overall_t0)/60:.1f} min")
    print(f">>> combined summary: {out_combined}")


if __name__ == '__main__':
    main()
