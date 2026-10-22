#!/usr/bin/env python3
"""Read-batch decoding (Figures 4 and 5).

Reads are taken in FASTQ order, R1 reads first and then R2 reads. For each batch
fraction f, only the first int(total * f) reads are read from disk and decoded with
run_final_pipeline() over the full index range. No random sampling is applied.

A batch counts as recovered when every non-degenerate facet of the original STL is
matched by the decoded mesh (coordinate tolerance 1e-3 after a fixed translation) and
no extra facets are present. The translation is taken from the 100% run, so the
first fraction must be 1.0 unless --shift-from points to an earlier output file.

Example:
  python decoders/run_read_batch.py --name Jikji6 \
      --r1 Jikji6_1.fastq --r2 Jikji6_2.fastq --stl Jikji6_6.stl \
      --fractions 1.0 0.1 0.01 0.004 --out read_batch_Jikji6.jsonl
"""
import argparse, contextlib, io, json, os, sys, time
import numpy as np
from scipy.spatial import cKDTree
import trimesh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dna_decode_baseline_260507_Final as DEC
import run_downsampling_experiment as EXP

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--name", required=True, help="dataset key used for primer selection (e.g. Jikji6, tripitaka)")
ap.add_argument("--r1", required=True)
ap.add_argument("--r2", required=True)
ap.add_argument("--stl", required=True, help="original STL file")
ap.add_argument("--fractions", type=float, nargs="+", default=[1.0, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01])
ap.add_argument("--reference", help="designed sequences, one per line (optional; enables consensus accuracy)")
ap.add_argument("--shift-from", help="earlier output file whose first line holds the 100%% translation")
ap.add_argument("--out", required=True, help="output JSONL (appended)")
args = ap.parse_args()
files = [args.r1, args.r2]

orig = trimesh.load(args.stl, process=False)
REF = np.asarray(orig.triangles).reshape(-1, 9)[orig.area_faces > 1e-12]
print(f"reference {args.stl}: {len(REF)} non-degenerate faces", flush=True)

# capture the mesh built by the decoder
CAPTURED = []
_orig_make = DEC._make_triangulation
def _capture(tc, tp):
    m = _orig_make(tc, tp); CAPTURED.append(m); return m
DEC._make_triangulation = _capture
EXP._make_triangulation = _capture

# consensus accuracy: consensus vs designed sequences at positions 21-130
REFBIN = None
if args.reference:
    REFBIN = np.array([DEC._dna_to_binary(l.strip()[20:130].replace("N", "A")) for l in open(args.reference)])
CONS = []
_bin = EXP.Binarization
def _capture_bin(cell_seq):
    b = _bin(cell_seq)
    if REFBIN is not None and len(cell_seq) == len(REFBIN) and not CONS: CONS.append(b)
    return b
EXP.Binarization = _capture_bin

def count_reads(fn):
    with open(fn, "rb") as f:
        return sum(1 for _ in f) // 4

def read_front(n):
    """First n reads of [R1 reads] + [R2 reads], reading no further than needed."""
    seqs = []
    for fn in files:
        if len(seqs) >= n: break
        with open(fn, "r") as f:
            while len(seqs) < n:
                h = f.readline()
                if not h: break
                seqs.append(f.readline().strip()); f.readline(); f.readline()
    return seqs

def compare(mesh, shift):
    if mesh is None:
        return dict(matched_pct=0.0, extra=None, perfect=False)
    tri = np.asarray(mesh.triangles)
    t = (tri + shift).reshape(-1, 9)[mesh.area_faces > 1e-12]
    if len(t) == 0:
        return dict(matched_pct=0.0, extra=0, perfect=False)
    d_ref, _ = cKDTree(t).query(REF, p=np.inf)
    d_tri, _ = cKDTree(REF).query(t, p=np.inf)
    matched = int((d_ref <= 1e-3).sum()); extra = int((d_tri > 1e-3).sum())
    return dict(matched_pct=round(matched / len(REF) * 100, 4), extra=extra,
                perfect=(matched == len(REF) and extra == 0))

counts = [count_reads(f) for f in files]
total = sum(counts)
print(f"{args.name}: R1 {counts[0]:,} + R2 {counts[1]:,} = {total:,} reads", flush=True)

shift = None
if args.shift_from:
    shift = np.array(json.loads(open(args.shift_from).readline())["shift"])
elif args.fractions[0] != 1.0:
    sys.exit("the first fraction must be 1.0 (or give --shift-from)")

for frac in args.fractions:
    n = int(total * frac)
    t0 = time.time()
    reads = read_front(n)
    read_s = time.time() - t0
    CAPTURED.clear(); CONS.clear()
    t1 = time.time()
    with contextlib.redirect_stdout(io.StringIO()):
        success, timings, results = EXP.run_final_pipeline(reads, args.name)
    decode_s = time.time() - t1
    mesh = CAPTURED[-1] if (success and CAPTURED) else None
    if shift is None:  # 100% run: translation between decoder output and original STL
        if mesh is None:
            sys.exit("100% run did not produce a mesh; cannot align to the original STL")
        shift = orig.bounds[0] - mesh.bounds[0]
    cmp = compare(mesh, shift)
    if CONS:
        cmp["consensus_accuracy_pct"] = round(float((CONS[0] == REFBIN).all(axis=1).mean() * 100), 3)
    rec = dict(dataset=args.name, fraction=frac, n_reads=len(reads),
               r1_reads=min(len(reads), counts[0]), r2_reads=max(0, len(reads) - counts[0]),
               read_s=round(read_s, 2), decode_s=round(decode_s, 2), e2e_s=round(read_s + decode_s, 2),
               stage_s=timings, mesh_produced=bool(success), recovery_stage=results.get("recovery_stage"),
               valid_reads=results.get("valid_reads"), mean_depth=results.get("mean_depth"),
               rs_uncorrectable=results.get("rs_uncorrectable"),
               shift=[round(float(x), 5) for x in shift], **cmp, error=results.get("error"))
    with open(args.out, "a") as f: f.write(json.dumps(rec) + "\n")
    print(f"[{frac:.4f}] reads {len(reads):,} read {read_s:.1f}s decode {decode_s:.1f}s "
          f"perfect={cmp['perfect']} matched={cmp['matched_pct']}% extra={cmp['extra']} "
          f"depth={results.get('mean_depth')}", flush=True)
    del reads
