#!/usr/bin/env python3
"""Full FASTQ versus x5.33 random down-sampling (Figure 3a).

Both conditions read all R1 and R2 reads. For down-sampling, a fixed fraction
(0.18762 = 1/5.33) of the reads is drawn at random without replacement (seed 42)
and decoded with run_final_pipeline(). pipeline_s excludes FASTQ loading.
Recovery is judged against the original STL as in run_read_batch.py.

Example:
  python decoders/run_downsampling_check.py --name Jikji6 \
      --r1 Jikji6_1.fastq --r2 Jikji6_2.fastq --stl Jikji6_6.stl --out downsampling_Jikji6.jsonl
"""
import argparse, contextlib, io, json, os, random, sys, time
import numpy as np
from scipy.spatial import cKDTree
import trimesh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dna_decode_baseline_260507_Final as DEC
import run_downsampling_experiment as EXP

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--name", required=True, help="dataset key used for primer selection (e.g. Jikji6)")
ap.add_argument("--r1", required=True)
ap.add_argument("--r2", required=True)
ap.add_argument("--stl", required=True, help="original STL file")
ap.add_argument("--ratio", type=float, default=0.18762, help="fraction of reads kept (default 1/5.33)")
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--out", required=True, help="output JSONL (appended)")
args = ap.parse_args()

orig = trimesh.load(args.stl, process=False)
REF = np.asarray(orig.triangles).reshape(-1, 9)[orig.area_faces > 1e-12]

CAPTURED = []
_orig_make = DEC._make_triangulation
def _capture(tc, tp):
    m = _orig_make(tc, tp); CAPTURED.append(m); return m
DEC._make_triangulation = _capture
EXP._make_triangulation = _capture

def compare(mesh, shift):
    if mesh is None:
        return dict(matched_pct=0.0, extra=None, perfect=False)
    t = (np.asarray(mesh.triangles) + shift).reshape(-1, 9)[mesh.area_faces > 1e-12]
    if len(t) == 0:
        return dict(matched_pct=0.0, extra=0, perfect=False)
    d_ref, _ = cKDTree(t).query(REF, p=np.inf)
    d_tri, _ = cKDTree(REF).query(t, p=np.inf)
    matched = int((d_ref <= 1e-3).sum()); extra = int((d_tri > 1e-3).sum())
    return dict(matched_pct=round(matched / len(REF) * 100, 4), extra=extra,
                perfect=(matched == len(REF) and extra == 0))

shift = None
for label, ratio in (("full", 1.0), ("downsampled", args.ratio)):
    t0 = time.time()
    seq = DEC.read_fastq_file(args.r1) + DEC.read_fastq_file(args.r2)
    read_s = time.time() - t0
    t1 = time.time()
    if ratio < 1.0:
        seq = random.Random(args.seed).sample(seq, int(len(seq) * ratio))
    sample_s = time.time() - t1
    CAPTURED.clear(); t2 = time.time()
    with contextlib.redirect_stdout(io.StringIO()):
        ok, tm, res = EXP.run_final_pipeline(seq, args.name)
    pipe_s = time.time() - t2
    mesh = CAPTURED[-1] if (ok and CAPTURED) else None
    if shift is None:  # full run: translation between decoder output and original STL
        if mesh is None:
            sys.exit("full run did not produce a mesh; cannot align to the original STL")
        shift = orig.bounds[0] - mesh.bounds[0]
    rec = dict(label=label, ratio=ratio, n_reads=len(seq), read_s=round(read_s, 2), sample_s=round(sample_s, 2),
               pipeline_s=round(pipe_s, 2), stage_s=tm, valid_reads=res.get("valid_reads"),
               mean_depth=res.get("mean_depth"), rs_uncorrectable=res.get("rs_uncorrectable"),
               **compare(mesh, shift))
    with open(args.out, "a") as f: f.write(json.dumps(rec) + "\n")
    print(f"{label}: reads {len(seq):,} pipeline {pipe_s:.1f}s depth={rec['mean_depth']} "
          f"perfect={rec['perfect']}", flush=True)
