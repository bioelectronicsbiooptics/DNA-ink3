# `decoders/` — iSCAN decoder family

Six Python decoder variants used to produce the §1 / §3 / §4 / §5 results of the manuscript. All variants share the same iSCAN 4-stage skeleton (primer sorting → KIS detection → index clustering → consensus + RS/XOR); the differences are which **decoding stage** or **experimental knob** each variant exposes.

---

## File map

| File | § in manuscript | Role |
|---|---|---|
| **`dna_decode_baseline_260507_Final.py`** | §1 | Pure iSCAN 4-stage decoder (positional decoding only — no optimization layers, no down-sampling, no patching). The validator used to generate the v8_11 results across all 13 datasets. **Reference point** for every other variant. |
| **`run_downsampling_experiment.py`** | §3 (3a, horizontal axis) | Drives `dna_decode_baseline_260507_Final.py` over a *range* of coverage targets (e.g. ×1421 → ×500 → ×100 → ×20). Random downsampling via reservoir sampling on the combined R1+R2 stream; one trial per target. Produces the **per-ratio time / accuracy curve**. |
| **`run_downsampling_n100.py`** | §3 (3a/3b/3c, error bars) | Drives `dna_decode_baseline_260507_Final.py` at a *single* operating point (default ratio = 1/5.33 → ×504 coverage) for **N = 100 independent seeds**. Produces the **statistical confidence intervals** underlying the "×504 / 2.47 s / 100 % reconstruction" claim, plus per-trial RS-fail counts and the depth-array dump that feeds Figure 3b. |
| **`Option1_front_batch_decode_260508.py`** | §4 (Option 1) | Patch-based decoding with a **front-load** schedule: the full FASTQ is parsed once into memory, then partitioned by KIS-derived index range into discrete patches and dispatched to multiprocessing workers. Lower latency, higher peak RAM. |
| **`Option2_streaming_batch_decode_260508.py`** | §4 (Option 2) | Patch-based decoding with a **streaming** schedule: the FASTQ is read chunk-by-chunk via `mmap` + `islice`, and each chunk is routed to its target patch worker. Higher latency, near-constant peak RAM (bounded by chunk size, not file size). |
| **`dna_decode_260328_0125.py`** | §5 | iSCAN baseline + the **four-layer optimization stack** (NumPy whole-array vectorization → uint8 byte-level encoding → lookup-table GF(2⁸) Reed–Solomon → Numba JIT consensus). This file produces the **headline "Jikji 6 109 s → 3.9 s (28×)"** benchmark cited throughout §5 and the conclusion. CPU-only; no GPU dependency. |

---

## Narrative chain — how to read this folder

```
§1   dna_decode_baseline_260507_Final.py
        │   └── "positional decoding works" (360× over de novo alignment)
        │
§3   dna_decode_baseline_260507_Final.py  +  run_downsampling_experiment.py  +  run_downsampling_n100.py
        │   └── "×504 coverage is sufficient" (×5.33 downsampling, 100% recovery)
        │
§4   dna_decode_baseline_260507_Final.py  +  Option1_front_batch_decode_260508.py  /  Option2_streaming_batch_decode_260508.py
        │   └── "patch-level parallelism preserves index purity" (99.87% at ×504)
        │
§5   dna_decode_260328_0125.py
            └── "algorithmic optimization alone reaches 28×" (109 s → 3.9 s, no GPU)
```

Each downstream variant **calls** the baseline decoder (`dna_decode_baseline_260507_Final.py`) for its core stages, so the four manuscript claims are produced from a single decoding kernel exercised in four different harnesses. `dna_decode_260328_0125.py` replaces that kernel with the four-layer optimized version while preserving identical input / output contracts.

---

## Quick reference — running each variant

```bash
# §1 — baseline iSCAN on a single FASTQ pair
python decoders/dna_decode_baseline_260507_Final.py --r1 example_data/J6_R1.example.fastq \
                                   --r2 example_data/J6_R2.example.fastq

# §3a — coverage sweep (multiple ratios, one trial each)
python decoders/run_downsampling_experiment.py --targets 1421,500,100,20

# §3b/c — single ratio, 100 seeded trials (statistics)
python decoders/run_downsampling_n100.py --n 100 --ratio 0.18762

# §4 (Option 1) — front-load patch decoding
python decoders/Option1_front_batch_decode_260508.py --patches 10 --workers 10

# §4 (Option 2) — streaming patch decoding
python decoders/Option2_streaming_batch_decode_260508.py --chunk 50000 --patches 10 --workers 10

# §5 — 4-layer optimized decoder (28× headline)
python decoders/dna_decode_260328_0125.py --dataset Jikji6 --full
```

Concrete dataset paths and the small reproducible slice for CI live in `../example_data/`.

---

## What is *not* in this repository

- **GPU stretch variant** (`dna_decode_260508_OptB_GPU_NullOpt_EvertOpt.py`) — Apple Silicon CPU+GPU heterogeneous prototype, withheld until the V15 G-series end-to-end timings are finalised. The §5 manuscript text describes the algorithmic stack only; GPU is reported as a "potential additional layer" but not relied upon for the headline 28× number.
- **Figure-generation code** (SVG fillers, Excel aggregators) — intentionally excluded: this repository carries only the decoding pipeline, not the manuscript figure factory.
- **MATLAB reference pipeline** — the v240108 reference decoder used to establish numerical parity is hosted with the companion wet-lab paper at [bioelectronicsbiooptics/3D-Culture-Heritage](https://github.com/bioelectronicsbiooptics/3D-Culture-Heritage); this repository carries the Python iSCAN3 stack only.
- **Raw sequencing data** — deposited at NCBI BioProject **PRJNA1091804**, not redistributed here.
