# `decoders/` — iSCAN decoder family

Five Python decoder variants used to produce the §1 / §3 / §4 results of the manuscript. All variants share the same iSCAN 4-stage skeleton (primer sorting → KIS detection → index clustering → consensus + RS/XOR); the differences are which **decoding stage** or **experimental knob** each variant exposes.

---

## File map

| File | § in manuscript | Role |
|---|---|---|
| **`dna_decode_baseline_260507_Final.py`** | §1 | Pure iSCAN 4-stage decoder (positional decoding only — no optimization layers, no down-sampling, no patching). The validator used to generate the v8_11 results across all 13 datasets. **Reference point** for every other variant. |
| **`run_downsampling_experiment.py`** | §3 (3a, horizontal axis) | Drives `dna_decode_baseline_260507_Final.py` over a *range* of coverage targets (e.g. ×1421 → ×500 → ×100 → ×20). Random downsampling via reservoir sampling on the combined R1+R2 stream; one trial per target. Produces the **per-ratio time / accuracy curve**. |
| **`run_downsampling_n100.py`** | §3 (3a/3b/3c, error bars) | Drives `dna_decode_baseline_260507_Final.py` at a *single* operating point (default ratio = 1/5.33 → ×504 coverage) for **N = 100 independent seeds**. Produces the **statistical confidence intervals** underlying the "×504 / 13.3 s / 100 % reconstruction" claim, plus per-trial RS-fail counts and the depth-array dump that feeds Figure 3b. |
| **`run_patch_experiment_v3_baseline.py`** ⭐ | §4 (full FASTQ) | **Patch-based parallel decoder (baseline iSCAN)** — partitions reads by oligo index range into `n_patches` chunks with ±100-index overlap, dispatches via `multiprocessing.Pool(fork)`, merges per-patch payloads. Sweeps `patch_{2,4,8,10}`. Operates on **full FASTQ** (no downsampling). Use this to characterise patch parallelism in isolation from §3. |
| **`run_patch_experiment_v4_baseline_downsampled.py`** ⭐ | §4 (manuscript condition) | **Patch-based parallel decoder + ×5.33 downsampling** — same patch architecture as v3 but applies the manuscript's ×5.33 downsampling (`ratio = 0.18762`) BEFORE preprocessing, matching the §4 condition "evaluated under the previously established the lowest downsampled conditions". Produces the **headline §4 patch-parallelism numbers** (e.g. Tripitaka: 11.9 min sequential → 1.8 min patch_4, **6.71× speedup**). |

---

## Narrative chain — how to read this folder

```
§1   dna_decode_baseline_260507_Final.py
        │   └── "positional decoding works"  (Apple M5 Max: J6 130 s end-to-end)
        │
§3   dna_decode_baseline_260507_Final.py  +  run_downsampling_experiment.py  +  run_downsampling_n100.py
        │   └── "×504 coverage is sufficient" (×5.33 downsampling, 100% recovery)
        │
§4   run_patch_experiment_v3_baseline.py             (full FASTQ baseline + patch)
     run_patch_experiment_v4_baseline_downsampled.py  ⭐ manuscript condition
            └── "patch_4 partition is optimal under multiprocessing.fork overhead"
            └── Tripitaka v4: 11.9 min seq → 1.8 min patch_4 (6.71×)
```

Each §3 / §4 variant imports the §1 baseline (`dna_decode_baseline_260507_Final.py`) for its core decoding primitives (Consensus / Binarization / RSdecoding), so the manuscript claims are produced from a single decoding kernel exercised in different harnesses.

---

## Quick reference — running each variant

```bash
# §1 — baseline iSCAN on a single FASTQ pair
python decoders/dna_decode_baseline_260507_Final.py Jikji6

# §3a — coverage sweep (multiple ratios, one trial each)
python decoders/run_downsampling_experiment.py --targets 1421,500,100,20

# §3b/c — single ratio, 100 seeded trials (statistics)
python decoders/run_downsampling_n100.py --n 100 --ratio 0.18762

# §4 — patch parallelism on full FASTQ (baseline + patch)
python decoders/run_patch_experiment_v3_baseline.py Jikji6 tripitaka --n-patches 2,4,8,10

# §4 ⭐ — patch parallelism under manuscript's ×5.33 downsampling condition
python decoders/run_patch_experiment_v4_baseline_downsampled.py Jikji6 tripitaka \
       --n-patches 2,4,8,10 --ratio 0.18762
```

Concrete dataset paths and the small reproducible slice for CI live in `../example_data/`.

---

## Configuring `DATA_ROOT`

All scripts share a single hard-coded constant `DATA_ROOT` pointing at the local FASTQ directory:

```python
DATA_ROOT = '/Users/.../3D culture heritage seq data'
```

Edit this constant in each script before running on a new system, or set up your local FASTQ folder to match the original path. Per-dataset subfolders (`Jikji/`, `Tripitaka_extracted/`, `Statue of Liberty/`, etc.) are listed in `resolve_paired_fastq()` of each variant.

---

## §4 patch-parallel decoder — measurement summary (Apple M5 Max)

The two patch-parallel variants (v3, v4) were benchmarked across the 13 cultural-heritage payloads on **Apple M5 Max (18 cores: 6P + 12E, 128 GB RAM)** running `multiprocessing.Pool(fork)`:

### v4 (manuscript condition: baseline + ×5.33 downsampling + patch) — Tripitaka highlight

| Method | Decode | End-to-end | Speedup (decode) |
|---|---|---|---|
| sequential | 11.9 min | 16.3 min | 1.00× |
| patch_2 | 17.7 min | 22.2 min | 0.67× (under-partitioned) |
| **patch_4** ⭐ | **1.8 min** | **6.2 min** | **6.71×** |
| patch_8 | 5.1 min | 9.6 min | 2.32× |
| patch_10 | 4.1 min | 8.6 min | 2.89× (fork overhead) |

→ **patch_4 is the optimal partition** for `multiprocessing.fork`-based patching on M5 Max.

### v3 (full FASTQ baseline + patch) — Tripitaka highlight

| Method | Decode | End-to-end | Speedup (decode) |
|---|---|---|---|
| sequential | 60.2 min | 75.3 min | 1.00× |
| **patch_4** ⭐ | **9.2 min** | **24.5 min** | **6.52×** |
| patch_10 | 18.2 min | 33.3 min | 3.31× |

→ patch speedup (6.52×) is preserved across input scales (full vs ×5.33 downsampled), confirming that the patch algorithm operates **independently of input size**.

Detailed per-dataset JSON results are produced under `patch_baseline_v3_out/` and `patch_baseline_v4_ds_out/` (parent of `decoders/` in the working directory).

---

## What is *not* in this repository

- **Figure-generation code** (SVG fillers, Excel aggregators) — intentionally excluded: this repository carries only the decoding pipeline, not the manuscript figure factory.
- **MATLAB reference pipeline** — the v240108 reference decoder used to establish numerical parity is hosted with the companion wet-lab paper at [bioelectronicsbiooptics/3D-Culture-Heritage](https://github.com/bioelectronicsbiooptics/3D-Culture-Heritage); this repository carries the Python iSCAN3 stack only.
- **Raw sequencing data** — deposited at NCBI BioProject **PRJNA1091804**, not redistributed here.
