# DNA-ink3 — iSCAN Positional Decoding & Optimization Stack

Implementation of the paper:
**"High data-density, high decoding-speed, and high decoding-accuracy DNA data ink for digital preservation — iSCAN3 optimization stack"** *(in preparation, target: Briefings in Bioinformatics, 2026)*

A position-based DNA data-storage decoder that replaces *de novo* sequence alignment with a localized index-scan, combined with a four-layer algorithmic optimization stack (NumPy vectorization → uint8 byte-level encoding → lookup-table GF(2⁸) Reed–Solomon → Numba JIT consensus). Achieves **28× end-to-end speed-up over the MATLAB reference** on the Jikji 6 payload (109 s → 3.9 s) while preserving 100 % bit-for-bit reconstruction.

---

## Pipeline overview

```
                ┌──────────────────────────────┐
   FASTQ (R1)  ─┤  uint8 ingestion              │
                │  np.frombuffer → N × 150 mat. │
   FASTQ (R2)  ─┤  reverse-complement (R2)      │
                └────────────┬─────────────────┘
                             │
              ┌──────────────▼──────────────┐
              │ Stage 1 — Primer Sorting    │  ≥ 18 nt forward + reverse match
              └──────────────┬──────────────┘
                             │
              ┌──────────────▼──────────────┐
              │ Stage 2 — KIS Detection     │  poly-A anchor at pos 98 (constant time)
              └──────────────┬──────────────┘
                             │
              ┌──────────────▼──────────────┐
              │ Stage 3 — Index Clustering  │  5–10 nt local scan,
              │  (optional patch parallel)  │  argsort + searchsorted
              └──────────────┬──────────────┘
                             │
              ┌──────────────▼──────────────┐
              │ Stage 4 — Consensus + RS    │  uint8 column-mode (Numba JIT)
              │  + outer XOR / dropout      │  reedsolo GF(2⁸) lookup-table
              └──────────────┬──────────────┘
                             ▼
                       Decoded STL
                  (bit-for-bit identical)
```

**Inputs**
- Paired-end FASTQ (Illumina iSeq 100 / HiSeq 2500, 150 nt paired-end).
- Per-dataset primer set definitions (`P1`–`P7`).

**Outputs**
- Decoded binary payload + STL mesh (vertex / facet tables verified against ground-truth hashes from NCBI BioProject **PRJNA1091804**).
- Per-stage wall-clock timing, per-cluster coverage histogram, dropout / XOR recovery report.

---

## Decoder family (decoders/)

Six decoder variants, one per manuscript section, all sharing the same iSCAN 4-stage skeleton:

| File | § | Role |
|---|---|---|
| `decoders/iscan_baseline.py` | §1 | Pure positional 4-stage decoder (no optimization layers). Reference validator used to generate v8_11 results across 13 datasets. |
| `decoders/downsampling_sweep.py` | §3 (3a) | Drives the baseline over a range of coverage targets (×1421 → ×20). Produces the per-ratio time/accuracy curve. |
| `decoders/downsampling_n100.py` | §3 (3b/3c) | Single ratio × N=100 seeds for statistical confidence intervals behind the "×504 / 2.47 s / 100 % reconstruction" claim. |
| `decoders/patch_front_batch.py` | §4 (Option 1) | Patch-based decoder with front-load schedule (parse full FASTQ → partition → dispatch). Lower latency. |
| `decoders/patch_streaming.py` | §4 (Option 2) | Patch-based decoder with streaming schedule (mmap chunk → route to patch worker). Near-constant peak RAM. |
| `decoders/iscan_optimized.py` | §5 | Baseline + four-layer optimization stack (NumPy / uint8 / LUT-RS / Numba JIT). Produces the **headline 109 s → 3.9 s (28×)** result. |

See `decoders/README.md` for the full narrative chain and per-variant CLI examples.

---

## Headline benchmarks

| Dataset | Payload | Oligos | MATLAB reference | iSCAN3 (4-layer) | Speed-up |
|---|---|---|---|---|---|
| Jikji 6 | 0.2 MB | 1,018 | 109 s | **3.9 s** | **28×** |
| Jikji 7 | — | 1,339 | 250 s | 17.5 s | 14× |
| Tripitaka Koreana | 16.6 MB | 92,626 | ~23 min | ~1–2 min | ≥10× |

Down-sampling (×5.33, mean coverage ×2,686 → ×504): J4 34.5 s → **2.47 s** (14×) at 100 % accuracy, NGS cost reduced ~80 %.

Patch-based parallel decoding: index match-rate **99.87 %** at ×504 mean coverage; Tripitaka 92,626-oligo payload partitioned into 10 patches reduces wall-clock from ~23 min to ~3 min on a 10-core CPU.

---

## Datasets

Thirteen cultural-heritage DNA-data-storage payloads (NCBI BioProject **PRJNA1091804**):

| Object group | Code names | Oligo range |
|---|---|---|
| Jikji metal type | `Jikji1` – `Jikji7` | 895 – 1,339 |
| Seokgatap pagoda | `seokga` | 2,586 |
| Cultural landmarks | `goldengate`, `pisa`, `liberty`, `venus` | 35,464 – 63,952 |
| Tripitaka Koreana | `tripitaka` | 92,626 |
| Overhang chemistry variant | `J6_OH` (SRR33084752) | — |

The full FASTQ collection (≈ 230 GB) is hosted at NCBI; a small reproducible slice for smoke testing lives in `example_data/`.

---

## Project structure

```
DNA-ink3/
├── README.md                          # this file
├── LICENSE                            # MIT
├── pyproject.toml                     # Python dependencies + metadata
├── .gitignore
│
├── decoders/                          # Decoder family — six §-aligned variants
│   ├── README.md                      # Per-file narrative + CLI examples
│   ├── iscan_baseline.py              # §1 — pure positional decoder
│   ├── downsampling_sweep.py          # §3 — coverage sweep
│   ├── downsampling_n100.py           # §3 — N=100 seeded trials
│   ├── patch_front_batch.py           # §4 — front-load patch decoder
│   ├── patch_streaming.py             # §4 — streaming patch decoder
│   └── iscan_optimized.py             # §5 — 4-layer optimization stack (28×)
│
├── matlab_ref/                        # MATLAB R2022a reference (v240108)
│   └── Single_Decoding_Time.m
│
├── Manuscript_Drafts_v2/              # Manuscript artifacts
│   ├── 260521_Methods.md              # Methods section draft
│   └── Manuscript_Full_Compiled_v15.md
│
└── example_data/                      # Small reproducible slice for CI / tutorials
    └── README.md                      # Layout + populate instructions
```

---

## Requirements

- **Python** ≥ 3.10 (3.14 recommended; tested on macOS arm64 / Linux x86_64)
- **MATLAB** ≥ R2022a (only required to reproduce the reference baseline in `matlab_ref/`)
- **GCC / Clang** with C99 support (for Numba native compilation)

Core Python dependencies (full list in `pyproject.toml`):

| Package | Min version | Role |
|---|---|---|
| numpy | 2.3 | Whole-array vectorization (Layer 1) |
| pandas | 2.3 | Per-dataset aggregation tables |
| matplotlib | 3.10 | Figure rendering (Agg backend) |
| biopython | 1.85 | FASTQ I/O |
| reedsolo | 1.7 | GF(2⁸) lookup-table Reed–Solomon (Layer 3) |
| numba | 0.60 | JIT consensus loop (Layer 4) |
| scipy | 1.16 | Statistical helpers |
| trimesh | 4.11 | STL re-export and mesh diff |
| xlsxwriter | 3.2 | Per-stage timing tables |

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/bioelectronicsbiooptics/DNA-ink3.git
cd DNA-ink3

# 2. Create and activate a virtual environment
#    --- option A: uv (recommended, byte-level reproducible) ---
uv sync

#    --- option B: pip / venv ---
python -m venv .venv
source .venv/bin/activate     # Linux/macOS
pip install -e .
```

The first import of `decoders/iscan_optimized.py` triggers a one-time Numba JIT compile (≈ 8 s) and a 64 KB GF(2⁸) lookup-table precompute (instantaneous); both are cached on disk and re-used on every subsequent run.

---

## Usage

```bash
# §1 — baseline iSCAN on a single FASTQ pair
python decoders/iscan_baseline.py --r1 example_data/J6_R1.example.fastq \
                                   --r2 example_data/J6_R2.example.fastq

# §3a — coverage sweep
python decoders/downsampling_sweep.py --targets 1421,500,100,20

# §3b/c — single ratio, 100 seeded trials (statistics)
python decoders/downsampling_n100.py --n 100 --ratio 0.18762

# §4 (Option 1) — front-load patch decoding
python decoders/patch_front_batch.py --patches 10 --workers 10

# §4 (Option 2) — streaming patch decoding
python decoders/patch_streaming.py --chunk 50000 --patches 10 --workers 10

# §5 — 4-layer optimized decoder (28× headline)
python decoders/iscan_optimized.py --dataset Jikji6 --full

# MATLAB reference cross-check
matlab -batch "addpath('matlab_ref'); run_smoke_all"
```

See `decoders/README.md` for the full per-variant narrative.

---

## Verification

- Python and MATLAB pipelines agree on all 13 public payloads to **bit-for-bit STL identity** (last verified 2026-05-07).
- A minimal smoke test (`example_data/`, ≈ 1,000 reads, runs in < 1 s) is exercised on every commit.

---

## Citation

```bibtex
@article{kang2026iscan3,
  title   = {High data-density, high decoding-speed, and high decoding-accuracy DNA data ink for digital preservation},
  author  = {Kang, Taeseok and Lim, Doyeon and Lee, Wonjin and Kim, Jungwoo and Huang, Xiaohua and Kim, Jinchul and Song, Youngjun},
  journal = {Briefings in Bioinformatics},
  year    = {2026},
  note    = {Manuscript in preparation},
  url     = {https://github.com/bioelectronicsbiooptics/DNA-ink3}
}
```

The companion ACS Nano article describing the wet-lab DNA-ink platform is:

> Kang, T.; Lim, D.; et al. *High Data-Density, High Decoding-Speed, and High Decoding-Accuracy DNA Data Ink for Digital Preservation.* ACS Nano **2026**, 20, 3444–3460. [DOI 10.1021/acsnano.5c13665](https://doi.org/10.1021/acsnano.5c13665)

---

## License

This repository is released under the [MIT License](LICENSE).
Raw sequencing data are subject to the NCBI SRA submission terms under BioProject **PRJNA1091804**.

---

## Contact

- Corresponding author (computational): **Youngjun Song** (Incheon National University) — yjunsong@inu.ac.kr
- Corresponding author (wet-lab): **Jinchul Kim** (KRIBB) — jinchul@kribb.re.kr
- Repository maintainer: **Taeseok Kang** — labmanager.bebo.inu@stdbioelec.com

Bug reports and feature requests are welcome via the GitHub issue tracker.
