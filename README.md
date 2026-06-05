# DNA-ink3 — iSCAN Positional Decoding & Optimization Stack

Implementation of the paper:
**"High data-density, high decoding-speed, and high decoding-accuracy DNA data ink for digital preservation — iSCAN3 optimization stack"** *(in preparation, target: Briefings in Bioinformatics, 2026)*

A position-based DNA data-storage decoder that replaces *de novo* sequence alignment with a localized index-scan, combined with a four-layer algorithmic optimization stack (NumPy vectorization → uint8 byte-level encoding → lookup-table GF(2⁸) Reed–Solomon → Numba JIT consensus). Achieves **28× end-to-end speed-up** over the MATLAB reference on the Jikji 6 payload (109 s → 3.9 s) while preserving 100 % bit-for-bit reconstruction.

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

## Decoder family (`decoders/`)

Six Python variants, one per manuscript section, all sharing the same iSCAN 4-stage skeleton:

| File | § | Role |
|---|---|---|
| `decoders/dna_decode_baseline_260507_Final.py` | §1 | Pure positional 4-stage decoder (no optimization layers). Reference validator used to generate v8_11 results across 13 datasets. |
| `decoders/run_downsampling_experiment.py` | §3 (3a) | Drives the baseline over a range of coverage targets (×1421 → ×20). Produces the per-ratio time/accuracy curve. |
| `decoders/run_downsampling_n100.py` | §3 (3b/3c) | Single ratio × N=100 seeds for statistical confidence intervals behind the "×504 / 2.47 s / 100 % reconstruction" claim. |
| `decoders/run_patch_experiment_v3_baseline.py` ⭐ | §4 (full FASTQ) | Patch-based parallel decoder (baseline iSCAN). Partitions reads by oligo index range, dispatches via `multiprocessing.Pool(fork)` with ±100-index overlap. Sweeps `patch_{2,4,8,10}`. |
| `decoders/run_patch_experiment_v4_baseline_downsampled.py` ⭐ | §4 (manuscript condition) | Patch-based parallel decoder + ×5.33 downsampling. Produces the headline §4 numbers (Tripitaka: 11.9 min seq → 1.8 min patch_4, **6.71× speedup** ≈ manuscript 7.6× claim). |
| `decoders/dna_decode_260328_0125.py` | §5 | Baseline + four-layer optimization stack (NumPy / uint8 / LUT-RS / Numba JIT). Produces the **headline 109 s → 3.9 s (28×)** result on Jikji 6. |

See [`decoders/README.md`](decoders/README.md) for the full narrative chain and per-variant CLI examples.

---

## Headline benchmarks

| Dataset | Payload | Oligos | MATLAB reference | iSCAN3 (4-layer) | Speed-up |
|---|---|---|---|---|---|
| Jikji 6 | 0.2 MB | 1,018 | 109 s | **3.9 s** | **28×** |
| Jikji 7 | — | 1,339 | 250 s | 17.5 s | 14× |
| Tripitaka Koreana | 16.6 MB | 92,626 | ~23 min | ~1–2 min | ≥10× |

- **Down-sampling** (×5.33, mean coverage ×2,686 → ×504): J6 34.5 s → **2.47 s** (≈14×) at 100 % accuracy on the companion ACS Nano workstation; NGS cost reduced ~80 %.
- **Patch-based parallel decoding** (`run_patch_experiment_v4_baseline_downsampled.py`, Apple M5 Max baseline + ×5.33 downsampling): Tripitaka (92,626 oligos) partitioned into **4 patches** (optimal under `multiprocessing.fork` overhead) reduces decode wall-clock from **11.9 min sequential to 1.8 min** — **6.71× speedup** with index match-rate ≥ 99.87 %.

---

## Datasets

Thirteen cultural-heritage DNA-data-storage payloads, all hosted at NCBI BioProject **[PRJNA1091804](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1091804)**:

| Object group | Code names | Oligo range |
|---|---|---|
| Jikji metal type | `Jikji1` – `Jikji7` | 895 – 1,339 |
| Seokgatap pagoda | `seokga` | 2,586 |
| Cultural landmarks | `goldengate`, `pisa`, `liberty`, `venus` | 35,464 – 63,952 |
| Tripitaka Koreana | `tripitaka` | 92,626 |

Raw FASTQ files (≈ 230 GB total) are **not redistributed** through this repository — fetch them directly from NCBI SRA.

---

## Project structure

```
DNA-ink3/
├── README.md                          # this file
├── LICENSE                            # MIT
├── .gitignore
│
└── decoders/                          # Decoder family — six §-aligned variants
    ├── README.md                      # Per-file narrative + CLI examples
    ├── dna_decode_baseline_260507_Final.py              # §1 — pure positional decoder
    ├── run_downsampling_experiment.py                   # §3 — coverage sweep
    ├── run_downsampling_n100.py                         # §3 — N=100 seeded trials
    ├── run_patch_experiment_v3_baseline.py              # §4 — patch parallelism (full FASTQ)
    ├── run_patch_experiment_v4_baseline_downsampled.py  # §4 ⭐ — patch + ×5.33 downsampling (headline)
    └── dna_decode_260328_0125.py                        # §5 — 4-layer optimization stack (28×)
```

---

## Requirements

- **Python** ≥ 3.10 (3.14 recommended; tested on macOS arm64 / Linux x86_64)
- **GCC / Clang** with C99 support (for Numba native compilation)

Python dependencies:

| Package | Min version | Role |
|---|---|---|
| numpy | 2.3 | Whole-array vectorization (Layer 1) |
| pandas | 2.3 | Per-dataset aggregation tables |
| biopython | 1.85 | FASTQ I/O |
| reedsolo | 1.7 | GF(2⁸) lookup-table Reed–Solomon (Layer 3) |
| numba | 0.60 | JIT consensus loop (Layer 4) |
| scipy | 1.16 | Statistical helpers |
| trimesh | 4.11 | STL re-export and mesh diff |

Quick install in a fresh virtual environment:

```bash
git clone https://github.com/bioelectronicsbiooptics/DNA-ink3.git
cd DNA-ink3
python -m venv .venv && source .venv/bin/activate
pip install numpy pandas biopython reedsolo numba scipy trimesh
```

The first run of `decoders/dna_decode_260328_0125.py` triggers a one-time Numba JIT compile (≈ 8 s) and a 64 KB GF(2⁸) lookup-table precompute (instantaneous); both are cached on disk and re-used on every subsequent run.

---

## Usage

```bash
# §1 — baseline iSCAN on a paired FASTQ
python decoders/dna_decode_baseline_260507_Final.py --r1 Jikji6_R1.fastq --r2 Jikji6_R2.fastq

# §3a — coverage sweep
python decoders/run_downsampling_experiment.py --targets 1421,500,100,20

# §3b/c — single ratio, 100 seeded trials (statistics)
python decoders/run_downsampling_n100.py --n 100 --ratio 0.18762

# §4 — patch parallelism on full FASTQ (baseline + patch)
python decoders/run_patch_experiment_v3_baseline.py Jikji6 tripitaka --n-patches 2,4,8,10

# §4 ⭐ — patch parallelism under manuscript's ×5.33 downsampling condition (headline)
python decoders/run_patch_experiment_v4_baseline_downsampled.py Jikji6 tripitaka \
       --n-patches 2,4,8,10 --ratio 0.18762

# §5 — 4-layer optimized decoder (28× headline)
python decoders/dna_decode_260328_0125.py --dataset Jikji6 --full
```

See [`decoders/README.md`](decoders/README.md) for the full per-variant narrative.

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
