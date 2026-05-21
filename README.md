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

## Headline benchmarks

| Dataset | Payload | Oligos | MATLAB reference | iSCAN3 (4-layer) | Speed-up |
|---|---|---|---|---|---|
| Jikji 6 | 0.2 MB | 1,018 | 109 s | **3.9 s** | **28×** |
| Jikji 7 | — | 1,339 | 250 s | 17.5 s | 14× |
| Tripitaka Koreana | 16.6 MB | 92,626 | ~23 min | ~1–2 min | ≥10× |

Down-sampling (×5.33, mean coverage ×2,686 → ×504): J4 34.5 s → **2.47 s** (14×) at 100 % accuracy, NGS cost reduced ~80 %.

Patch-based parallel decoding: index match-rate 99.87 % at ×504 mean coverage; Tripitaka 92,626-oligo payload partitioned into 10 patches reduces wall-clock from ~23 min to ~3 min on a 10-core CPU.

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

The corresponding STL ground-truth meshes are bundled with this repository under `example_data/stl/`; the full FASTQ collection (≈ 230 GB) is hosted at NCBI.

---

## Project structure

```
DNA-ink3/
├── README.md                          # this file
├── LICENSE                            # MIT
├── pyproject.toml                     # Python dependencies + project metadata
├── uv.lock                            # Pinned dependency tree (uv)
│
├── python_pipeline/                   # Core decoder module (Python 3.14)
│   ├── dna_decode_260328_0125.py      # 28× headline reference (V14)
│   ├── dna_decode_260328_0132.py      # numba prange 16-core parallel variant
│   └── ...                            # incremental optimization history
│
├── Figure1/                           # Figure 1 — MATLAB baseline + Python conversion
│   ├── matlab/                        # MATLAB v240108 reference (15 .m files)
│   └── python/                        # Python translation of the baseline
│
├── Figure2/                           # Figure 2 — Data volume scaling
│   ├── run_volume_scaling.py          # 13-dataset driver
│   ├── aggregate_v8_11.py             # Per-dataset progress aggregator
│   └── fill_figure2_svgs_v2.py        # Filled SVG generator
│
├── Figure3/                           # Figure 3 — Downsampling efficiency
│   └── run_downsampling_experiment.py
│
├── Figure4/                           # Figure 4 — Patch-based decoding
│   └── run_patch_experiment.py
│
├── Figure5/                           # Figure 5 — Scalable optimization (4-layer stack)
│   └── dna_decode_*.py                # Per-layer ablation variants
│
├── Figure6/                           # Figure 6 — Architectural generalization
│   └── papers/                        # 16 external reference encodings (PDFs)
│
├── GPU_Bench/                         # Optional MLX vs. CPU benchmark (Apple Silicon)
├── GPU_RS_Prototype/                  # GPU-accelerated Reed–Solomon syndrome prototype
│
├── matlab_ref/                        # MATLAB R2022a reference (v240108)
│   └── Single_Decoding_Time.m
│
└── example_data/                      # Small reproducible slice for CI / tutorials
    ├── J6_R1.example.fastq            # 1,000 paired reads from Jikji 6
    ├── J6_R2.example.fastq
    └── stl/J6_reference.stl           # Ground-truth STL mesh
```

---

## Requirements

- **Python** ≥ 3.10 (3.14 recommended; tested on macOS arm64 / Linux x86_64)
- **MATLAB** ≥ R2022a (only required to regenerate the wet-lab FASTQ used in §1 reference)
- **GCC / Clang** with C99 support (for Numba native compilation)
- *(Optional)* **Apple Silicon** with MLX 0.31.2 for the GPU RS prototype

Core Python dependencies (full list pinned in `pyproject.toml` / `uv.lock`):

| Package | Version | Role |
|---|---|---|
| numpy | 2.3.4 | Whole-array vectorization (Layer 1) |
| pandas | 2.3.3 | Per-dataset aggregation tables |
| matplotlib | 3.10.7 | Figure rendering (Agg backend) |
| biopython | 1.85 | FASTQ I/O |
| reedsolo | latest | GF(2⁸) lookup-table Reed–Solomon (Layer 3) |
| numba | 0.65.0 | JIT consensus loop (Layer 4) |
| scipy | 1.16.3 | Statistical helpers |
| trimesh | 4.11.5 | STL re-export and mesh diff |
| xlsxwriter | 3.2.9 | Per-stage timing tables |

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

A first import of `dna_decode_260328_0125` will trigger the one-time Numba JIT compile (≈ 8 s) and a 64 KB GF(2⁸) lookup-table precompute (instantaneous); both artifacts are cached on disk and re-used on every subsequent run.

---

## Usage

### 1. End-to-end decode of the example payload

```bash
python python_pipeline/dna_decode_260328_0125.py \
    --r1 example_data/J6_R1.example.fastq \
    --r2 example_data/J6_R2.example.fastq \
    --primer-set P6 \
    --output decoded_J6.stl
```

Expected output on the M2 Studio reference machine:

```
[Stage 1 — Primer sorting]     0.04 s
[Stage 2 — KIS detection]      0.01 s
[Stage 3 — Index clustering]   0.03 s
[Stage 4 — Consensus + RS]     0.02 s
[Total wall-clock]             0.10 s   (example slice, full J6: 3.9 s)
[Reconstruction fidelity]      100 %    (bit-for-bit identical to ground-truth STL)
```

### 2. Reproduce the V14 28× headline benchmark

```bash
python Figure5/dna_decode_260328_0125.py --dataset Jikji6 --full
```

### 3. Down-sampling sweep (Figure 3)

```bash
python Figure3/run_downsampling_experiment.py --rates 1.0,0.50,0.25,0.188,0.10
```

### 4. Patch-based parallel decode (Figure 4)

```bash
python Figure4/run_patch_experiment.py \
    --dataset tripitaka --patches 10 --workers 10
```

### 5. Volume-scaling sweep across all 13 datasets (Figure 2)

```bash
export DNA_DATA_DIR=/path/to/fastq      # required: directory containing all *_R1/R2.fastq
python Figure2/run_volume_scaling.py
python Figure2/aggregate_v8_11.py       # generates v8_11_results.xlsx
```

### 6. MATLAB reference cross-check

```matlab
% From the matlab_ref/ directory in MATLAB R2022a
addpath('matlab_ref')
run_smoke_all                            % 13-dataset reference timing
```

---

## Reproducing each figure

| Manuscript figure | Script | Output |
|---|---|---|
| Figure 1 (iSCAN pipeline diagram + 360× headline) | `Figure1/python/dna_decode_baseline.py` | `Figure1/Outputs/` |
| Figure 2 (volume scaling, 13 datasets) | `Figure2/run_volume_scaling.py` + `aggregate_v8_11.py` + `fill_figure2_svgs_v2.py` | `Manuscript_Drafts_v2/Figure_2/*.svg` |
| Figure 3 (down-sampling efficiency) | `Figure3/run_downsampling_experiment.py` | `Figure3/Outputs/` |
| Figure 4 (patch-based parallel) | `Figure4/run_patch_experiment.py` | `Figure4/Outputs/` |
| Figure 5 (4-layer optimization stack) | `Figure5/dna_decode_260328_0125.py` | `Figure5/Outputs/` |
| Figure 6 (21-method architectural generalization) | `Figure6/Figure_6a_Option1_table.svg` (manual edit) | `Figure6/Figure_6a_Option1_table.png` |

Per-figure SVG fillers parse `*_results.xlsx` and substitute placeholder text directly in the SVG, so re-running the analysis pipeline plus the corresponding filler is sufficient to regenerate every published panel without leaving Python.

---

## Verification and tests

- A wet-lab cross-check is provided as `matlab_ref/Single_Decoding_Time.m`; Python and MATLAB pipelines agree on all 13 public payloads to bit-for-bit STL identity (last verified 2026-05-07).
- A minimal smoke test (`example_data/`, ≈ 1,000 reads, runs in < 1 s) is exercised on every commit.
- Coverage and dropout integrity for the down-sampling sweep are emitted to JSON for inspection (`Figure3/Outputs/downsampling_results/*.json`).

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
