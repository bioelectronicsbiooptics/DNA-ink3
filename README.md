# DNA-ink3

Code for the manuscript **"A Design-Guided Framework for Efficient and Scalable DNA Data Storage Decoding"** (Kang T, Song Y; submitted to *Briefings in Bioinformatics*).

The decoder assigns reads to oligo groups by the designed index region instead of full-read alignment, then reconstructs the payload with consensus, inner Reed–Solomon (RS) and outer XOR error correction. This repository contains the baseline decoder and the scripts for down-sampling and read-batch decoding.

## Files

| File | Used for | Description |
|---|---|---|
| `decoders/dna_decode_baseline_260507_Final.py` | Figures 1–2 | Baseline decoder: primer sorting, metadata detection, index sorting, consensus, inner RS, outer XOR, STL reconstruction. |
| `decoders/run_downsampling_experiment.py` | — | Runs the baseline pipeline on a given read list (`run_final_pipeline`); imported by the scripts below. |
| `decoders/run_downsampling_check.py` | Figure 3a | Full FASTQ versus ×5.33 random down-sampling (1/5.33 of reads drawn without replacement, seed 42). |
| `decoders/run_read_batch.py` | Figures 4–5 | Read-batch decoding: decodes the first N% of reads in FASTQ order (R1 then R2), without random sampling. |

## Requirements

Tested with Python 3.13 on macOS (Apple M5 Max).

```bash
pip install numpy scipy biopython galois trimesh pandas matplotlib
```

RS decoding uses `galois` (GF(2⁸), shortened RS(164,156)).

## Data

Sequencing data: NCBI BioProject [PRJNA1091804](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1091804). FASTQ files are not included in this repository.

## Usage

Baseline decoder (edit `DATA_ROOT` in `main()` to point to the FASTQ folder):

```bash
python decoders/dna_decode_baseline_260507_Final.py Jikji6
```

Down-sampling (Figure 3a):

```bash
python decoders/run_downsampling_check.py --name Jikji6 \
    --r1 Jikji6_1.fastq --r2 Jikji6_2.fastq --stl Jikji6_6.stl \
    --out downsampling_Jikji6.jsonl
```

Read-batch decoding (Figures 4–5). The first fraction must be 1.0; it sets the alignment to the original STL.

```bash
python decoders/run_read_batch.py --name Jikji6 \
    --r1 Jikji6_1.fastq --r2 Jikji6_2.fastq --stl Jikji6_6.stl \
    --fractions 1.0 0.1 0.01 0.004 --out read_batch_Jikji6.jsonl
```

Add `--reference designed_sequences.txt` (one designed oligo per line) to report consensus accuracy. A batch is counted as recovered when every facet of the original STL is matched within 10⁻³ coordinate units and no extra facets are present.

## License

MIT License (see `LICENSE`).

## Contact

Youngjun Song (corresponding author) — yjunsong.inu@stdbioelec.com
