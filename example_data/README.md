# example_data/

Placeholder for a small reproducible slice (≈ 1,000 paired reads from the J6
Jikji 6 payload) plus the ground-truth STL mesh, used by the README smoke
test and by `Figure3/run_downsampling_experiment.py` in dry-run mode.

The full FASTQ collection (≈ 230 GB across 13 datasets) is hosted at NCBI
under BioProject **PRJNA1091804** and is **not** redistributed through this
repository.

## Expected layout once populated

```
example_data/
├── J6_R1.example.fastq          # 1,000 paired reads
├── J6_R2.example.fastq
└── stl/
    └── J6_reference.stl         # ground-truth mesh
```

## How to populate

Run the helper script (added in a follow-up commit), or extract the slice
manually:

```bash
seqtk sample -s42 Jikji6_R1.fastq 1000 > example_data/J6_R1.example.fastq
seqtk sample -s42 Jikji6_R2.fastq 1000 > example_data/J6_R2.example.fastq
```
