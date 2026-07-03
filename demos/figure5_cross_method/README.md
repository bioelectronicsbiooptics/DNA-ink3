# Figure 5a — Cross-method read-grouping substitution demos

Reproducible scripts for **Figure 5a** of the iSCAN manuscript: substituting
iSCAN's positional index read for the **content-based read-grouping stage**
of representative external DNA-storage codecs, measured on **real sequencing
reads**.

The core claim these scripts support: for codecs that recover oligos by
content-based clustering/alignment over a **fixed positional index**, iSCAN's
positional read replaces that grouping step and accelerates it by one to two
orders of magnitude at ≥ 99 % oligo recovery — **without a fixed-offset
assumption**, because iSCAN *searches* for a conserved landmark (anchor /
primer) exactly as the core pipeline scans for its KIS poly-A anchor.

---

## Scope (read before interpreting the numbers)

- **Only the read-grouping stage is substituted.** The shared consensus and
  error-correction (RS + XOR) stages are *not* re-run; "recovery" here means
  the fraction of oligos recovered into stable read groups, **not** a full
  file decode.
- Consequently the **end-to-end** gain is bounded by the shared consensus
  stage (≈ 0.9×–4.6×) and is **not** equal to the grouping-stage speedup.
  Do not quote a single end-to-end multiplier.
- The grouping-stage speedup is robust; the exact multiplier varies with
  (a) the baseline (in-memory LSH vs an external aligner process),
  (b) the iSCAN variant (anchor-search / nearest / exact-slice), and
  (c) data scale.

---

## Scripts

| Script | Codec | What it demonstrates |
|---|---|---|
| `antkowiak_realreads_demo.py` | Antkowiak 2020 | **Real Illumina reads**: `ACAAC` anchor search → 7-nt index. ≈ 10× grouping speedup, 98.7 % oligo recovery (89.9 % anchor detection). |
| `antkowiak_demo.py` | Antkowiak 2020 | Real oligos + the authors' own LSH code as baseline (+ consensus over-partition analysis). |
| `yazdi_demo.py` | Yazdi 2017 | **Real MinION reads** (1,360 reads, ground-truth labels): BWA-MEM baseline vs 5′-address nearest-match. ≈ 8×, 100 % (17/17 oligos). |
| `yinyang_realreads_demo.py` | Ping 2022 (Yin-Yang) | **Real DNBSEQ reads** (CNGB CNP0001650): 20-nt primer PrimerSort (96.5 % detection). ≈ 124×, ≈ 100 % recovery. |
| `yinyang_demo.py` | Ping 2022 (Yin-Yang) | Real Yin-Yang oligos + generic LSH baseline (illustrative — the codec has no in-codec grouping step). |
| `organick_demo.py` | Organick 2018 | **Boundary case** — high-error clustered Nanopore reads with random cluster centres (no designed error-correcting index). iSCAN exact-slice is fast but recovery drops to 86 %; nearest-match reaches 99 % but is slow. Shows where the clean substitution does **not** hold. |

Headline three (Antkowiak, Yazdi, Yin-Yang) are the codecs plotted in
Figure 5a. Organick is documented as the honest boundary case.

---

## External data (public sources — not redistributed here)

The scripts read third-party data from the original studies. Obtain each
dataset from its public source and set the data path at the top of the
corresponding script:

| Codec | Data used | Public source |
|---|---|---|
| Antkowiak 2020 | Illumina reads (`File3`) + reference oligos (`File1_ODNA.txt`) | Reconstruction code + oligo references: `MLI-lab/noisy_dna_data_storage` (GitHub); raw reads from the study deposition. Nat Commun 2020;11:5345, doi:10.1038/s41467-020-19148-3 |
| Yazdi 2017 | MinION reads + reference sequences | `smhty/MATLAB_MinION` (GitHub). Sci Rep 2017;7:5011, doi:10.1038/s41598-017-05188-1 |
| Ping 2022 (Yin-Yang) | DNBSEQ reads (File 3) + encoded oligos | Reads: CNGB / CNSA accession **CNP0001650**; encoded oligos: `BGI-SynBio/YinYangCode` (GitHub). Nat Comput Sci 2022;2:234–242, doi:10.1038/s43588-022-00231-2 |
| Organick 2018 | Clustered Nanopore reads (`Centers.txt`, `Clusters.txt`) | `microsoft/clustered-nanopore-reads-dataset` (GitHub). Nat Biotechnol 2018;36:242–248, doi:10.1038/nbt.4079 |

> The raw sequencing datasets are large and belong to their original authors;
> only the substitution **scripts** are provided here.

### Configuring data paths

Scripts read their inputs from `$ISCAN_DEMO_DATA` (default: `./data`). Download
each dataset from the sources above and arrange it as:

```
$ISCAN_DEMO_DATA/
  antkowiak/File1_ODNA.txt          reference oligos (MLI-lab repo)
  antkowiak/File3.fq.gz             Illumina reads
  yazdi/ref.fasta                   17 reference codewords
  yazdi/reads.fasta                 MinION reads (ref{N}_read{M} labels)
  yinyang/encoded.dna               encoded oligos (BGI-SynBio repo)
  yinyang/DS_Sep_S2_F3_2_1.fq.gz    DNBSEQ reads (CNGB CNP0001650)
  organick/Centers.txt              cluster centres
  organick/Clusters.txt             clustered Nanopore reads
```

External tools `bwa` and MUSCLE are resolved from `PATH`; override with the
`BWA` and `MUSCLE` environment variables if needed.

---

## Dependencies

- Python ≥ 3.9 with `numpy`
- `bwa` (for the Yazdi BWA-MEM baseline)
- MUSCLE (for consensus in the Antkowiak/Yazdi over-partition analysis)

---

## Usage

```bash
export ISCAN_DEMO_DATA=/path/to/data      # location of the datasets above (default: ./data)

# Real-read validations (Figure 5a headline)
python3 antkowiak_realreads_demo.py [N]   # optional N = reads to sample
python3 yinyang_realreads_demo.py  [N]
python3 yazdi_demo.py                      # Yazdi is already real MinION reads

# Supporting / boundary
python3 antkowiak_demo.py                  # authors' LSH baseline + over-partition analysis
python3 yinyang_demo.py                    # illustrative LSH baseline
python3 organick_demo.py                   # boundary case (high-error Nanopore)
```

---

## Honesty notes

- **Antkowiak end-to-end** depends on over-partitioning: the authors' LSH
  fragments 500 true oligos into ~2,923 (≥ 5-read) clusters, whereas iSCAN
  yields ~513. Depending on whether over-partitioned clusters are all sent to
  consensus (~4.6×) or merged downstream to a 500 cap (~1.0×), the end-to-end
  gain spans ≈ 1.0×–4.6×. These demos measure the **grouping stage**, not the
  authors' full pipeline (which adds SSW cluster refinement); a single
  end-to-end multiplier is therefore **not** reported.
- **Yin-Yang baseline is illustrative**: the codec exposes no in-codec
  grouping step (the paper uses external assembly), so a generic LSH stands in
  as the baseline. Its index is the decoding-binary leading bits, not a
  DNA positional field.
- **Organick is a genuine limit**: high-error Nanopore reads with random
  cluster centres (no designed, hashable index) break the clean substitution.
- **CHN and OneJoin/DBSCAN** are substitutable in principle but their reference
  data are not publicly available, so they are not reproduced here.
