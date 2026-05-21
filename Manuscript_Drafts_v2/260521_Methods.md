# Materials and Methods

> **iSCAN3 manuscript — Methods section (V0.1)**
> 작성일: 2026-05-21
> 출처: (1) ACS Nano 2026 (선행 wet-lab) — overhang PCR 및 VR 섹션 제외, 변형 사용 / (2) BIB bbag224 — inline parameter 명시 스타일 + Data Availability GitHub 분리 적용

---

## Materials

### 3D model source and pre-processing

Three-dimensional models of cultural heritage objects were obtained from the Cultural Heritage Administration of the Republic of Korea and public 3D model repositories (Supplementary Note 1). Each STL mesh was inspected and trimmed using ZBrush (Maxon Computer, CA, USA) to remove non-manifold faces and isolated vertices prior to digital-to-DNA translation. The trimmed STL files were converted to oligonucleotide pools through the in-house *3D DTL* (3D Data Translation Layer) encoder, executed on MATLAB R2022a (MathWorks, MA, USA). The same MATLAB toolchain was used downstream for FASTQ generation and decoding cross-checks throughout this work.

### Reagents and primer sets

Primer sets used to define seven independently addressable DNA pools (P1–P7, paired with payload pools J1–J7) were adapted from a prior random-access study and procured from Integrated DNA Technologies (Coralville, IA, USA). The encoded payload pools themselves were ordered as SurePrint custom oligonucleotide pools (Agilent Technologies, Santa Clara, CA, USA). Stock DNA pools were resuspended in nuclease-free water and quantified by Qubit dsDNA HS assay on a Qubit 4 Fluorometer (Thermo Fisher Scientific, Waltham, MA, USA).

---

## Wet-lab procedures

### Random-access PCR amplification

To selectively retrieve a single payload from the multiplexed pool, random-access PCR was carried out per target index. Each 20-µL reaction contained 0.424 ng of the multiplexed J1–J7 synthetic pool, 0.5 µM of the corresponding paired primer set (one of P1–P7), and KAPA HiFi HotStart ReadyMix (Roche, Basel, Switzerland). The thermocycling profile was: initial denaturation at 98 °C / 30 s, followed by 40 cycles of 98 °C / 20 s, 58 °C / 20 s, 72 °C / 15 s, and a final hold at 4 °C. Amplification specificity was verified by two-step real-time PCR (RT-PCR) on a QuantStudio 5 (Thermo Fisher Scientific) with KAPA SYBR FAST qPCR Master Mix; the program comprised 40 cycles of 95 °C / 15 s, 58 °C / 15 s, and 72 °C / 40 s. Ct values were analyzed in QuantStudio Design and Analysis Desktop Software v1.5.2.

### Gel-based size selection and PCR clean-up

PCR products were resolved on 2 % agarose gels (Sigma-Aldrich, St. Louis, MO, USA) containing SYBR Gold nucleic acid stain (Thermo Fisher Scientific). Bands of the expected length were visualized under LED transillumination (KA.33-61, Korea Ace Scientific, Republic of Korea), excised, and purified using the NucleoSpin Gel and PCR Clean-up Kit (Macherey-Nagel, Düren, Germany) according to the manufacturer's instructions. Purified DNA was re-quantified on the Qubit 4 prior to any downstream step.

### DNA-ink formulation and stamping

DNA ink was prepared by mixing purified DNA amplicons with traditional Indian ink (Korean Mukwoondang, Seoul, Republic of Korea) so that the final composition was 2.709 ng/µL DNA in 10 % (v/v) Indian ink; the optimal DNA-to-ink ratio was determined empirically by stamping titrations (data summarized in Figure 6). For physical stamping, scaled aluminum replicas of the Jikji metal type J1–J7 were fabricated on an SLM 3D printer 280 2.0 (SLM Solutions Group AG, Lübeck, Germany), inked, and pressed onto Whatman filter paper (GE Healthcare Korea, Seoul, Republic of Korea). Stamping concentration for the ink-on-paper readout was independently optimized to 2.35 ng/µL DNA in a 10 µL ink suspension dispensed within a 2 µL printing volume.

### Stamped-pattern imaging

Stamped DNA-ink patterns were characterized at multiple length scales. Macro-scale fluorescence was acquired under LED transillumination (KA.33-61, Korea Ace Scientific). Cellular-scale fluorescence imaging used a Leica DM RBE upright microscope coupled to a Prime 95B sCMOS camera (Teledyne Technologies, Thousand Oaks, CA, USA). Nano-scale morphology and elemental distribution were obtained on a JSM-7800F field-emission scanning electron microscope (JEOL, Akishima, Tokyo, Japan) equipped with EDS. For DNA recovery from stamped patterns, the inked paper was cut, eluted in TE buffer (10 mM Tris-HCl, 1 mM EDTA, pH 8.0), and incubated for 30 min at 30 °C. The eluate was processed with the QIAquick PCR Purification Kit (Qiagen, Hilden, Germany) and re-amplified using the random-access PCR program above.

### Low-copy DNA-ink preparation

To test ultra-low-input recovery, twelve independent dilution samples (A–L) were prepared from the J6 payload after random-access amplification with primer pair P6. Samples A–J were each prepared by 10-fold serial dilution from the upstream stock; samples K and L were prepared by 4-fold dilution. Samples A and B were quantified directly on the Qubit 4 (212 pg/µL and 22.4 pg/µL, respectively). Quantification of the remaining ten samples (C–L), whose concentrations were below the Qubit linear range, was performed by qPCR against a standard curve generated from A and B, with copy number computed using the DNA Copy Number and Dilution Calculator (Thermo Fisher Scientific). All twelve samples were sequenced in six iSeq 100 runs, two samples per run, after independent library preparation. Coverage and recovery for every fabrication step (original stock, first PCR, first clean-up, ink mixing, paper elution, second PCR, second clean-up) were tracked by Qubit dsDNA HS and reported in Figure 6h.

### 3D printing of decoded cultural heritage replicas

To physically demonstrate that decoded payloads could be re-instantiated, decoded STL files were printed in diverse filaments. The Golden Gate miniature used red polylactic acid (PLA) filament (Four Way Inc., Republic of Korea); the Venus Statue miniature used marble-filled PLA (SUNLU Inc., Zhuhai, China); the Tripitaka Koreana replica was printed at native scale on PLA wood filament (SUNLU Inc.) using both an Ender-3 Pro (Creality 3D Technology, Shenzhen, China) and equipment provided by Mankato Graphics Inc. (Incheon, Republic of Korea). Aluminum replicas of the Jikji metal type were printed on the SLM 3D printer 280 2.0 (SLM Solutions Group AG). To showcase material versatility, additional Jikji replicas were printed in white-rubber ABS (eSUN, Shenzhen, China), PLA (Four Way Inc.), and wood filament (SUNLU Inc.) on an X-Marker 3D printer (QIDI Tech, Wenzhou, China).

### Sequencing library preparation and Illumina iSeq 100 sequencing

A modified adapter-ligation protocol — derived from Illumina TruSeq Nano and TruSeq ChIP — was used throughout. PCR-amplified pools were first end-repaired with End Repair buffer (ERP2) and cleaned with AMPure XP beads (Beckman Coulter, Brea, CA, USA). A single-base 'A' overhang was added by A-tailing ligase, after which TruSeq DNA CD Index adapters (Illumina, San Diego, CA, USA) were ligated, each sample receiving a unique index combination. Adapter-ligated 25-µL reactions were cleaned with Illumina purification beads, enriched in a 5-µL PCR Primer Cocktail + 20-µL Enhanced PCR Mix, and amplified for 8 PCR cycles. Final library concentration was measured by Qubit dsDNA HS. Libraries were diluted to 50 pM, spiked with 10–20 % PhiX, and run on the Illumina iSeq 100 platform under the 150-nt paired-end protocol per the iSeq 100 System Denature and Dilute Library Guide. Selected libraries were additionally sequenced on the Illumina HiSeq 2500 platform via Macrogen Inc. (Seoul, Republic of Korea), again under the 150-nt paired-end protocol. FASTQ files used for decoding were generated by an in-house MATLAB R2022a pipeline.

---

## Computational methods

The remainder of this section describes the iSCAN decoding pipeline and the four-layer algorithmic optimization stack analyzed in §1–§5. Implementation follows the structure popularized by recent bioinformatics method papers: all hyperparameters, dimensions, and threshold values are reported inline; the corresponding source code is released through a public GitHub repository (see *Data and Code Availability*).

### Software and hardware environment

All computational experiments were executed on Apple M2 Studio workstations (Apple Inc., Cupertino, CA, USA) running macOS 26.3. The decoder is implemented in Python 3.14.3. Numerical routines depend on NumPy 2.3.4, SciPy 1.16.3, and Biopython 1.85 for FASTQ I/O; Reed–Solomon arithmetic uses the `reedsolo` package configured with the GF(2⁸) primitive polynomial `0x11d`; mode-based consensus is JIT-compiled through Numba 0.65.0; STL re-export uses trimesh 4.11.5. Plotting and report rendering use Matplotlib 3.10.7 (Agg backend) and pandas 2.3.3. A MATLAB R2022a reference implementation (see `matlab_ref/Single_Decoding_Time.m`, v240108) is included to validate numerical parity. Random number generation, where used, was seeded at 42 for reproducibility.

### iSCAN positional decoding pipeline

The iSCAN decoder operates in four sequential stages, each restricted to the index-bearing positional regions of the read rather than performing global alignment.

**Stage 1 — Primer sorting.** Paired-end reads (R1, R2, 150 nt each) are tested against the forward primer over positions 1–20 nt and the reverse primer over positions 131–150 nt. A read pair is retained only when both primers achieve an exact match of ≥ 18 contiguous nucleotides. Reverse reads are reverse-complemented prior to all downstream comparisons so that subsequent stages can operate on a directionally homogeneous read pool.

**Stage 2 — Key Information Sequence (KIS) detection.** A conserved poly-A marker (`AAAAA`) embedded at position 98 of each 150-nt read serves as a constant-time positional anchor. Detection is performed by a windowed binary mask over the byte-encoded read matrix (see Layer 1 below). Once anchored, the downstream three-dimensional matrix metadata — index length, vertex count, facet count, and null-strand count — is read out from the immediately adjacent positions; the operation completes in O(N) over the read count with a single pass and is independent of payload size.

**Stage 3 — Index clustering.** Reads are clustered exclusively on the basis of a localized 5–10 nt index scan immediately following the KIS anchor. Cluster assignment is performed by stable index sorting followed by group-boundary detection via `numpy.searchsorted`; no global pairwise alignment is performed. The patch-based variant of this stage is described in *Patch-based parallel decoding* below.

**Stage 4 — Consensus and error correction.** Within each index cluster, the consensus payload is computed column-wise by majority mode over the uint8 codeword matrix (Layer 4). Residual per-oligo substitutions are repaired by inner Reed–Solomon decoding on the consensus; oligos that the inner decoder cannot recover are explicitly *dropped* and reconstructed by the outer XOR layer (see *Down-sampling and recovery strategy*).

### Four-layer algorithmic optimization stack

The pipeline above is implemented through a stack of four orthogonal optimization layers, each targeting a distinct bottleneck identified during stage-wise profiling.

**Layer 1 — Whole-array vectorization.** All per-element loops over reads, nucleotides, and primers are reformulated as whole-array NumPy operations executed in the compiled BLAS/SIMD backend. Primer matching becomes a broadcast matrix equality test over the byte-encoded read matrix; KIS scanning becomes a strided slice over the same buffer; cluster assignment uses `numpy.argsort` + `numpy.searchsorted` for group-boundary detection. Per-stage measured speed-ups under this layer alone are 11× (primer sorting), 122× (KIS scan), 73× (consensus), 10× (Evert direction correction), and 16× (null-strand detection).

**Layer 2 — uint8 byte-level encoding.** Each FASTQ read is converted at ingestion from ASCII text into an N × 150 `uint8` matrix through `numpy.frombuffer`, after which all downstream stages operate on dense integer buffers. The mapping {A, C, G, T} ↔ {0, 1, 2, 3} is implemented as a 256-entry lookup table indexed directly by the ASCII byte value. This conversion is lossless and eliminates the repeated character-to-integer translation that otherwise dominates the binarization stage; the measured per-stage speed-up is 71× for binarization.

**Layer 3 — Lookup-table GF(2⁸) Reed–Solomon arithmetic.** Reed–Solomon decoding requires multiplication and inversion in GF(2⁸) under the primitive polynomial `0x11d`. Rather than computing each field multiply via bitwise polynomial multiplication, the full 256 × 256 GF(2⁸) multiplication table (64 KB) and the corresponding log / antilog tables are precomputed once at program start. Each subsequent field multiplication then reduces to a single 2-D array lookup, while field addition remains a single bitwise XOR. The complete decoder is supplied by the `reedsolo` Python package configured at `n = 255, k = 223` with shortening to match the per-oligo payload length; batch decoding over all clusters is dispatched as a single Python call. The measured speed-up over the previous galois-based implementation is 175×.

**Layer 4 — Just-in-time native compilation of the consensus hot loop.** After Layers 1–3 are applied, the residual bottleneck is the column-wise majority-mode computation over the uint8 codeword matrix during consensus; for the Tripitaka-scale dataset, approximately 5 × 10⁸ individual reductions are performed. This function is annotated with `@numba.njit(parallel=True)` and is specialized at first call to the observed `uint8` dtype and 2-D shape, after which it is cached and re-used. The measured incremental speed-up is 53× for binarization-stage consensus.

**Aggregate.** With all four layers active, the end-to-end wall-clock for the J6 Jikji 6 payload decreases from 109 s (MATLAB reference) to 3.9 s (28×), and for the J7 payload from 250 s to 17.5 s (14×); the Tripitaka payload (92,626 oligos, ~22 GB FASTQ) decreases from ~23 min to approximately 1–2 min. 100 % bit-for-bit reconstruction is preserved at every operating point. Per-dataset wall-clock timings are tabulated in Supplementary Table S2.

### Patch-based parallel decoding

For payloads larger than ~10⁴ oligos, the cluster matrix is partitioned along the index axis into discrete patches. Patch boundaries are chosen using the upstream KIS index map so that each patch contains a contiguous range of oligo indices, and adjacent patches share a ± 100-index overlap region. The overlap is required because outer XOR parity may span the boundary, and its width was determined empirically as the smallest value at which boundary RS / XOR recovery did not degrade below the no-partition baseline. Default patch granularity for the experiments reported here is 10 patches per Tripitaka-scale dataset (≈ 9,263 oligos per patch). Patches are dispatched to worker processes through Python's `multiprocessing.Pool`, one worker per physical CPU core, with patch results merged in the parent process. Index purity at the operating point of ×504 average coverage is 99.87 %, measured as the proportion of reads whose post-merge index assignment matches the ground-truth index.

### Down-sampling and recovery strategy

Down-sampling is applied uniformly at the FASTQ level prior to any decoding stage. Reads are retained independently with probability `p = 1 / k`, where `k` is the down-sampling factor; for the headline experiment, `k = 5.33` corresponds to a transition from ×2,686 to ×504 mean per-strand coverage. After down-sampling, the four-stage iSCAN pipeline runs unchanged; the integrated inner-outer error-correction layer absorbs the additional sparsity. Specifically, for each cluster the inner Reed–Solomon decoder reports a pass/fail flag based on its syndrome check; clusters that fail and whose per-strand coverage falls below the dropout threshold of ×16 are *removed* from the consensus matrix and flagged as missing oligos. Missing oligos are then reconstructed by the outer XOR layer using the parity structure embedded at encoding time, which can recover up to *m* lost oligos per XOR group given *m* surviving parity strands.

### Performance benchmarking

End-to-end timing was measured by wrapping each pipeline stage with `time.perf_counter()` calls; reported numbers are means over five independent runs with cold OS file-system cache (cleared between runs) on the workstation described above. Reconstruction fidelity was checked at three layers: (1) bitwise comparison of the decoded payload against the original encoded binary, (2) STL mesh-level equality of vertex and facet tables after re-export with trimesh, and (3) NCBI-listed dataset hashes for the 13 public payloads (PRJNA1091804). Where present, deviations were re-traced to the input FASTQ and inspected manually; no deviations remained after all four optimization layers were active.

---

## Data and Code Availability

**Sequencing data.** All raw FASTQ files generated and analyzed in this study are deposited at the NCBI Sequence Read Archive under BioProject **PRJNA1091804** (3D Culture Heritage DNA data storage), covering the thirteen payloads Jikji 1–7, Seokgatap, Tower of Pisa, Statue of Liberty, Venus Statue, Golden Gate, and Tripitaka Koreana, plus the J6 overhang-chemistry variant J6_OH (SRR33084752).

**Decoded artifacts.** Decoded STL models, depth and perfect-match CSV reports, and the run-level FASTQ-to-mesh provenance log are bundled with the source release described below.

**Source code.** All code used to produce the figures and tables in this manuscript — the iSCAN decoder, the four-layer optimization stack, the patch-based parallel driver, the down-sampling and dropout / XOR-recovery routines, the MATLAB reference implementation, and the benchmarking harness — is openly available at:

> **https://github.com/bioelectronicsbiooptics/DNA-ink3**

The repository is released under the MIT License and includes a `README.md` with environment-setup instructions (`pyproject.toml`, `uv.lock`), an `example_data/` directory containing a 1,000-read FASTQ slice and the expected decoded STL, and pre-frozen pip / uv lockfiles for byte-level reproducibility. Specific commit hashes used to produce each figure are recorded in `Supplementary Note S3`.

---

## Notes on this draft

- Wet-lab section is transformed from ACS Nano 2026 *Materials and Methods*; the **Overhang PCR with Long Primers** subsection and the **Online Upload and VR** subsection have been intentionally excluded as instructed.
- Inline-parameter style mirrors *Briefings in Bioinformatics* article bbag224 (BIB 27(3), May 2026), with code referenced through a public GitHub repository and not reproduced in-text.
- All numerical values reported here are *placeholders that match the V15 manuscript narrative*. Final figures may shift slightly once the V16 Python-baseline qualifier is applied to the §2 timing claims; see `260518_Project_Status_Summary.md §5.A` for the open discrepancy.
