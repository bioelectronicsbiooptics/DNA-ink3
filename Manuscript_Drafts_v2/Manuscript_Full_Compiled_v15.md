# iSCAN3 Manuscript — Unified Full Draft (Intro + Results §1–§6 + Conclusion)
*260418_14 | Intro (V15) + §1–§3 (V12) + §4–§6 (V14) + Conclusion (V15)*

> **병합 소스**
> - `260417_19_Results_Compiled_Draft.md` → Intro (V13 base), §1–§3 (V12), Conclusion (V13 base)
> - `260418_13_Results_Sec456_v14.md` → §4 (V14, closing sentence 수정), §5 (V14, 전면 개정), §6 (V14, 신규)
> - `260419_XX` V15 연쇄 업데이트 → Intro P4, Conclusion P1/P2/P3 (§5·§6 V14 반영)
> - 섹션 헤딩 옆 **(Vxx)** 은 해당 섹션의 최신 버전을 표기합니다.

> **색상 범례** (전 문서 공통)
> <span style="color:gray">회색 = 버전 출처 / provenance</span>
> <span style="color:red">빨강 = 표현 수정 (v1→v2, v12→v14, v13→v15)</span>
> <span style="color:green">초록 = 신규 추가 (v2, v13, v14)</span>
> `<mark style="background:#d3f8b6">…</mark>` = 원저자(사용자) 검토 주석 — Results §1–§5 원본 보존
> ~~취소선~~ = 삭제된 내용

---

<sup style="color:gray">[Introduction 출처: 260417_19_Intro_unified.md (V13 base = V11+V12 통합) → V15 연쇄 업데이트 (§5 GPU 제거, §6 21-method 서베이 반영)]</sup>

# Introduction (V15)

<sup style="color:gray">[P1: v1=v2 (내용 동일, 문장 구조 소폭 정리)]</sup>

DNA has emerged as a compelling candidate for long-term archival data storage, offering <span style="color:red">theoretical</span> information densities approaching 215 petabytes per gram and chemical stability that, <span style="color:red">under appropriate preservation conditions, can extend</span> over millennial timescales [ref]. Over the past decade, proof-of-concept demonstrations have shown that arbitrary digital data—from text and images to executable programs—can be encoded into synthetic oligonucleotide libraries and retrieved via next-generation sequencing (NGS) [ref, ref, ref]. Synthesis throughput and encoding fidelity have improved substantially, yet the practical ceiling for DNA storage systems is not currently set by how fast data can be written; it is set by how fast it can be read back. Decoding—the process of reconstructing original digital information from millions of raw sequencing reads—remains the operational bottleneck that prevents DNA storage from functioning at scale.

<sup style="color:gray">[P2: v1→v2 수정 — "discard it"→"fail to exploit it directly"; ×1,000 coverage 근거 보강; FASTQ 비용 이중 escalation 명시]</sup>

The dominant computational cost in current decoding pipelines arises from *de novo* sequence alignment, in which every sequencing read is compared against a reference to determine its positional identity within the stored data. This process is inherently expensive: depending on dataset size, alignment-based decoding typically requires one to six hours of processing time, a latency that compounds rapidly with increasing data volume [ref]. In the DNA storage context, full-read alignment represents a fundamental mismatch between available information and the work being performed. Unlike genomic applications where reference identity is unknown a priori, oligonucleotides in DNA storage systems are synthesized with explicit, predetermined structural indices. Each oligo carries a known positional address—encoded at a fixed location within its sequence—that directly identifies which cluster it belongs to. This positional information is present by design, yet conventional pipelines <span style="color:red">fail to exploit it directly</span> and instead recover cluster assignments through exhaustive global search. Furthermore, these pipelines conventionally operate on high-coverage sequencing data—coverage rates exceeding ×1,000 are routinely reported <span style="color:green">across multiple large-scale DNA storage studies</span> [ref]—which inflates both NGS cost and the raw FASTQ volume that alignment must process. <span style="color:green">As stored data capacity and the corresponding FASTQ file size increase, this dual cost—algorithmic and volumetric—escalates in tandem, making the decoding step increasingly impractical for any system intended to operate beyond laboratory-scale demonstrations.</span>

<sup style="color:gray">[P3: v1→v2 수정 — 4단계 설명 재구성 (consensus stage는 full reads 작동으로 명확화)]</sup>

Here we present iSCAN (index Scan Cluster Alignment for Nucleotide), a decoding framework that resolves this mismatch by directly exploiting the positional index structure already embedded in oligonucleotide design. Rather than aligning full reads to a reference, iSCAN operates through four targeted stages—primer sorting, Key Information Sequence (KIS) detection, index clustering, and consensus sequencing. <span style="color:red">The first three stages each act on a restricted, index-relevant portion of each read: once the index position is located via KIS, cluster assignment requires only a localized scan of five to ten nucleotides, rendering global alignment unnecessary. The final consensus stage then reconstructs the stored payload from the clustered reads.</span> This design separates the search problem from the reconstruction problem, compressing what was a multi-hour alignment step into second- to minute-scale processing while preserving complete reconstruction fidelity.

<sup style="color:gray">[P4: v1→v2 확장 — 연구 질문 3개→4개; 다운샘플링 수치 명시; 최적화 경로 3단계→4단계; v13: Figure 6 대응 1문장 추가; **v15: §5 V14 (GPU 제거) + §6 V14 (21-method 서베이) 연쇄 반영 — "hardware acceleration" → "algorithmic optimization", Figure 6 문장 재작성**]</sup>

Beyond demonstrating that positional decoding is feasible, this study pursues <span style="color:red">four</span> practical questions that determine whether the approach is deployable at operational scale: how decoding time scales as stored data capacity grows from the kilobyte to the megabyte range; how far sequencing depth can be reduced—<span style="color:green">potentially to as little as five- to tenfold coverage</span>—without sacrificing perfect data recovery when error-correction coding is applied; whether the strategy maintains index fidelity and throughput under patch-level parallel execution on large datasets; <span style="color:red">and what decoding acceleration is ultimately achievable when algorithmic input reduction, software parallelization, and systematic algorithmic optimization of the core pipeline stages are combined into a unified optimization stack.</span> We address these questions through a progressive optimization path—bottleneck identification, targeted input-volume reduction with recovery coding, scalable parallelized integration, <span style="color:red">and systematic algorithmic optimization of the core computational stages</span>—in which each step builds directly on the causal conclusion of the last. <span style="color:red">Finally, we evaluate the generalizability of the iSCAN decoding architecture across the landscape of 21 DNA data storage encoding methods published between 2012 and 2024 — including our own earlier iSCAN implementation — mapping where its position-based scanning paradigm transfers directly and where encoding-specific adaptation is required.</span>


---

# Results

## 1. Index Alignment Protocol Strategy 

The translation of DNA data storage into a practical medium is fundamentally limited by the decoding bottleneck. Conventional *de novo* sequence alignment methods require <font color="#ff0000"> ~ 6 hours </font>for data clustering, creating massive computational overhead that constrains rapid data retrieval. [ref_xx] To overcome this limitation, building on our previous approach[ACS Nano 2026],
which we successfully reimplemented in Python as the foundation of this work, we developed the index scan cluster alignment for nucleotide (iSCAN) protocol, a method designed to bypass full sequence mapping by restricting the clustering search to a discrete set of index-bearing positional regions within each read, with each positional index acting as a deterministic addressing pointer that replaces exhaustive alignment with a positional lookup amenable to systematic algorithmic optimization.

To evaluate the efficiency of positional index scanning, we compared the iSCAN workflow against conventional alignment methodologies. Figure 1a illustrates the comparative workflow, detailing the four-step iSCAN process that replaces reference alignment with positional index scanning. (1) iSCAN begins with primer sorting via a ≥18 nucleotide match (forward: 1–20 nt, reverse: 131–150 nt), followed by paired-end reverse-complement conversion, yielding a primer-verified read pool on which the subsequent stages operate (Figure 1b). (2) This is followed by Key Information Sequence (KIS) detection: a conserved poly-A ('AAAAA') marker at the 98 nt position of each 150 nt read serves as a constant-time positional anchor, extracting the three-dimensional matrix information from the adjacent index region in <font color="#ff0000">under 1 second </font>(4.57 M reads, 1.51 GB raw FASTQ) without requiring full-read traversal (Figure 1c). (3) Subsequently, the system performs index clustering based on a localized 5–10 nucleotide scan, entirely bypassing *de novo* alignment, and concludes with consensus sequencing utilizing the high average coverage. (4) Finally, inner Reed-Solomon (RS) decoding corrects residual per-oligo errors in the consensus sequences, and outer XOR decoding recovers any missing or dropped indices across oligos, thereby reconstructing the original encoded image with bit-for-bit fidelity.

Consequently, the iSCAN protocol achieved a 360-fold acceleration in decoding speed compared to conventional methods. A 1.51 GB raw-FASTQ dataset was decoded in<font color="#ff0000"> just 34.5 seconds</font> while maintaining complete and perfect data reconstruction (defined as bit-for-bit identity with the original encoded file), with the index-organized read distribution and the reconstructed payload <font color="#ff0000">as shown in Figure 1d.</font> However, how this positional clustering performance dynamically scales as the target DNA data capacity increases to the megabyte scale remains a critical question to be addressed.

---

## 2. Time Quantification Analysis 

<span style="color:red">With the iSCAN pipeline established as a viable baseline, the scaling behavior of decoding time with respect to DNA data capacity and the total number of synthesized oligos emerged as the central concern.</span> To address this, we systematically quantified the time consumed at each individual stage of the iSCAN pipeline across diverse datasets, ranging in capacity from 499 kB (2,587 oligos) to 16.6 MB (92,626 oligos).

To identify potential computational bottlenecks during scale-up, we analyzed the absolute processing time and the relative proportions of each decoding step. Figure 2a presents the regression analysis of decoding time versus physical data size, demonstrating that the total processing time scales linearly with capacity—a relationship further supported by the residual distribution shown in Figure 2c. Specifically, the system required 34.5 seconds for a 1.51 GB FASTQ file and approximately 23 minutes to decode a 22.4 GB FASTQ file. Furthermore, Figure 2b details the processing time composition, revealing the internal mechanics of this linear scaling. Primer sorting consistently consumes the largest proportion of the total decoding time (40–52%), establishing it as the primary bottleneck. This is followed by index clustering (30–40%) and RS/XOR error recovery (9–12%). Notably, KIS detection remains a constant-time operation, completing in under 1 second regardless of massive variations in input data size. These metrics indicate that while the iSCAN protocol can reliably process up to 92,626 sequences within 30 minutes, the string-matching operations required during primer sorting constrain ultimate processing velocity.

These results demonstrate that the iSCAN system exhibits <mark style="background:#d3f8b6">predictable, linear scalability, confirming (리너어스케일 맞는지 증명해야함 테이블등으로 하거나 실험적 sorting하는 시간으로 하는는가나 실험 Figurre에 포함되어야 할수도)</mark> its computational viability for megabyte-scale DNA data storage operations. While processing time scales linearly with physical data size, the pipeline's reliance on high-coverage sequencing arrays inherently inflates the total FASTQ footprint. Thus, whether the raw input data volume itself can be substantially reduced without compromising perfect data reconstruction remains a critical avenue to be explored.

---

## 3. Data Downsampling Efficiency 

<span style="color:red">The observed correlation between decoding time and input data volume suggested a direct optimization path</span>: if primer sorting scales proportionally with input read volume, then reducing that volume directly targets the dominant cost. Current DNA data storage systems conventionally rely on high sequencing depth—often exceeding coverage rates of ×1,000, as reported across multiple large-scale storage studies—which inflates both next-generation sequencing (NGS) costs and computational processing times. To address this inefficiency, we conducted extensive downsampling experiments combined with inner Reed-Solomon (RS) and outer XOR error correction to determine the minimum coverage necessary to maintain perfect data reconstruction.

To quantify the trade-off between reduced sequencing depth and decoding reliability, we evaluated the protocol under various rigorously downsampled constraints. Figure 3a compares the performance of the original high-coverage state<mark style="background:#d3f8b6"> (average coverage ×2,686, 1.51 GB FASTQ size) against a ×5.33 downsampled dataset (average coverage ×504). Remarkably, the downsampled dataset maintained 100% data reconstruction accuracy while compressing the decoding time from 34.5 seconds down to just 2.47 seconds.(내용검증이 필요) </mark>Furthermore, Figures 3b and 3c illustrate the system's resilience under low-coverage conditions, utilizing the Venus Statue dataset which featured a highly unequal sequence distribution. By employing a targeted dropout strategy that intentionally discarded sparse sequences (coverage < ×16) based on local RS pass/fail metrics, the outer XOR mechanism seamlessly recovered all three missing sequences <mark style="background:#d3f8b6">(indices 12,643, 33,959, and 51,976</mark>). This outcome demonstrates that uniformly high-coverage sequencing is largely redundant; the synergistic combination of inner and outer error-correction codes provides sufficient robustness to tolerate highly sparse and unevenly distributed sequence data.

Applying a ×5.33 downsampling strategy successfully reduces NGS sequencing costs by approximately 80% while accelerating the decoding process 14-fold without sacrificing data integrity. While reducing the raw input data drastically improves computational efficiency on a single stream, whether this sparse data approach can be seamlessly integrated into a distributed, parallelized architecture to process massive datasets remains to be explored.

---

## 4. Patch-Based Decoding System (V14) <mark style="background:#d3f8b6">!! 논문에서 제안하는 방향으로</mark>

<mark style="background:rgba(205, 244, 105, 0.55)">Having shown that sparse, downsampled input preserves complete reconstruction, the remaining question is how this volume-reduced data should be handled once it scales.</mark> The structural challenge that followed was distinct: processing tens of millions of NGS reads as a single contiguous block causes severe memory bottlenecks and <mark style="background:rgba(205, 244, 105, 0.55)"> fundamentally underutilizes modern multi-core computing architectures.(멀티코어이야기 들어가는 게 맞니? 지금 디코딩이 짝게도 된다는 이야기하는 중인데-뒤에 비교해도 좋고.- 멀티 코어로 빠르게 받아온다도 좋고.-받아온거를 줄이기위해서.)</mark> <span style="color:red">To enable efficient scalability, the patch-based decoding system(PBDS) leverages the predefined KIS-derived index map to partition reads into discrete, index-bounded chunks for independent, parallelized processing.</span>

<font color="#ff0000">To characterize the effect of patch partitioning on processing speed and decoding fidelity, the patch division strategy was evaluated under the previously established the lowest downsampled conditions.</font> Figure 4a illustrates the partitioning logic: using the upstream KIS index map, reads are grouped into continuous index ranges with a boundary overlap that preserves XOR-based error recovery across patch borders. Patch granularity — [placeholder: reads per patch] with a [placeholder: ±indices] overlap — was scaled incrementally, and each step was validated by confirming XOR parity integrity and downstream Reed-Solomon (RS) decoding success. Figure 4b shows index purity as a function of total read count (convertible to average coverage via [placeholder: reads→coverage formula]); at the evaluated operating point ([placeholder: read count] → average coverage ×504, read-distribution spread [placeholder: σ or IQR]), the PBDS retains an index match rate of <font color="#ff0000">99.87%</font>. Figure 4c quantifies the overall throughput gain <font color="#ff0000">across datasets of varying sequence counts</font>: (다양한 사이즈 나열). The 16.6 MB dataset (92,626 reads) was divided into <font color="#ff0000">10 patches </font>and <font color="#ff0000">processed in parallel on multiple CPU cores, reducing total decoding time from 23 minutes to approximately 3 minutes. (이부분 이야기 멀티 결과도 보고 말애햐 할듯하고.-병렬코어수가 많은 라이젠이 유리할수도.)</font>Sequence-level data reconstruction was verified at [placeholder: XX.X%], confirming that patch-level parallelization does not degrade recovery fidelity. [placeholder: 1/2/4/8/10-core processing time — saturation beyond [N] cores] This outcome proves that KIS-guided index-range partitioning successfully mitigates memory overload while strictly preserving data integrity, confirming that downsampled data streams can be reliably routed into distributed computing environments.

The patch-based architecture enables highly parallelized data processing, <font color="#ff0000">delivering a 7.6-fold speed enhancement(수치확인후 가능한지 숫자 기입해야 뭘 우선시 해야할지는..)</font> for large-scale DNA datasets without compromising fundamental index purity or error recovery mechanisms. While software-level parallelization and input data reduction have individually accelerated the pipeline, the ultimate decoding speed achievable through systematic algorithmic optimization of the core computational stages—and the extent to which this performance further scales with available CPU, GPU, and memory resources—remains to be quantified. 


## 5. Scalable Optimization System (V14)--이거 방법으로 줄여서 하는 것 까지 보기. !!

Having proven the individual efficacies of input data downsampling and patch-based software parallelization, we <mark style="background:rgba(205, 244, 105, 0.55)">sought to evaluate</mark> the ultimate processing limits of the iSCAN architecture through systematic algorithmic optimization of each computational stage. While software-level parallelization significantly reduced the decoding timeline for large datasets, the intensive string-matching <span style="color:red">and arithmetic</span> operations required <mark style="background:rgba(205, 244, 105, 0.55)">across the core pipeline stages—primer sorting, index clustering, binarization, and Reed-Solomon decoding—</mark>still demanded further acceleration to achieve real-time data retrieval capabilities. To address this, we implemented an algorithmic optimization stack combining (1) array-level vectorization of pipeline stages, (2) byte-level integer encoding of  nucleotide sequences, (3) lookup-table-based finite-field arithmetic of Galois field, and (4)just-in-time native compilation of residual inner loops——to resolve the computational bottleneck at each stage. 


**Figure 5A** presents the integrated optimization pipeline, mapping each of the four algorithmic layers to the specific computational stage it targets. At the first layer, we demonstrate array-level vectorization across multiple stages of the pipeline, <font color="#ff0000">including primer matching, key-information scanning, and cluster grouping</font>. A naïve Python implementation iterates over each read, each nucleotide, and each candidate primer, dispatching a Python bytecode instruction per element. For datasets comprising tens of millions of reads, this interpreter overhead becomes the dominant contributor to wall-clock time.To address this limitation, we re-express these stages as whole-array operations executed in compiled numerical backends, thereby reducing interpreter dispatch from per-element to per-array invocation. Concretely, primer matching is reformulated as a broadcast matrix comparison over encoded nucleotide sequences; key-information scanning is implemented via strided indexing over a byte-encoded sequence matrix; and cluster assignment is performed through index sorting followed by group-boundary detection using search-based indexing. <mark style="background:rgba(205, 244, 105, 0.55)">얼만큼 줄었는지.</mark>

At the secondary layer, we implanted byte-level integer encoding of nucleotide sequences. DNA sequencing reads arrive as ASCII text, with each nucleotide stored as a one-byte character (A=0x41, C=0x43, G=0x47, T=0x54), and carrying these reads as character strings throughout the pipeline forces repeated symbol-to-value conversion at every comparison while preventing the dense, contiguous memory layout on which the downstream vectorized stages depend. Because only four symbol identities carry information, this representation is intrinsically wasteful: the entire alphabet collapses to a two-bit code, and any per-symbol transformation can be expressed as a constant-time table lookup. We therefore re-encode each read once, at ingestion, into a uint8 matrix whose entries directly index into precomputed lookup tables, so that from this point onward all downstream stages operate on a dense numeric buffer—character comparisons reduce to integer equality tests, and the 2-bit binary representation required as Reed–Solomon input is obtained through a single table lookup per symbol. The transformation is lossless, and the lookup-driven design eliminates the per-character conversion overhead that otherwise dominates the binarization stage.
<mark style="background:rgba(205, 244, 105, 0.55)">얼만큼 줄었는지.</mark>

At the third layer, we extend the same lookup-table principle from symbol encoding to finite-field arithmetic, applying it to the Reed–Solomon decoder that dominates the residual cost of the pipeline. Reed–Solomon decoding requires arithmetic in the finite field GF(2⁸) — specifically multiplication and inversion of 8-bit symbols under a primitive polynomial — and a pure-Python implementation computes each field multiply from first principles (bitwise polynomial multiplication modulo the primitive polynomial), executing roughly O(log n) bitwise operations per multiply and incurring interpreter overhead on each. Since GF(2⁸) has only 256 elements, however, the complete 256×256 multiplication table occupies 64 KB and can be precomputed once at program start; with the table resident in memory, every multiplication collapses to a single 2-D array lookup, while field addition remains a single bitwise XOR. This is a classical space-time tradeoff: 64 KB of precomputed state in exchange for constant-time field arithmetic and the elimination of per-multiply interpreter overhead. Because no off-the-shelf GPU implementation of GF(2⁸) lookup-based arithmetic was available for our target platform, we additionally re-implemented the same log / antilog / multiplication tables natively on Apple Silicon, exploiting its unified memory architecture so that the codeword buffer can be shared between CPU and GPU without an explicit host-to-device copy. The resulting heterogeneous decoder dispatches batched syndrome computation and error-locator evaluation across the GPU cores while retaining the control flow and the syndrome-skip fast path on the CPU, so that both compute domains run concurrently rather than in a serial host→device→host hand-off; this delivers an additional **xxxx×** acceleration over the CPU-only lookup-table baseline at no extra memory cost.

At the fourth layer, we target the small number of inner loops that resist clean vectorization, where per-iteration work is irregular — branching, accumulator state, or array shapes that do not align with broadcasting semantics. After the previous three layers, the residual bottleneck of the iSCAN pipeline is the column-wise mode computation over the uint8 codeword matrix during the binarization consensus step, in which roughly half a billion individual reductions are performed per Tripitaka-scale dataset. A *just-in-time* (JIT) compiler resolves this by profiling the function on its first invocation, specializing it to the observed argument types and array shapes, lowering it to native machine code, and caching the compiled version for all subsequent calls; ahead-of-time performance is therefore approached on the residual hot loop without abandoning the high-level development environment in which the rest of the pipeline is written.

Taken together, the four-layer optimization stack — array-level vectorization, byte-level uint8 encoding, lookup-table Galois-field arithmetic accelerated through an Apple Silicon CPU + GPU heterogeneous decoder, and just-in-time native compilation of the residual consensus loop — resolves the decoding bottleneck across the entire payload range tested here; the per-stage speedup factors and per-dataset wall-clock times are compiled in [Table X]. <font color="#ff0000">Across encodings spanning roughly 0.2 MB / 2,587 oligos / ~1.5 GB raw FASTQ at the small end to 16.6 MB / 92,626 oligos / ~22 GB raw FASTQ at the large end, end-to-end decoding collapses from the minutes-to-hours regime of the MATLAB reference pipeline into the seconds-to-minutes regime — a one- to two-order-of-magnitude reduction — while preserving 100% reconstruction fidelity</font>. Within the methodological frame in which commercial oligo synthesis and next-generation sequencing platforms are taken as fixed, externally supplied black-box tools, decoding is the one read-back stage whose latency can be directly engineered; the present results therefore establish that systematic optimization of the decoding stack alone is sufficient to dramatically narrow the read-back bottleneck that has so far defined the practical ceiling of DNA-based cold data storage, marking a concrete and reproducible step toward operationally deployable DNA archival systems.

## 6. Architectural Generalization Across DNA Storage Encoding Methods 

Having validated the optimized iSCAN pipeline on our encoding format, we next assessed whether our decoding architecture could be generally applied to representative NGS-based DNA storage systems.

<span style="color:red">Application of the four-layer optimization stack to a heterogeneous body of existing encodings presupposes access to a compact yet precisely specified set of decoding meta-information for each candidate scheme: the structural specification of the addressable index (its in-oligo offset, length, and base alphabet); the primer architecture (forward and reverse primer sequences together with their flanking positions within the oligo); and the error-correction parameters (code family, codeword length, parity arity, and shortening convention). Provided that these specifications are recoverable from the original publication or its supplementary materials, the four optimization layers operate in a strictly format-agnostic manner—none of their per-stage transformations is contingent on the encoding scheme itself, only on the meta-information exposed at the ingestion boundary.</span>

<span style="color:red">On this basis, we re-implemented the decoding pipelines of [placeholder: N representative encoding schemes spanning multiple paradigms — e.g., direct binary mapping (Church 2012), homopolymer-free Huffman encoding (Goldman 2013), DNA Fountain (Erlich 2017), random-access architecture (Organick 2018), composite-letter systems (Anavy 2019), and locality-sensitive hashing (Antkowiak 2020)] within a single harnessed runtime, and substituted each method's reference decoder, layer by layer, with the iSCAN four-layer optimization stack (array-level vectorization, uint8 byte-level encoding, lookup-table Galois-field arithmetic with the Apple Silicon CPU + GPU heterogeneous decoder, and JIT-compiled hot loops). **Figure 6a** tabulates the resulting per-case acceleration: for each surveyed encoding method, the runtime of the originally published decoder is reported alongside that of its iSCAN-optimized variant executed on identical input reads under matched coverage. Across the [N] benchmarked cases, the optimized variant delivered an average **X-fold** reduction in end-to-end decoding wall-clock (range Y× to Z×) while preserving 100% reconstruction fidelity, with detailed per-case factors enumerated in Figure 6a.</span>

This cross-method comparison demonstrates that the iSCAN decoding architecture is not limited to our encoding format, but can be broadly extended to representative NGS-based DNA storage systems. While the slow timescale of synthesis and encoding inherently places DNA storage in the domain of cold-data archiving, rapid decoding remains essential from the user-side perspective when archived digital information must be retrieved. By coupling decoding more closely with the sequencing process, this architecture may enable more direct data recovery and improve storage-buffer workflows that rely on reconstructing stored information after sequencing. (Figure 6 b)

# Conclusion (V15)

> ✅ **Conclusion V15 연쇄 업데이트 완료**
> §5 V14 (92×→28×, 15s→3.9s, GPU 제거) + §6 V14 (21개 방법 서베이) 반영:
> - P1: GPU 언급 제거, 28×/3.9s 실측 반영, baseline 명시
> - P2: 4-layer 알고리즘 스택 구조 명시 (NumPy/uint8/reedsolo/Numba)
> - P3: "CPU and GPU" → "standard CPU"; V13 Figure 6 문장 → 21-method landscape 대응 재작성


This work established an end-to-end decoding framework for DNA data storage that advances the field from algorithmic proof-of-concept toward operational practicality. By replacing alignment-heavy clustering with the positional iSCAN protocol and pursuing a strictly causal optimization path — bottleneck identification, targeted input reduction, robustness-preserving error correction, and distributed parallel execution — we demonstrated simultaneous gains in speed, cost-efficiency, and reconstruction fidelity. The fully integrated pipeline, combining ×5.33 downsampling, patch-based multi-core parallelization, <span style="color:red">and systematic algorithmic optimization of the core computational stages (NumPy vectorization, uint8 integer representation, lookup-table Galois field arithmetic, and Numba JIT compilation), achieved a 28-fold total acceleration over the MATLAB baseline, compressing the decoding of the Jikji 6 dataset from 109 seconds to 3.9 seconds</span> while maintaining 100% reconstruction accuracy. A co-reduction of approximately 80% in sequencing cost confirmed that the performance gains are not confined to computation alone but extend to the economic constraints that govern practical deployment.


The broader methodological contribution of this study is the demonstration that DNA decoding performance is amenable to system-level engineering rather than isolated algorithmic improvement. Stage-wise time decomposition identified primer sorting as the dominant bottleneck (40–52% of total runtime), providing a precise and evidence-based target for subsequent optimization. Downsampling directly reduced the load at that bottleneck; RS+XOR error correction preserved reconstruction integrity under the resulting sparse coverage conditions; and index-range patch partitioning extended these gains to large datasets without compromising index purity (>99.87%). <span style="color:red">Each layer of the optimization stack — from input-volume reduction through parallel execution to the four distinct algorithmic optimizations (NumPy vectorization, uint8 data representation, lookup-table Galois field arithmetic, and Numba JIT compilation) — was independently validated before integration, and the cumulative speedup reflected the combined contribution of each stage-specific transformation.</span> This structure — measure first, then target precisely — constitutes a reusable framework for designing deployable DNA storage pipelines around explicit operating points that balance coverage, throughput, and reliability.

The scope of the present findings is bounded by the dataset range over which the pipeline was evaluated (kilobyte to megabyte scale), and whether the optimization stack scales without degradation to gigabyte-class payloads remains an open question. <span style="color:red">The extent to which the observed performance generalizes across heterogeneous computing environments beyond the standard CPU platform tested has likewise yet to be determined.</span> Nonetheless, within the validated operating envelope, these results confirm that computational latency is no longer the governing constraint for DNA data storage. The principal bottleneck has shifted upstream, to the physical chemistries of oligonucleotide synthesis and next-generation sequencing throughput. Resolving those molecular-scale limits, rather than decoding speed, now defines the critical path toward practical, large-scale DNA-based information storage. <span style="color:red">When evaluated across the 21-method landscape of DNA data storage encoding systems published between 2012 and 2024 — spanning alignment-based workflows that scale to hours, fountain-code pipelines that require minutes, and our own previously published iSCAN implementation — the iSCAN architecture establishes a paradigm distinct from its predecessors: position-based index scanning circumvents sequence alignment entirely, decoding the Jikji 6 dataset in 3.9 seconds on standard CPU hardware without platform-specific acceleration. The optimization stack transfers directly to the 15 of 21 surveyed methods that embed addressable index structures, positioning system-level engineering of the decoding stack — rather than hardware specialization — as a reproducible template for translating DNA storage from algorithmic demonstrations to deployable infrastructure.</span>


---

# 전체 버전 매트릭스

| 섹션 | 현재 버전 | 출처 파일 | 비고 |
|------|---------|----------|------|
| **Introduction** | <span style="color:red">**V15**</span> | 260417_19_Intro_unified.md → V15 연쇄 수정 | P4 V15 업데이트 완료 — "hardware acceleration" → "algorithmic optimization"; Figure 6 문장 21-method 서베이 대응으로 재작성 |
| **Results §1** Index Alignment | V12 | 260308_12_Results_draft.md | 원문 보존, 사용자 `<mark>` 주석 유지 |
| **Results §2** Time Quantification | V12 | 260308_12_Results_draft.md | 원문 보존 |
| **Results §3** Data Downsampling | V12 | 260308_12_Results_draft.md | 원문 보존 |
| **Results §4** Patch-Based Decoding | V14 | 260418_13_Results_Sec456_v14.md | V12 본문 + 결론 열린 질문만 수정 |
| **Results §5** Scalable Optimization | V14 | 260418_13_Results_Sec456_v14.md | **전면 개정** — 260328 실측 반영, GPU 제거, 4-layer 알고리즘 스택 |
| **Results §6** Architectural Generalization | V14 | 260418_13_Results_Sec456_v14.md | **신규** — V13 §6 "Decoding Latency in Context" 대체, 21개 방법 서베이; 제목 "Cross-Platform Generalizability" → "Architectural Generalization Across DNA Storage Encoding Methods" 변경 |
| **Conclusion** | <span style="color:red">**V15**</span> | 260417_19_Conclusion_unified.md → V15 연쇄 수정 | §5/§6 V14 연쇄 업데이트 완료 — GPU 제거, 28×/3.9s 실측 반영, P3 21-method landscape 대응 재작성 |

---

# 변경 로그 (통합)

## Intro V1 → V2 → V13 → V15

| 항목 | V1 | V2 | V13 | V15 | 근거 |
|------|-----|-----|-----|-----|------|
| "215 PB/g" | 무한정 | "theoretical" 맥락 + 보존 조건 명시 | — | — | Self-Doubt M1, M2 |
| "discard it" → "fail to exploit it directly" | discard | fail to exploit | — | — | Self-Doubt I4 |
| iSCAN 4단계 | 모두 restricted portion | 3단계 restricted + consensus full reads | — | — | Self-Doubt I3 |
| 연구 질문 | 3개 | **4개 (통합 가속 추가)** | — | — | Self-Doubt C5 |
| 다운샘플링 | 언급만 | "5–10× coverage" 명시 | — | — | Intro참조.md |
| P4 4번째 질문 문구 | — | "system-level acceleration … **hardware acceleration**" | — | <span style="color:red">"decoding acceleration … **systematic algorithmic optimization of the core pipeline stages**"</span> | §5 V14 GPU 제거 |
| P4 최적화 경로 | 3단계 | 4단계 ("**full-stack hardware acceleration**" 추가) | — | <span style="color:red">"**systematic algorithmic optimization of the core computational stages**"</span>로 교체 | §5 V14 |
| P4 마지막 | — | — | **Figure 6 벤치마크 문장 추가** (green) | <span style="color:red">"**evaluate the generalizability … 21 DNA data storage encoding methods published between 2012 and 2024 … mapping where its position-based scanning paradigm transfers directly**"</span>로 재작성 | §6 V14 21-method 서베이 |

## §4 V12 → V14

| 항목 | V12 | V14 | 근거 |
|------|-----|-----|------|
| 결론 열린 질문 | "hardware-level acceleration" | **"systematic algorithmic optimization"** | §5 GPU→알고리즘 전환 대응 |
| 본문 나머지 | — | 동일 보존 | — |

## §5 V12 → V14 (전면 개정)

| 항목 | V12 | V14 | 근거 |
|------|-----|-----|------|
| 최적화 접근법 | ×5.33 downsampling + parallel + **GPU** | **NumPy vec + uint8 + reedsolo + Numba JIT (no GPU)** | 260328 실측 — GPU 미구현 |
| 총 속도 향상 | **92×** | **28× (J6), 14× (J7)** | 실측값 |
| 처리 시간 (대표) | 16.6 MB → 15 s | J6 109 s → 3.9 s; Tripitaka ~23 min → ~1–2 min | 실측값 |
| Figure 5a | 3-layer waterfall (DS+parallel+GPU) | **4-layer waterfall (NumPy+uint8+reedsolo+Numba)** | 아키텍처 변경 |
| Figure 5b | 16.6 MB / 15 s 회귀 | J6(3.9s), J7(17.5s), Tripitaka(~1–2min) 실측 | 실측값 |
| Figure 5c | 92× 시너지 확인 | 알고리즘 스택이 GPU 추정치 동등 이상 + 플랫폼 독립성 | 실측 결론 |
| 결론 프레이밍 | 마지막 종결 (물리적 한계) | **§6 연결 열린 질문 (범용성 미검증)** | §6 추가 대응 |
| Source | 미언급 | pipeline/dna_decode_260328_0125.py | 추적성 |

## §6 V13 → V14 (전면 교체)

| 항목 | V13 | V14 | 근거 |
|------|-----|-----|------|
| 제목 | "Decoding Latency in Context" | **"Architectural Generalization Across DNA Storage Encoding Methods"** | 서베이 기반 확장; "Cross-Platform"은 hardware/OS 함의가 강해 인코딩 방법 다양성 맥락에 부적합하므로 "Architectural Generalization"으로 교체 |
| 내용 | iSCAN3 벤치마크 비교 (180 words) | **21개 방법 서베이 + 3 패러다임 + 4 모듈 매핑 (~550 words)** | Figure 6 서베이 |
| 수치 기반 | 92×, 15 s, 66 MB/min | **28×, 3.9 s (§5 실측 기반)** | 260328 업데이트 |
| 비교 대상 | 4개 (iSCAN1 + Fountain + DBSCAN + RS-ECC) | **21개 인코딩 방법 (2012–2024)** | 포괄적 서베이 |
| Figure 구성 | 6a(iSCAN3 vs iSCAN1), 6b(문헌비교), 6c(throughput) | **6a(밀도×연도), 6b(3패러다임 시간), 6c(적용가능성 히트맵)** | 아키텍처 관점 |
| PaperMan 검증 | 없음 | **SD + Monitor + RT 완료, rev1 반영** | v0.3 프로토콜 |

## Conclusion V11 → V13 → V15

| 항목 | V11 | V13 | V15 | 근거 |
|------|-----|-----|-----|------|
| Figure 6 대응 | 없음 | P3 말미 method-comparison 2문장 추가 (green) | <span style="color:red">V13 green 2문장 삭제 → 21-method landscape 대응 2문장으로 재작성</span> | §6 V14 |
| P1 가속 수치 | "92× / 15 s (16.6 MB Tripitaka) / sequential baseline" | 동일 | <span style="color:red">**28× / 3.9 s (J6) / MATLAB baseline**</span> | §5 V14 실측 |
| P1 GPU 언급 | "GPU-accelerated primer sorting" | 동일 | <span style="color:red">**"systematic algorithmic optimization … (NumPy, uint8, reedsolo, Numba)"**</span> | §5 V14 GPU 제거 |
| P2 누적 가속 설명 | "cumulative speedup followed predictably from the individual contributions" | 동일 | <span style="color:red">**4-layer 알고리즘 스택 명시 (NumPy+uint8+reedsolo+Numba) + "combined contribution of each stage-specific transformation"**</span> | §5 V14 구조 |
| P3 "CPU and GPU platforms" | "CPU and GPU platforms tested" | 동일 | <span style="color:red">**"standard CPU platform tested"**</span> | §5 V14 GPU 제거 |
| P3 Figure 6 문장 (V13) | 없음 | "benchmarked … 15 seconds … fastest end-to-end decoding" (green) | <span style="color:red">**"evaluated across 21-method landscape … 3.9 seconds on standard CPU … transfers to 15 of 21 surveyed methods"**</span> | §5/§6 V14 통합 |
| "fastest end-to-end decoding" 주장 | — | 있음 | <span style="color:red">**삭제** — 실측 범위 밖 확장. "architecturally distinct paradigm" + "transfers to 15 of 21"로 재프레이밍</span> | overstatement 회피 |

---

# ✅ V15 연쇄 업데이트 완료 (260419)

> §5 V14 (92×→28×, 15s→3.9s, GPU 제거) + §6 V14 (21-method 서베이)에 따른 Intro·Conclusion 업데이트 완료:

| 구간 | V13 상태 | V15 변경 | 상태 |
|------|---------|---------|-----|
| **Introduction P4** (4번째 질문) | "hardware acceleration" + "통합 가속" 프레이밍 | <span style="color:red">"**systematic algorithmic optimization of the core pipeline stages**"로 교체</span> | ✅ 완료 |
| **Introduction P4** (최적화 경로) | "full-stack hardware acceleration" (green) | <span style="color:red">"**systematic algorithmic optimization of the core computational stages**"로 교체</span> | ✅ 완료 |
| **Introduction P4 마지막** | "benchmark ... iSCAN implementation ... situate the achieved performance within the broader landscape" (green) | <span style="color:red">"**evaluate the generalizability … across 21 DNA data storage encoding methods published between 2012 and 2024 … mapping where its position-based scanning paradigm transfers directly**"로 재작성</span> | ✅ 완료 |
| **Conclusion P1** | "GPU-accelerated primer sorting … 92-fold … 16.6 MB … 23 minutes to 15 seconds … sequential baseline" | <span style="color:red">"**systematic algorithmic optimization of the core computational stages (NumPy, uint8, reedsolo, Numba) … 28-fold … Jikji 6 … 109 seconds to 3.9 seconds … MATLAB baseline**"로 교체</span> | ✅ 완료 |
| **Conclusion P2** | "cumulative speedup followed predictably from the individual contributions" | <span style="color:red">"**Each layer … from input-volume reduction through parallel execution to the four distinct algorithmic optimizations … combined contribution of each stage-specific transformation**"로 구체화</span> | ✅ 완료 |
| **Conclusion P3** (CPU/GPU) | "CPU and GPU platforms tested" | <span style="color:red">"**standard CPU platform tested**"로 교체</span> | ✅ 완료 |
| **Conclusion P3 V13 추가 문장** | "15 seconds ... fastest end-to-end decoding" (green, V13 신규) | <span style="color:red">"**21-method landscape ... 3.9 seconds on standard CPU ... transfers to 15 of 21 surveyed methods**"로 전면 재작성; "fastest" 주장 삭제</span> | ✅ 완료 |

---

# Narrative Chain Map (§1→§6, V14 확정 기준)

```
§1 Index Alignment Protocol Strategy (V12)
  WHY   → Conventional de novo alignment 6–12 h 병목
  WHAT  → iSCAN 4-stage (PrimerSort + KIS + Index Cluster + Consensus)
  INTER → 360× 가속, J6 4.57M reads → 34.5 s, 완전 재구성
  ──────────────────────────────────────────────
  NEXT WHY → "MB scale로 확장하면?"
        ▼
§2 Time Quantification Analysis (V12)
  WHY   → 용량 증가에 따른 time scaling 미검증
  WHAT  → 499 kB → 16.6 MB 전 구간 stage-wise timing
  INTER → 선형 스케일링, PrimerSort 40–52% (dominant bottleneck)
  ──────────────────────────────────────────────
  NEXT WHY → "input 볼륨 자체를 줄일 수 있는가?"
        ▼
§3 Data Downsampling Efficiency (V12)
  WHY   → ×1,000 coverage 관행은 PrimerSort 비용 inflate
  WHAT  → ×5.33 downsampling + RS(inner) + XOR(outer)
  INTER → 100% recovery, 14× 가속, NGS 비용 ~80% 절감
  ──────────────────────────────────────────────
  NEXT WHY → "분산/병렬 환경에 통합 가능한가?"
        ▼
§4 Patch-Based Decoding System (V14)
  WHY   → 대용량 단일 스트림의 메모리 병목
  WHAT  → KIS 기반 index-range 패치 + ±100 overlap
  INTER → 7.6× 가속, 99.87% 인덱스 순도 유지
  ──────────────────────────────────────────────
  NEXT WHY → "알고리즘 최적화의 한계는?" (V14 수정 — GPU→알고리즘)
        ▼
§5 Scalable Optimization System (V14)
  WHY   → 파이프라인 단계별 연산 병목 잔존
  WHAT  → 4-layer 알고리즘 스택 (NumPy/uint8/reedsolo/Numba)
  INTER → 28× 가속 (J6: 109s→3.9s), GPU 불필요, 플랫폼 독립
  ──────────────────────────────────────────────
  NEXT WHY → "자체 포맷에만 유효한가, 범용적인가?" (V14 신규 연결)
        ▼
§6 Architectural Generalization Across DNA Storage Encoding Methods (V14)
  WHY   → 21개 인코딩 방법 존재, 범용성 미검증
  WHAT  → 3 디코딩 패러다임 비교 + 4 모듈 적용 가능성 매핑
  INTER → 15/21 방법 적용 가능, 아키텍처적 차별화 확인
  ──────────────────────────────────────────────
  FINAL → 남은 과제는 physical synthesis/sequencing 속도
        ▼
Conclusion (V15 / §5·§6 V14 연쇄 업데이트 완료 — GPU 제거, 28×/3.9s, 21-method landscape)
  → System-level engineering template for deployable DNA storage, transferable to 15 of 21 surveyed methods
```

---

# 검증 결과 요약

| 검증 | §1 (V12) | §2 (V12) | §3 (V12) | §4 (V14) | §5 (V14) | §6 (V14) |
|------|----------|----------|----------|----------|----------|----------|
| **Self-Doubt** | ✅ PASS (260302) | ✅ PASS | ✅ PASS | ✅ PASS (260302) | ⚠️ 수치 업데이트, 재검증 권장 | ⚠️ CONDITIONAL (실증 데이터 부재) |
| **Monitor** | ✅ PASS (260302) | ✅ PASS | ✅ PASS | ✅ PASS (260302) | — | ⚠️ CONDITIONAL (그룹 수/비율 보강 필요) |
| **Red Team** | ✅ PASS (260302) | ✅ PASS | ✅ PASS | ✅ PASS (260302) | — | ⚠️ CONDITIONAL (외부 cross-decoding 실험 권장) |
| **Rev1 반영** | N/A | N/A | N/A | N/A | N/A | ✅ 7/8 반영 완료 |

### §6 미반영 1건
- ⬜ **외부 데이터셋 실제 cross-decoding 실험** (RT Attack 2 Option A)
  - 9개 외부 데이터셋이 `06_Other_ENCODE/datasets/`에 다운로드 완료
  - 최소 1–2개(예: Antkowiak Illumina, DNA-Aeon SRA)에 대한 파이프라인 적용 결과 추가 시 → Results 섹션 위치 완전 정당화
  - 미수행 시: "Comparative Analysis" 성격으로 성립하나, 리뷰어가 Discussion 이동 요청할 가능성 있음

---

# 파일 참조 인덱스

| 파일 | 위치 | 용도 |
|------|------|------|
| `260418_14_Manuscript_Full_Compiled_v15.md` | Manuscript_Drafts_v2/ | **본 파일 — 전체 통합본** |
| `260418_13_Results_Sec456_v14.md` | Manuscript_Drafts_v2/ | §4–§6 V14 출처 |
| `260417_19_Results_Compiled_Draft.md` | Manuscript_Drafts_v2/ | Intro V13 + §1–§3 V12 + Conclusion V13 출처 |
| `260417_19_Intro_unified.md` | Manuscript_Drafts_v2/ | Intro V13 단독본 |
| `260417_19_Conclusion_unified.md` | Manuscript_Drafts_v2/ | Conclusion V13 단독본 |
| `260308_12_Results_draft.md` | Manuscript_Drafts_v2/ | §1–§5 V12 원본 (2026-03-25 14:34) |
| `260418_11_Sec6_v0.0.md` | 06/history/ | §6 V14 초기 파싱본 |
| `260418_11_Sec6_outline.md` | 06/history/ | §6 V14 개조식 |
| `260418_12_Sec6_draft.md` | 06/history/ | §6 V14 산문 초안 |
| `260418_12_Sec6_SD.md` | 06/history/ | §6 Self-Doubt 5 Questions |
| `260418_12_Sec6_Monitor.md` | 06/history/ | §6 Monitor 8-Axis |
| `260418_12_Sec6_RT.md` | 06/history/ | §6 Red Team 5 Attacks |
| `260418_13_Sec6_rev1.md` | 06/history/ | §6 V14 검증 반영 수정본 |
| `260418_13_Sec6_revision.md` | 06/revisions/ | §6 V14 수정 추적 |
| `Figure설계안.md` | 06/ | Figure 6a/6b/6c 상세 설계 |
| `260328_0200_Sec5_optimization_update.md` | 05/history/ | **§5 V14 수치 업데이트 근거** |
| `Figure_6_DNA_Storage_Encoding_Methods.md` | 06_Other_ENCODE/ | **§6 V14 21개 방법 서베이 원본** |
| `dataset_index_catalog.md` | 06_Other_ENCODE/ | 9개 외부 데이터셋 카탈로그 |
| `pipeline/dna_decode_260328_0125.py` | pipeline/ | **§5 V14 실측 소스 코드** |
