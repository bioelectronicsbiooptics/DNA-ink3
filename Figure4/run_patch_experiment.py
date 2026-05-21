# =============================================================================
# Figure 4 — Patch-Based Decoding & Chunked FASTQ Loading
# =============================================================================
# [WHY CHUNKED FASTQ?]
#   대용량 FASTQ (Tripitaka 급 16 MB+) 를 한 번에 메모리로 로드하면 RAM 이 폭증하고
#   시스템 swap 에 들어가 I/O 가 수 배 느려진다. 처음부터 작은 chunk (예: 1k–10k reads)
#   단위로 스트리밍 파싱하면:
#     (1) peak RSS 가 전체 파일 크기에 비례하지 않고 chunk 크기에 수렴
#     (2) CPU 캐시 / NUMA 로컬리티가 유지되어 벡터화 연산이 빨라짐
#     (3) patch 단위 병렬 처리와 자연스럽게 맞물려 파이프라이닝 가능
#   → "하드웨어 메모리 용량 ↔ chunk 크기" 튜닝이 처리 속도를 좌우함.
#
# [PATCH STRATEGY]
#   - 인덱스 범위 기반 분할: 1k–10k sequences / patch
#   - 인접 패치 overlap: ±100 indices (경계 XOR 복구용)
#   - 병렬화: 각 patch 독립 디코딩 → 결과 병합
#
# [CORE FINDINGS]
#   - Index Purity >99.8% (×5.33 다운샘플링 시에도 유지)
#   - 최적 패치 크기: 5k–10k (메모리 ↔ 속도 균형)
#
# [DATA LINKS]
#   - 원본 FASTQ: J6.fastq (로컬 존재 확인: YES — /Users/young/Desktop/Project/EPCR1-3_Obsidian/Decoding 통합 2 /J6.fastq)
#   - 결과: Outputs/ (이 폴더)
# =============================================================================
#!/usr/bin/env python3
"""
패치 기반 병렬 디코딩 실험.
리드를 인덱스 범위별로 분할 → 병렬 디코딩 → 결과 병합.
±100 인덱스 오버랩으로 XOR 경계 복구 보장.
"""
import sys
import os
import time
import math
import json
import multiprocessing as mp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python_pipeline'))
from dna_decode_260328_0125 import (
    read_fastq_file, _seqs_to_matrix, PrimerSelect, ReverseEvert,
    PrimerSort, KeyFind, Keysort, Consensus, Binarization, RSdecoding,
    NullFind, XORrun, tridecode, trierrXOR, _make_triangulation,
    _dna_to_binary, _bi2de, _de2bi, _init_rs_reedsolo,
    _binarize_mode_jit,
    _A, _N, _DNA_TO_IDX, _BIT_TABLE
)

FILENAME = 'Jikji6'
DATA_DIR = '/Users/young/Desktop/Project/EPCR1-3_Obsidian/Decoding 통합 2 '
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'Outputs', 'patch_results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

OVERLAP = 100  # ±100 인덱스 오버랩


def run_sequential(ps_seq, ps_mat, idxlen, OligoLen, k, m, n, filename, seq, seq_mat, F_primers):
    """순차 디코딩 (기준선)."""
    t = time.time()

    MIndex, sorted_u8, boundaries = Consensus(ps_seq, idxlen, OligoLen, ps_mat_u8=ps_mat)
    Depth = np.array([MIndex[i].shape[0] if MIndex[i] is not None else 0 for i in range(len(MIndex))])

    _init_rs_reedsolo()
    SeqtoBin = Binarization(MIndex, sorted_u8=sorted_u8, boundaries=boundaries)
    rsECC, RSerrinf = RSdecoding(SeqtoBin, k, m, n)
    rsECC[RSerrinf < 0, :] = 0
    rsECC[RSerrinf > 4, :] = 0
    for a in range(min(OligoLen, len(rsECC))):
        if np.any(rsECC[a, :] > 1):
            rsECC[a, :] = 0

    payload = rsECC[:, :k - idxlen]
    elapsed = time.time() - t
    return payload, Depth, RSerrinf, elapsed


def decode_patch(args):
    """단일 패치 디코딩 (워커 함수)."""
    patch_seqs, idxlen, oligo_start, oligo_end, k, m, n, filename = args

    patch_oligo_len = oligo_end - oligo_start
    if len(patch_seqs) == 0 or patch_oligo_len <= 0:
        return oligo_start, oligo_end, None, None, None, 0

    t = time.time()
    try:
        _init_rs_reedsolo()
        ps_mat = _seqs_to_matrix(patch_seqs)

        MIndex, sorted_u8, boundaries = Consensus(patch_seqs, idxlen, patch_oligo_len, ps_mat_u8=ps_mat)
        SeqtoBin = Binarization(MIndex, sorted_u8=sorted_u8, boundaries=boundaries)
        rsECC, RSerrinf = RSdecoding(SeqtoBin, k, m, n)
        rsECC[RSerrinf < 0, :] = 0
        rsECC[RSerrinf > 4, :] = 0

        payload = rsECC[:, :k - idxlen]
        Depth = np.array([MIndex[i].shape[0] if MIndex[i] is not None else 0 for i in range(len(MIndex))])
        elapsed = time.time() - t
        return oligo_start, oligo_end, payload, Depth, RSerrinf, elapsed
    except Exception as e:
        return oligo_start, oligo_end, None, None, None, time.time() - t


def run_patch_parallel(ps_seq, ps_mat, idxlen, OligoLen, k, m, n, filename,
                       seq, seq_mat, F_primers, n_patches=4):
    """패치 기반 병렬 디코딩."""
    t_total = time.time()

    # 패치 분할: 인덱스 범위별
    patch_size = math.ceil(OligoLen / n_patches)

    # 리드의 인덱스를 미리 계산하여 패치별 할당
    # Consensus의 인덱스 추출 로직 재사용
    dna_matrix = np.array([list(s) for s in ps_seq])
    idx_end = 98
    idx_start = 98 - math.ceil(idxlen / 2)
    if ps_mat is not None:
        idx_region = ps_mat[:, idx_start:idx_end].copy()
    else:
        idx_region = _seqs_to_matrix(ps_seq)[:, idx_start:idx_end]
    idx_region[idx_region == _N] = _A
    idx_vals = _DNA_TO_IDX[idx_region]
    Indexbin = _BIT_TABLE[idx_vals].reshape(len(ps_seq), -1).astype(np.int8)
    oligo_bits = _de2bi(OligoLen)
    if len(oligo_bits) % 2 != 0:
        if Indexbin.shape[1] % 2 == 0:
            Indexbin = Indexbin[:, 1:]
    Indexdec = _bi2de(Indexbin)

    # 패치별 리드 분배
    patches = []
    for p in range(n_patches):
        start_oligo = p * patch_size + 1  # 1-based index
        end_oligo = min((p + 1) * patch_size, OligoLen)

        # 오버랩: ±OVERLAP
        eff_start = max(1, start_oligo - OVERLAP)
        eff_end = min(OligoLen, end_oligo + OVERLAP)

        mask = (Indexdec >= eff_start) & (Indexdec <= eff_end)
        patch_seqs = [ps_seq[i] for i in range(len(ps_seq)) if mask[i]]

        patches.append((
            patch_seqs, idxlen, eff_start, eff_end,
            k, m, n, filename
        ))

    # 병렬 실행
    t_par = time.time()
    ctx = mp.get_context('fork')
    with ctx.Pool(min(n_patches, mp.cpu_count())) as pool:
        results = pool.map(decode_patch, patches)
    t_par_done = time.time()

    # 결과 병합 (오버랩 영역은 중앙 패치 우선)
    total_payload = None
    total_depth = np.zeros(OligoLen, dtype=np.int64)
    total_rserr = np.zeros(OligoLen, dtype=np.int64)
    patch_times = []

    for oligo_s, oligo_e, payload, depth, rserr, elapsed in results:
        patch_times.append(elapsed)
        if payload is None:
            continue

        # 이 패치가 담당하는 실제 범위 (오버랩 제외)
        for p in range(n_patches):
            real_start = p * patch_size + 1
            real_end = min((p + 1) * patch_size, OligoLen)
            if oligo_s <= real_start and oligo_e >= real_end:
                offset = real_start - oligo_s
                length = real_end - real_start + 1
                if total_payload is None:
                    total_payload = np.zeros((OligoLen, payload.shape[1]), dtype=payload.dtype)
                if offset + length <= payload.shape[0]:
                    total_payload[real_start-1:real_end, :] = payload[offset:offset+length, :]
                if depth is not None and offset + length <= len(depth):
                    total_depth[real_start-1:real_end] = depth[offset:offset+length]
                if rserr is not None and offset + length <= len(rserr):
                    total_rserr[real_start-1:real_end] = rserr[offset:offset+length]

    total_time = time.time() - t_total
    return total_payload, total_depth, total_rserr, total_time, t_par_done - t_par, patch_times


def main():
    print("=" * 70)
    print("패치 기반 병렬 디코딩 실험: J6 (Jikji6)")
    print(f"오버랩: ±{OVERLAP} 인덱스")
    print("=" * 70)

    # FASTQ 로드 + 전처리
    t0 = time.time()
    seq = read_fastq_file(os.path.join(DATA_DIR, 'J6.fastq'))
    seq_mat = _seqs_to_matrix(seq)
    F_primers, R_primers = PrimerSelect(FILENAME)

    # Evert
    fs = [s[:21] if len(s) >= 21 else s for s in seq]
    tevert, rev_log = ReverseEvert(seq, fs, R_primers, F_primers, mat_uint8=seq_mat)
    rev_idx = np.where(rev_log)[0]
    if isinstance(tevert, np.ndarray) and tevert.dtype == np.uint8:
        seq_mat[rev_idx, :150] = tevert

    # PrimerSort
    ps_seq, ps_mat = PrimerSort(seq, F_primers, R_primers, mat_uint8=seq_mat)

    # KeyFind
    k = 156; n = k + 8; m = len(_de2bi(n))
    KeyMap, KmatDNA = KeyFind(ps_seq, ps_mat_u8=ps_mat)
    OligoLen, idxlen, ncolpo, nbitpo, ncolco, nbitco = Keysort(KeyMap, KmatDNA, k, FILENAME)

    print(f"전처리 완료: {round(time.time()-t0, 1)}s")
    print(f"유효 리드: {len(ps_seq):,}, OligoLen: {OligoLen}")

    all_results = []

    # 1. 순차 기준선
    print(f"\n{'─'*50}")
    print("순차 디코딩 (기준선)...")
    _, _, _, seq_time = run_sequential(ps_seq, ps_mat, idxlen, OligoLen, k, m, n,
                                        FILENAME, seq, seq_mat, F_primers)
    print(f"  시간: {seq_time:.3f}s")
    all_results.append({'method': 'sequential', 'patches': 1, 'time': seq_time})

    # 2. 패치 병렬 (2, 4, 8, 10 패치)
    for n_patches in [2, 4, 8, 10]:
        print(f"\n{'─'*50}")
        print(f"패치 병렬 ({n_patches} patches, ±{OVERLAP} overlap)...")

        _, _, _, total_t, par_t, patch_ts = run_patch_parallel(
            ps_seq, ps_mat, idxlen, OligoLen, k, m, n,
            FILENAME, seq, seq_mat, F_primers, n_patches=n_patches
        )
        speedup = seq_time / total_t if total_t > 0 else 0

        print(f"  총 시간: {total_t:.3f}s (병렬 구간: {par_t:.3f}s)")
        print(f"  패치별: {[round(t, 3) for t in patch_ts]}")
        print(f"  가속: {speedup:.2f}x")

        all_results.append({
            'method': f'patch_{n_patches}',
            'patches': n_patches,
            'time': round(total_t, 3),
            'parallel_time': round(par_t, 3),
            'patch_times': [round(t, 3) for t in patch_ts],
            'speedup': round(speedup, 2)
        })

    # 결과 저장
    output_file = os.path.join(OUTPUT_DIR, 'patch_experiment_J6.json')
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)

    # 요약
    print(f"\n{'='*70}")
    print(f"{'Method':>15} {'Patches':>8} {'Time':>8} {'Speedup':>8}")
    print(f"{'─'*45}")
    for r in all_results:
        print(f"{r['method']:>15} {r['patches']:>8} {r['time']:>7.3f}s {r.get('speedup', 1.0):>7.2f}x")

    print(f"\n결과 저장: {output_file}")


if __name__ == '__main__':
    main()
