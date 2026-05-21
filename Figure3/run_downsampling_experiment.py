#!/usr/bin/env python3
# =============================================================================
# Figure 3 — Downsampling Efficiency (다운샘플링 효율성 결과)
# =============================================================================
# [OBJECTIVE] 커버리지 감소가 복호화 정확도/속도에 미치는 영향 측정.
# [CORE FINDING] ×504 (×5.33 다운샘플링) 에서도 100% 정확도 유지, NGS 비용 1/5 감소.
# [MEASURED] J4 기준 34.5s → 2.47s (다운샘플링 적용).
# [DATA LINKS]
#   - 원본 FASTQ: 4_S2_L001_R1_001.fastq + 4_S2_L001_R2_001.fastq (로컬 부재 — 루트에 없음)
#   - 기존 결과: Downsampling_Results/Jikji4/ (존재 확인됨)
#   - Venus 데이터셋: Venus1_R1/R2, Venus2_R1/R2 (로컬 부재 가능)
# [NOTES] Venus 는 불균등 커버리지 → dropout + XOR 복구 적용.
#         PCR random access 로 0.52× 초저커버리지에서도 복구 성공.
# =============================================================================
"""
다운샘플링 실험: J6 데이터에서 coverage별 디코딩 성능 측정.
원본 파이프라인(dna_decode_260328_0125.py)을 사용.
"""
import sys
import os
import time
import random
import json
import numpy as np

# 파이프라인 모듈 경로
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python_pipeline'))
from dna_decode_260328_0125 import (
    read_fastq_file, _seqs_to_matrix, PrimerSelect, ReverseEvert,
    PrimerSort, KeyFind, Keysort, Consensus, Binarization, RSdecoding,
    NullFind, XORrun, nestedXOR, reRSdec, tridecode, trierrXOR,
    _make_triangulation, _dna_to_binary, _bi2de, _de2bi,
    _init_rs_reedsolo, _binarize_mode_jit,
    _A, _N, _DNA_TO_IDX, _BIT_TABLE, _COMP_TABLE,
    visualize_mesh
)

# J6 설정
FILENAME = 'Jikji6'
DATA_DIR = '/Users/young/Desktop/Project/EPCR1-3_Obsidian/Decoding 통합 2 '
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'Outputs', 'downsampling_results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 다운샘플 비율 (원본 depth ~1421x 기준)
DOWNSAMPLE_TARGETS = [1421, 500, 200, 100, 50, 20, 10]  # coverage 목표


def downsample_reads(seq_list, target_coverage, original_coverage):
    """리드를 랜덤 다운샘플링."""
    ratio = target_coverage / original_coverage
    if ratio >= 1.0:
        return seq_list
    n_target = max(1, int(len(seq_list) * ratio))
    random.seed(42)  # 재현성
    return random.sample(seq_list, n_target)


def run_pipeline(seq, filename):
    """파이프라인 실행 → (성공여부, 시간dict, 결과dict)"""
    results = {}
    timings = {}

    try:
        # 행렬 구성
        seq_mat = _seqs_to_matrix(seq)

        # Evert
        t = time.time()
        F_primers, R_primers = PrimerSelect(filename)
        fs = [s[:21] if len(s) >= 21 else s for s in seq]
        tevert, rev_log = ReverseEvert(seq, fs, R_primers, F_primers, mat_uint8=seq_mat)
        rev_idx = np.where(rev_log)[0]
        if isinstance(tevert, np.ndarray) and tevert.dtype == np.uint8:
            seq_mat[rev_idx, :150] = tevert
        timings['Evert'] = round(time.time() - t, 4)

        # PrimerSort
        t = time.time()
        ps_seq, ps_mat = PrimerSort(seq, F_primers, R_primers, mat_uint8=seq_mat)
        timings['PrimerSort'] = round(time.time() - t, 4)
        results['valid_reads'] = len(ps_seq)

        if len(ps_seq) == 0:
            return False, timings, {'error': 'no valid reads'}

        # KeyFind
        t = time.time()
        k = 156; n = k + 8; m = len(_de2bi(n))
        KeyMap, KmatDNA = KeyFind(ps_seq, ps_mat_u8=ps_mat)
        OligoLen, idxlen, ncolpo, nbitpo, ncolco, nbitco = Keysort(KeyMap, KmatDNA, k, filename)
        timings['KeyFind'] = round(time.time() - t, 4)
        results['OligoLen'] = OligoLen

        # Consensus
        t = time.time()
        MIndex, sorted_u8, boundaries = Consensus(ps_seq, idxlen, OligoLen, ps_mat_u8=ps_mat)
        Depth = np.array([MIndex[i].shape[0] if MIndex[i] is not None else 0 for i in range(len(MIndex))])
        timings['Consensus'] = round(time.time() - t, 4)
        results['mean_depth'] = round(float(np.mean(Depth)), 1)
        results['min_depth'] = int(np.min(Depth))

        # Binarization + RS
        t = time.time()
        _init_rs_reedsolo()
        SeqtoBin = Binarization(MIndex, sorted_u8=sorted_u8, boundaries=boundaries)
        rsECC, RSerrinf = RSdecoding(SeqtoBin, k, m, n)
        rsECC[RSerrinf < 0, :] = 0
        rsECC[RSerrinf > 4, :] = 0
        for a in range(min(OligoLen, len(rsECC))):
            if np.any(rsECC[a, :] > 1):
                rsECC[a, :] = 0
        payload = rsECC[:, :k - idxlen]
        Temp_mat = payload[:-1, :].copy()
        timings['RS'] = round(time.time() - t, 4)
        results['rs_uncorrectable'] = int(np.sum(RSerrinf < 0))
        results['rs_corrected'] = int(np.sum(RSerrinf > 0))

        # Null + XOR
        t = time.time()
        term = ncolpo * 3 * nbitpo + ncolco * 3 * nbitco
        Nullkey_1 = term + 1 if term % 2 == 1 else term
        s_val = Temp_mat.shape[1]
        Nullkey_2 = Nullkey_1 // 2
        Null_count = 0; nk1 = Nullkey_1
        while nk1 % s_val != 0: nk1 += 1; Null_count += 1
        Null_xor_count = 0; nk2 = Nullkey_2
        while nk2 % s_val != 0: nk2 += 1; Null_xor_count += 1

        data_Nd, xor_Nd, dNc, xNc = NullFind(seq, Null_count, Null_xor_count, s_val, F_primers, seq_mat_u8=seq_mat)
        bindataN = _dna_to_binary(data_Nd[20:])
        binxorN = _dna_to_binary(xor_Nd[20:])

        matsz = k - idxlen
        Temp_mat = payload[:-1, :].copy()
        data_row = Temp_mat.size * 2 // (3 * s_val) - 1
        if data_row < Temp_mat.shape[0]:
            il = min(len(bindataN), Temp_mat.shape[1])
            Temp_mat[data_row, :il] = bindataN[:il]
        if Temp_mat.shape[0] > 0:
            il = min(len(binxorN), Temp_mat.shape[1])
            Temp_mat[-1, :il] = binxorN[:il]

        # Recovery
        success = False
        for stage in ['Global', 'Depth', 'Local']:
            try:
                if stage == 'Global':
                    ECC_tri = Temp_mat[:Temp_mat.size * 2 // (3 * s_val), :]
                    ECC_flat = ECC_tri.flatten(order='C')
                    if Null_count > 0: ECC_flat = ECC_flat[:-Null_count]
                    tc, tp = tridecode(ECC_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
                    trid = _make_triangulation(tc, tp)
                    success = True
                elif stage == 'Depth':
                    Depthcut = 16
                    ErrD = np.where(Depth <= Depthcut)[0]
                    if len(ErrD) > 0:
                        err_mask = RSerrinf[ErrD] != 0
                        err_indices = ErrD[err_mask]
                        if len(err_indices) > 0:
                            Temp_mat[err_indices, :] = 0
                    xtp, md, tri_p, xor_p, XA, XB, XC = XORrun(Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)
                    ECC_flat = xtp.flatten(order='C')
                    if Null_count > 0: ECC_flat = ECC_flat[:-Null_count]
                    tc, tp = tridecode(ECC_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
                    trid = trierrXOR(Temp_mat, tc, tp, Null_count, Null_xor_count, ncolpo, nbitpo, ncolco, nbitco)
                    success = True
                elif stage == 'Local':
                    ErrD = np.where(Depth <= 16)[0]
                    err_mask = RSerrinf[ErrD] != 0 if len(ErrD) > 0 else np.array([], dtype=bool)
                    err_indices = ErrD[err_mask] if len(err_mask) > 0 else np.array([], dtype=int)
                    if len(err_indices) > 0:
                        Temp_mat[err_indices, :] = 0
                    xtp, md, tri_p, xor_p, XA, XB, XC = XORrun(Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)
                    if len(ErrD) > 0:
                        XORMIndex, Del_xor = nestedXOR(MIndex, xtp, ErrD, RSerrinf, XA, XB, XC, payload, tri_p, xor_p, matsz, Null_count)
                        Temp_mat = reRSdec(XORMIndex, Temp_mat, k, m, n, Del_xor, matsz)
                    xtp2, _, _, _, _, _, _ = XORrun(Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)
                    ECC_flat = xtp2.flatten(order='C')
                    if Null_count > 0: ECC_flat = ECC_flat[:-Null_count]
                    tc, tp = tridecode(ECC_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
                    trid = trierrXOR(Temp_mat, tc, tp, Null_count, Null_xor_count, ncolpo, nbitpo, ncolco, nbitco)
                    success = True
                if success:
                    results['recovery_stage'] = stage
                    break
            except Exception:
                continue

        timings['XOR'] = round(time.time() - t, 4)

        if success:
            if hasattr(trid, 'vertices'):
                results['vertices'] = len(trid.vertices)
                results['faces'] = len(trid.faces)
            elif isinstance(trid, dict):
                results['vertices'] = len(trid.get('vertices', []))
                results['faces'] = len(trid.get('faces', []))
        results['success'] = success

        timings['total'] = round(sum(timings.values()), 3)
        return success, timings, results

    except Exception as e:
        results['error'] = str(e)
        results['success'] = False
        timings['total'] = round(sum(timings.values()), 3)
        return False, timings, results


def main():
    print("=" * 70)
    print("다운샘플링 실험: J6 (Jikji6)")
    print("=" * 70)

    # FASTQ 로드
    t0 = time.time()
    seq_all = read_fastq_file(os.path.join(DATA_DIR, 'J6.fastq'))
    print(f"FASTQ 로드: {len(seq_all):,} reads ({round(time.time()-t0, 1)}s)")

    # 원본 coverage 추정 (full run에서 mean depth 사용)
    original_coverage = 1421  # J6 실측값

    all_results = []

    for target_cov in DOWNSAMPLE_TARGETS:
        print(f"\n{'─'*50}")
        print(f"Coverage: {target_cov}x")

        # 다운샘플
        ds_seq = downsample_reads(seq_all, target_cov, original_coverage)
        print(f"  리드: {len(ds_seq):,} ({len(ds_seq)/len(seq_all)*100:.1f}%)")

        # 실행
        ok, timings, results = run_pipeline(ds_seq, FILENAME)

        results['target_coverage'] = target_cov
        results['n_reads'] = len(ds_seq)
        results['timings'] = timings

        print(f"  성공: {ok}")
        if ok:
            print(f"  실제 Depth: {results.get('mean_depth', '?')}x")
            print(f"  RS 수정불가: {results.get('rs_uncorrectable', '?')}")
            print(f"  복구 단계: {results.get('recovery_stage', '?')}")
            print(f"  꼭짓점/면: {results.get('vertices', '?')}/{results.get('faces', '?')}")
            print(f"  시간: {timings.get('total', '?')}s")
        else:
            print(f"  오류: {results.get('error', 'unknown')}")

        all_results.append(results)

    # 결과 저장
    output_file = os.path.join(OUTPUT_DIR, 'downsampling_experiment_J6.json')
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n결과 저장: {output_file}")

    # 요약 테이블
    print(f"\n{'='*70}")
    print(f"{'Coverage':>10} {'Reads':>12} {'Depth':>8} {'RS_err':>8} {'Stage':>8} {'V/F':>12} {'Time':>8}")
    print(f"{'─'*70}")
    for r in all_results:
        cov = r.get('target_coverage', '?')
        reads = r.get('n_reads', 0)
        depth = r.get('mean_depth', '?')
        rs_err = r.get('rs_uncorrectable', '?')
        stage = r.get('recovery_stage', 'FAIL')
        verts = r.get('vertices', '?')
        faces = r.get('faces', '?')
        t = r.get('timings', {}).get('total', '?')
        print(f"{cov:>10}x {reads:>12,} {depth:>8} {rs_err:>8} {stage:>8} {verts}/{faces:>6} {t:>7}s")


if __name__ == '__main__':
    main()
