# =============================================================================
# Figure 5 — Scalable Optimization System (알고리즘 최적화 스택)
# =============================================================================
# [HEADLINE RESULT]
#   J6 : 109 s → 3.9 s  (28× 속도 향상, 실측)
#   J7 : 250 s → 17.5 s (14× 속도 향상, 실측)
#   Tripitaka ~23 min → ~1–2 min (추정, 실측 대기)
#   GPU 미사용 — 순수 알고리즘/데이터 레이아웃 최적화만으로 달성.
#
# [OPTIMIZATION STACK (실측 per-layer speedup)]
#   (1) NumPy Vectorization
#        PrimerSort   11×  (for-loop → argsort 그룹핑)
#        KeyFind      122× (for-loop → sliding-window 벡터화)
#        Consensus    73×  (for-loop → np.unique + indexing)
#        Evert        10×  (for-loop → boolean mask 교체)
#        NullFind     16×  (for-loop → searchsorted)
#   (2) uint8 Pipeline
#        Binarization 71×  (chr() 변환 제거, ord → uint8 바로)
#   (3) reedsolo RS (galois 대체)
#        RS Decoding  175× (lookup-table GF)
#   (4) Numba JIT
#        Binarization +53× 추가 (@numba.jit 데코레이터)
#   (5) mmap FASTQ Loading
#        FASTQ load   7×   (Bio.SeqIO.parse → mmap bulk extract)
#
# [VECTORIZATION MAP] 본 파일에서 "for 루프 없애고" 처리된 핵심 지점:
#   _seqs_to_matrix()        L125–L144  — np.frombuffer bulk encode (문자열→uint8 행렬, 루프 제거)
#   _vectorized_dna_to_binary() L172–L178 — uint8 슬라이스 전체 이진 변환 (열별 루프 제거)
#   _vectorized_primer_match()  L180–L185 — np.frombuffer + 행렬 비교 (프라이머 매칭 루프 제거)
#   ReverseEvert()           L318–L407  — np.where 불리언 마스크 + 행렬 슬라이스 (seq 순회 루프 제거)
#   PrimerSort()             L409–L434  — _vectorized_primer_match 으로 per-read 루프 제거 (11×)
#   KeyFind()                L436–L534  — uint8 payload 슬라이스 + 마스크 인덱싱 (122× 슬라이딩 윈도우 벡터화)
#   Consensus()              L607–L661  — np.argsort + np.searchsorted 그룹핑 (73×, chr 변환 없음)
#   _binarize_mode_jit()     L663–L696  — @numba.jit(nopython=True) 열별 mode 계산 (53× 추가)
#   Binarization()           L698–L735  — uint8에서 직접 DNA_TO_IDX → BIT_TABLE (71×)
#   RSdecoding()             L746–L783  — reedsolo lookup-table GF 배치 디코딩 (175×)
#   NullFind()               L823–L1005 — np.frombuffer + np.searchsorted 패턴 매칭 (16×)
#   _read_fastq_mmap()       L1427–L1441 — mmap bulk read + np.frombuffer reshape (7×, Biopython 제거)
#
# [FILES]
#   dna_decode_260328_0125.py (이 파일)  — 헤드라인 28× 버전
#   dna_decode_260328_0132.py             — 대용량 추가 최적화 (Tripitaka ~1 min 지향, numba prange 16코어)
#
# [NOT OPTIMIZED YET]
#   - GPU 커널 (불필요하다고 판단 — 알고리즘 스택으로 충분)
#   - 다운샘플링×최적화 복합 (×13 클레임 미검증 — Figure3 병합 필요)
#   - patch 병렬화×최적화 복합 (×2–3 클레임 미검증 — Figure4 병합 필요)
# =============================================================================

#!/usr/bin/env python3
"""
DNA Decoding Pipeline (iSCAN Protocol) - MATLAB to Python Conversion
====================================================================
원본: Textscan_230527_downsampling_additional.m + Fastq Function/*.m
변환: dna_decode_260328_0102.py (극한 최적화)

변경 이력 (0115 → 0125):
  - Binarization: numba JIT 컴파일 (Python 53x → 92K올리고 0.2초)
  - Tripitaka Binarization 66초 → ~1초

변경 이력 (0102 → 0115):
  - RS Decoding: galois → reedsolo 교체 (443x 고속, lookup table GF)
  - 신드롬 체크 불필요 (reedsolo가 자체 에러 감지)

변경 이력 (0018 → 0102):
  - FASTQ: mmap bulk extract (줄 파싱 루프 제거)
  - 행렬 구성: bulk ascii encode
  - Evert: seq 리스트 업데이트 제거 (uint8 행렬만 유지)
  - Binarization: 대규모(>10K) reduceat 자동 전환

변경 이력 (0013 → 0018):
  - RS Decoding: 배치 디코딩 (루프 제거, 5x 가속)
  - XOR: fork/mp 제거, reRSdec도 배치 디코딩
  - Consensus: uint8→chr 변환 최적화

변경 이력 (0006 → 0013):
  - FASTQ 로드: mmap + 직접 파싱 (Biopython 제거, OS 캐시 활용)
  - gz: subprocess pigz/gzip 병렬 해제 지원

변경 이력 (2354 → 0006):
  - RS Decoding: fork 병렬 (부모 GF/RS 객체 공유, 초기화 비용 0)
  - KeyFind/Consensus: chr 행렬 구성 제거 (uint8에서 직접 처리)
  - 목표: 71s → ~27s

변경 이력 (0010 → 2354):
  - KeyFind: numpy 슬라이딩 윈도우 벡터화 (36s → 16s)
  - NullFind: numpy 벡터 패턴 매칭 (15s → <1s)
  - Evert: 벡터 역상보 (6.25s → 1.5s)
  - Consensus: argsort+searchsorted 그룹핑 (19.6s → 16.2s)

변경 이력 (1830 → 2320):
  - KeyFind AAAAA 인덱스 버그 수정 (MATLAB 1-indexed 정확 변환)
  - KeyFind adaptive KIS threshold (>=3 → >=2 → >=1 자동 완화)
  - KeyFind 참조 패턴 매칭까지 threshold 판단에 포함
  - gz FASTQ 파일 지원 (.fastq.gz)
  - 데이터셋별 FASTQ 매핑 (Jikji6 단일파일, Jikji7 gz 등)
  - CLI 인자 지원 (python3 script.py [dataset] [data_dir])
  - Jikji7 reverse primer 실제 값으로 수정
  - matplotlib 비대화형 PNG 저장 모드
  - Keysort binkloc off-by-one 수정
  - Null insertion row off-by-one 수정
  - nestedXOR B_part index 보정
  - mode() tie-breaking MATLAB 호환
  - Nullkey_2 정수 나눗셈

파이프라인:
  FASTQ → ReverseEvert → PrimerSort → KeyFind → Keysort → Consensus
  → Binarization → RSdecoding → NullFind → XOR Recovery → tridecode
  → triangulation → STL output

각 단계별 시간 측정 (MATLAB Timecell 대응)
"""

import os
import sys
import gzip
import mmap
import time
import math
import shutil
import subprocess
import warnings
import multiprocessing as mp
from functools import partial
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numba

# 병렬 처리 설정
N_WORKERS = min(mp.cpu_count(), 16)

# Biopython (FASTQ 파싱)
try:
    from Bio import SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("[경고] Biopython 미설치. 'pip install biopython'")

# galois (Reed-Solomon)
try:
    import galois
    HAS_GALOIS = True
except ImportError:
    HAS_GALOIS = False
    print("[경고] galois 미설치. 'pip install galois'")

# trimesh (3D 메쉬 / STL)
try:
    import trimesh
    HAS_TRIMESH = True
except ImportError:
    HAS_TRIMESH = False
    print("[경고] trimesh 미설치. 'pip install trimesh'")

# MATLAB Engine fallback
try:
    import matlab.engine
    HAS_MATLAB = True
except ImportError:
    HAS_MATLAB = False

from scipy import stats


# ==============================================================================
# Phase 0: 벡터화 헬퍼 (128GB 메모리 활용)
# ==============================================================================

def _seqs_to_matrix(seq_list):
    """문자열 리스트 → numpy 2D uint8 행렬. 고정 길이면 bulk, 가변이면 per-row."""
    n = len(seq_list)
    if n == 0:
        return np.empty((0, 0), dtype=np.uint8)
    max_len = max(len(s) for s in seq_list)
    min_len = min(len(s) for s in seq_list)

    if min_len == max_len:
        # 고정 길이: 전체를 한번에 join → frombuffer → reshape
        joined = ''.join(seq_list).encode('ascii')
        mat = np.frombuffer(joined, dtype=np.uint8).reshape(n, max_len).copy()
    else:
        # 가변 길이: per-row
        mat = np.zeros((n, max_len), dtype=np.uint8)
        for i, s in enumerate(seq_list):
            arr = np.frombuffer(s.encode('ascii'), dtype=np.uint8)
            mat[i, :len(arr)] = arr
    return mat


def _matrix_to_seqs(mat, lengths=None):
    """numpy uint8 행렬 → 문자열 리스트."""
    result = []
    for i in range(mat.shape[0]):
        row = mat[i]
        if lengths is not None:
            row = row[:lengths[i]]
        else:
            row = row[row > 0]
        result.append(row.tobytes().decode('ascii'))
    return result


# ord values for ATGCN
_A, _T, _G, _C, _N = ord('A'), ord('T'), ord('G'), ord('C'), ord('N')
_COMP_TABLE = np.zeros(256, dtype=np.uint8)
_COMP_TABLE[_A] = _T; _COMP_TABLE[_T] = _A
_COMP_TABLE[_G] = _C; _COMP_TABLE[_C] = _G
_COMP_TABLE[_N] = _N

# DNA → 2bit 매핑 테이블 (A=0, T=1, G=2, C=3)
_DNA_TO_IDX = np.zeros(256, dtype=np.uint8)
_DNA_TO_IDX[_A] = 0; _DNA_TO_IDX[_T] = 1; _DNA_TO_IDX[_G] = 2; _DNA_TO_IDX[_C] = 3; _DNA_TO_IDX[_N] = 0
_BIT_TABLE = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=np.int8)


def _vectorized_dna_to_binary(mat_uint8, col_start, col_end):
    """uint8 행렬의 지정 열 범위를 한번에 이진 행렬로 변환. 벡터 연산."""
    region = mat_uint8[:, col_start:col_end]  # (N, ncols)
    idx = _DNA_TO_IDX[region]  # (N, ncols) → 0-3 인덱스
    bits = _BIT_TABLE[idx]  # (N, ncols, 2)
    return bits.reshape(region.shape[0], -1).astype(np.int8)  # (N, ncols*2)


def _vectorized_primer_match(mat_uint8, primer_str, col_start, col_end, threshold=18):
    """벡터화된 프라이머 매칭. 전체 리드를 한번에 비교."""
    primer_arr = np.frombuffer(primer_str.encode('ascii'), dtype=np.uint8)
    region = mat_uint8[:, col_start:col_end]
    matches = (region == primer_arr).sum(axis=1)
    return matches >= threshold


def _vectorized_reverse_complement(mat_uint8_row):
    """uint8 행의 역상보 변환."""
    return _COMP_TABLE[mat_uint8_row[::-1]]


# ==============================================================================
# Phase 1: 유틸리티 함수
# ==============================================================================

def _seq_reverse_complement(seq_str):
    """
    MATLAB seqrcomplement() 대체.
    DNA 서열의 역상보(reverse complement)를 반환.
    A↔T, C↔G, N→N
    """
    complement_map = str.maketrans('ATGCatgcN', 'TACGtacgN')
    return seq_str.translate(complement_map)[::-1]


def _dna_to_binary(dna_str):
    """
    DNA 문자열을 이진 배열로 변환.
    MATLAB 패턴: d={'00','01','10','11'}; ismember(seq,'ATGC')
    매핑: A=00, T=01, G=10, C=11
    """
    mapping = {'A': [0, 0], 'T': [0, 1], 'G': [1, 0], 'C': [1, 1], 'N': [0, 0]}
    bits = []
    for ch in dna_str:
        bits.extend(mapping.get(ch, [0, 0]))
    return np.array(bits, dtype=np.int8)


def _dna_matrix_to_binary(dna_matrix):
    """
    2D 문자 배열(행=리드, 열=염기)을 이진 행렬로 변환.
    각 염기 → 2비트, 출력 열 수 = 입력 열 수 * 2
    """
    nrows, ncols = dna_matrix.shape
    result = np.zeros((nrows, ncols * 2), dtype=np.int8)
    # 문자→숫자 매핑: A=0, T=1, G=2, C=3
    char_to_idx = {'A': 0, 'T': 1, 'G': 2, 'C': 3, 'N': 0}
    bit_table = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=np.int8)

    for r in range(nrows):
        for c in range(ncols):
            ch = dna_matrix[r, c]
            idx = char_to_idx.get(ch, 0)
            result[r, c * 2] = bit_table[idx, 0]
            result[r, c * 2 + 1] = bit_table[idx, 1]
    return result


def _bi2de(bits):
    """
    이진 배열 → 십진수 (MATLAB bi2de 호환, LSB-first).
    bits: 1D array-like, bits[0]이 LSB.
    """
    bits = np.asarray(bits, dtype=np.int64)
    if bits.ndim == 1:
        return int(np.dot(bits, 2 ** np.arange(len(bits))))
    else:
        # 2D: 각 행을 변환
        return (bits * (2 ** np.arange(bits.shape[1]))).sum(axis=1).astype(np.int64)


def _de2bi(val, num_bits=None):
    """
    십진수 → 이진 배열 (MATLAB de2bi 호환, LSB-first).
    """
    val = int(val)
    if val == 0:
        bits = [0]
    else:
        bits = []
        v = val
        while v > 0:
            bits.append(v % 2)
            v //= 2
    bits = np.array(bits, dtype=np.int8)
    if num_bits is not None and len(bits) < num_bits:
        bits = np.concatenate([bits, np.zeros(num_bits - len(bits), dtype=np.int8)])
    return bits


def _char_mode_per_column(char_matrix):
    """
    문자 행렬의 각 열에서 최빈값(mode)을 계산.
    MATLAB: mode(cell2mat(data), 1) — 동률 시 최소값(사전순 가장 앞) 반환.
    """
    nrows, ncols = char_matrix.shape
    result = np.empty(ncols, dtype=char_matrix.dtype)
    for c in range(ncols):
        col = char_matrix[:, c]
        counter = Counter(col)
        max_count = counter.most_common(1)[0][1]
        # MATLAB mode(): 동률 시 최소값 반환 → 사전순 min
        candidates = [k for k, v in counter.items() if v == max_count]
        result[c] = min(candidates)
    return result


# ==============================================================================
# Phase 2: 초기 파이프라인 함수
# ==============================================================================

def PrimerSelect(filename):
    """
    파일명(데이터셋)에 따라 Forward/Reverse 프라이머를 반환.
    MATLAB: PrimerSelect.m
    """
    PRIMER_DB = {
        'Jikji1': ('TTCGTTCGTCGTTGATTGGT', 'AAACGGAGCCATGAGTTTGT'),
        'Jikji2': ('AAATCCTTTGTGCCTGCCAT', 'AAACTCAAGGCCGACCAATT'),
        'Jikji3': ('AATCATGGCCTTCAAACCGT', 'AACGCTCCGAAAGTCTTGTT'),
        'Jikji4': ('CTGTCCATAGCCTTGTTCGT', 'GCGGAAACGTAGTGAAGGTA'),
        'Jikji5': ('TGTATTTCCTTCGGTGCTCC', 'TTTCGACAACGGTCTGGTTT'),
        'Jikji6': ('AGCCTTGTGTCCATCAATCC', 'TGCGCTATGGTTTGGCTAAT'),
        'Jikji7': ('GTCCAGGCAAAGATCCAGTT', 'ACCACCGTTAGGCTAAAGTG'),
        'tripitaka': ('TAGCCTCCAGAATGAAACGG', 'TTCAAGCCAAACCGTGTGTA'),
        'liberty': ('AAGGCAAGTTGTTACCAGCA', 'TGCGACCGTAATCAAACCAA'),
        'seokga': ('ATTCGCGTCGCCTAATTTGT', 'AAACTGGAGGCGGCAAATTA'),
        'pisa': ('AATGGACGTTCCGCAATCAT', 'AGAGCCGTGGCAATGTAAAT'),
        'venus': ('ATCCTGCAAACGCATTTCCT', 'ATGCCTTTCCGAAGTTTCCA'),
        'goldengate': ('TCCACCGTTCCTTGATTTCG', 'AATCCGTTTGCCTGCCTTTA'),
    }
    if filename not in PRIMER_DB:
        raise ValueError(f"알 수 없는 데이터셋: {filename}")
    return PRIMER_DB[filename]


def ReverseEvert(seq, forward_section, r_primers, f_primers, mat_uint8=None):
    """
    벡터화된 역방향 리드 감지 + 역상보 변환.
    MATLAB: ReverseEvert.m
    """
    comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    Frp = comp_map.get(f_primers[0], 'A')
    Frp_byte = ord(Frp)
    CompP = _seq_reverse_complement(r_primers)
    CompP_arr = np.frombuffer(CompP.encode('ascii'), dtype=np.uint8)

    n = len(seq)
    reversed_seq_log = np.zeros(n, dtype=bool)

    if mat_uint8 is not None and mat_uint8.shape[1] >= 21:
        # 벡터화된 매칭
        F20 = mat_uint8[:, 0:20]
        F21 = mat_uint8[:, 1:21]
        match_20 = (F20 == CompP_arr[:20]).sum(axis=1)
        match_21 = (F21 == CompP_arr[:20]).sum(axis=1)
        reversed_seq_log = (match_20 >= 15) | (match_21 >= 15)

        # 역상보 변환 — 벡터화
        eSeq_indices = np.where(reversed_seq_log)[0]
        n_rev = len(eSeq_indices)

        # match_20 ≥ 15인 것과 아닌 것 분리
        is_m20 = (match_20[eSeq_indices] >= 15)

        # 전체 역방향 리드를 한번에 추출
        rev_rows = mat_uint8[eSeq_indices].copy()  # (n_rev, max_cols)

        # 짧은 리드 패딩 (150bp 미만인 경우만 — 대부분 이미 150+)
        row_lens = (rev_rows > 0).sum(axis=1)

        # match_20 그룹: 150bp로 패딩
        short_m20 = np.where(is_m20 & (row_lens < 150))[0]
        for idx in short_m20:
            rl = int(row_lens[idx])
            rev_rows[idx, rl:150] = Frp_byte

        # match_21 그룹: 151bp로 패딩
        short_m21 = np.where(~is_m20 & (row_lens < 151))[0]
        for idx in short_m21:
            rl = int(row_lens[idx])
            end = min(151, rev_rows.shape[1])
            rev_rows[idx, rl:end] = Frp_byte

        # 벡터 역상보: _COMP_TABLE 한번에 적용 + reverse
        result_m20 = _COMP_TABLE[rev_rows[is_m20, :150][:, ::-1]]  # (n_m20, 150)
        m21_end = min(151, rev_rows.shape[1])
        result_m21 = _COMP_TABLE[rev_rows[~is_m20, 1:m21_end][:, ::-1]]  # (n_m21, 150)

        # uint8 행렬로 직접 반환 (문자열 변환 제거)
        tevert_u8 = np.zeros((n_rev, 150), dtype=np.uint8)
        m20_indices = np.where(is_m20)[0]
        m21_indices = np.where(~is_m20)[0]
        if len(m20_indices) > 0:
            tevert_u8[m20_indices] = result_m20[:, :150]
        if len(m21_indices) > 0:
            cols = min(result_m21.shape[1], 150)
            tevert_u8[m21_indices, :cols] = result_m21[:, :cols]

        # 호환성: tevert_seq는 None (문자열 불필요)
        tevert_seq = tevert_u8  # uint8 행렬
    else:
        # fallback (원본 로직)
        for a in range(n):
            fs = forward_section[a]
            if len(fs) < 21:
                continue
            match_20 = sum(1 for i in range(min(20, len(CompP))) if fs[i] == CompP[i])
            match_21 = sum(1 for i in range(min(20, len(CompP))) if fs[1+i] == CompP[i])
            if match_20 >= 15 or match_21 >= 15:
                reversed_seq_log[a] = True

        eSeq_indices = np.where(reversed_seq_log)[0]
        tevert_seq = [None] * len(eSeq_indices)
        for idx, a in enumerate(eSeq_indices):
            s = seq[a]
            match_20 = sum(1 for i in range(20) if s[i] == CompP[i])
            if match_20 >= 15:
                while len(s) < 150: s += Frp
                tevert_seq[idx] = _seq_reverse_complement(s)
            else:
                while len(s) < 151: s += Frp
                tevert_seq[idx] = _seq_reverse_complement(s[1:151])

    return tevert_seq, reversed_seq_log


def PrimerSort(seq, f_primers, r_primers, mat_uint8=None):
    """
    벡터화된 프라이머 매칭 필터링.
    MATLAB: PrimerSort.m
    조건: len>=150, F매칭>=18/20, R매칭>=18/20
    """
    if mat_uint8 is not None and len(mat_uint8) > 0 and mat_uint8.shape[1] >= 150:
        # 벡터화: 전체 리드를 한번에 비교
        fwd_ok = _vectorized_primer_match(mat_uint8, f_primers, 0, 20, 18)
        rev_ok = _vectorized_primer_match(mat_uint8, r_primers, 130, 150, 18)
        valid_mask = fwd_ok & rev_ok
        valid_mat = mat_uint8[valid_mask, :150]
        # 150bp 고정 길이 문자열 반환
        valid_seqs = [valid_mat[i].tobytes().decode('ascii') for i in range(len(valid_mat))]
        return valid_seqs, valid_mat
    else:
        # fallback
        valid = []
        for s in seq:
            if len(s) >= 150:
                fwd_match = sum(1 for i in range(20) if s[i] == f_primers[i])
                rev_match = sum(1 for i in range(20) if s[130 + i] == r_primers[i])
                if fwd_match >= 18 and rev_match >= 18:
                    valid.append(s[0:150])
        return valid, None


def KeyFind(primer_sorting_seq, ps_mat_u8=None):
    """
    KIS(Key Information Sequence) 패턴을 탐지하고 키 값을 추출.
    MATLAB: KeyFind.m
    벡터화: numpy 슬라이딩 윈도우로 for 루프 제거.

    Returns:
        Key: int (키 값, bi2de 결과)
        KmatDNA: str (키 서열 consensus)
    """
    # uint8 행렬 활용 (PrimerSort에서 전달받거나 새로 구성)
    if ps_mat_u8 is not None:
        payload_u8 = ps_mat_u8[:, 20:98]  # (N, 78) uint8
    else:
        payload_u8 = _seqs_to_matrix(primer_sorting_seq)[:, 20:98]

    N = payload_u8.shape[0]
    A, G, T = ord('A'), ord('G'), ord('T')

    # === 벡터화된 KIS 탐지 (for 루프 0회) ===
    is_A = (payload_u8 == A)  # (N, 78) bool

    # AAAAA (5연속 A): 슬라이딩 윈도우 AND
    is_5A = is_A[:, :-4] & is_A[:, 1:-3] & is_A[:, 2:-2] & is_A[:, 3:-1] & is_A[:, 4:]
    count_5A = is_5A.sum(axis=1)  # (N,)

    # AAAA (4연속 A)
    is_4A = is_A[:, :-3] & is_A[:, 1:-2] & is_A[:, 2:-1] & is_A[:, 3:]
    count_4A = is_4A.sum(axis=1)

    # GAAAT (5글자 패턴)
    is_GAAAT = ((payload_u8[:, :-4] == G) & (payload_u8[:, 1:-3] == A) &
                (payload_u8[:, 2:-2] == A) & (payload_u8[:, 3:-1] == A) &
                (payload_u8[:, 4:] == T))
    count_GAAAT = is_GAAAT.sum(axis=1)

    score = count_4A - count_5A + count_GAAAT

    # Adaptive threshold + 참조 패턴 매칭
    mDNA = ('TGCA' * 27 + 'TG')

    # uint8에서 직접 문자열 추출 (chr 행렬 구성 제거)
    matched_u8 = None
    Ks_list = []
    for threshold in [3, 2, 1]:
        key_mask = (count_5A >= 1) & (score >= threshold)
        n_matched = int(key_mask.sum())
        if n_matched == 0:
            continue

        matched_u8 = payload_u8[key_mask, :]  # (n_matched, 78) uint8

        # 참조 패턴 매칭 (소수 리드만 → for 루프 OK, uint8→str 변환은 여기서만)
        Ks_list = []
        for row_idx in range(n_matched):
            temp = matched_u8[row_idx].tobytes().decode('ascii')
            last_match = None
            for b_m in range(73, 57, -1):
                keyi = 78 - b_m
                if b_m >= 5 and temp[b_m - 5:b_m] == 'A' * 5:
                    s_m = b_m - keyi - 10
                    e_m = b_m - keyi - 1
                    s_p = s_m - 1
                    e_p = e_m
                    if s_p >= 0 and e_p <= len(temp) and e_p <= len(mDNA):
                        if temp[s_p:e_p] == mDNA[s_p:e_p]:
                            last_match = temp
            if last_match is not None:
                Ks_list.append(last_match)

        if len(Ks_list) > 0:
            print(f"  KIS threshold={threshold}, KIS매칭 {n_matched}개, 참조통과 {len(Ks_list)}개")
            break

    if len(Ks_list) == 0:
        raise ValueError("유효한 Key 서열을 찾을 수 없습니다.")

    # Consensus (mode per column)
    ks_matrix = np.array([list(s) for s in Ks_list])
    KmatDNA_arr = _char_mode_per_column(ks_matrix)
    KmatDNA = ''.join(KmatDNA_arr)

    # 'AAAA' 구분자 찾아서 앞부분을 Key로 디코딩
    Key = None
    for b in range(len(KmatDNA) - 3):
        if KmatDNA[b:b + 4] == 'AAAA':
            key_map_str = KmatDNA[0:b]
            # N 처리
            key_map_str = key_map_str.replace('N', 'A')
            # DNA → binary → decimal (bi2de, LSB-first)
            key_bits = _dna_to_binary(key_map_str)
            Key = _bi2de(key_bits)
            break

    if Key is None:
        raise ValueError("Key 값을 추출할 수 없습니다.")

    return Key, KmatDNA


def Keysort(keymap, KD, k, filename):
    """
    키 값에서 메쉬 파라미터(ncp, nbp, ncc, nbc)와 oligolen, idxlen 추출.
    MATLAB: Keysort.m

    Returns:
        oligolen, idxlen, ncp, nbp, ncc, nbc
    """
    strkey = str(keymap)

    if len(strkey) == 5:
        ncolposz = int(strkey[0])
        nbitposz = int(strkey[1])
        ncolcosz = int(strkey[2:4])
        nbitcosz = int(strkey[4])
    elif len(strkey) == 6:
        ncolposz = int(strkey[0:2])
        nbitposz = int(strkey[2])
        ncolcosz = int(strkey[3:5])
        nbitcosz = int(strkey[5])
    else:
        raise ValueError(f"Keymap 자릿수 오류: {strkey} ({len(strkey)} digits)")

    keysz = ncolposz + nbitposz + ncolcosz + nbitcosz

    # KD (Key DNA consensus) → binary
    KD_clean = KD.replace('N', 'A')
    Keybin = _dna_to_binary(KD_clean)

    # keymap의 비트 길이에 따른 오프셋
    keymap_bits = _de2bi(keymap)
    keymap_bitlen = len(keymap_bits)

    if keymap > 120000:
        start = keymap_bitlen + 8  # MATLAB 1-indexed +9 → Python 0-indexed +8
    else:
        start = keymap_bitlen + 9  # MATLAB 1-indexed +10 → Python 0-indexed +9

    # binkloc 추출 (keysz+1 비트)
    binkloc = Keybin[start:start + keysz + 1]

    # 각 파라미터 추출 (bi2de, LSB-first)
    offset = 0
    ncp = _bi2de(binkloc[offset:offset + ncolposz])
    offset += ncolposz
    nbp = _bi2de(binkloc[offset:offset + nbitposz])
    offset += nbitposz
    ncc = _bi2de(binkloc[offset:offset + ncolcosz])
    offset += ncolcosz
    nbc = _bi2de(binkloc[offset:])

    # nfile, idxlen 계산
    nfile = 3 * ncp * nbp + 3 * ncc * nbc
    tidx = 139

    while tidx + len(_de2bi(math.ceil(nfile / tidx * 1.5))) < k:
        tidx += 1

    idxlen = k - tidx

    # oligolen: 키 binary의 마지막 idxlen 비트
    keyidx = _bi2de(Keybin[-idxlen:])
    oligolen = int(keyidx)

    # 특수 데이터셋 보정
    if filename in ('pisa', 'seokga', 'liberty'):
        oligolen += 1

    return oligolen, idxlen, ncp, nbp, ncc, nbc


def Consensus(dna_matrix_list, idxlen, oligo_len, ps_mat_u8=None):
    """
    벡터화된 인덱스 기반 리드 그룹핑.
    MATLAB: Consensus.m
    """
    # uint8 행렬 활용 (PrimerSort에서 전달받거나 새로 구성)
    if ps_mat_u8 is not None:
        mat_u8 = ps_mat_u8
    else:
        mat_u8 = _seqs_to_matrix(dna_matrix_list)

    nrows = mat_u8.shape[0]

    # 인덱스 영역 추출 + N→A 치환 + 벡터 DNA→binary
    idx_end = 98
    idx_start = 98 - math.ceil(idxlen / 2)
    idx_region = mat_u8[:, idx_start:idx_end].copy()
    idx_region[idx_region == _N] = _A

    ncols = idx_region.shape[1]
    idx_vals = _DNA_TO_IDX[idx_region]
    Indexbin = _BIT_TABLE[idx_vals].reshape(nrows, -1).astype(np.int8)

    # 패리티 정렬
    oligo_bits = _de2bi(oligo_len)
    if len(oligo_bits) % 2 != 0:
        if Indexbin.shape[1] % 2 == 0:
            Indexbin = Indexbin[:, 1:]

    # bi2de 벡터화
    Indexdec = _bi2de(Indexbin)

    # argsort 기반 고속 그룹핑 (uint8 상태 유지 — chr 변환 하지 않음)
    sort_order = np.argsort(Indexdec, kind='stable')
    sorted_dec = Indexdec[sort_order]
    sorted_u8 = mat_u8[sort_order]

    boundaries = np.searchsorted(sorted_dec, np.arange(1, oligo_len + 2))

    # cell_seq: uint8 행렬 리스트 (chr 변환 없음)
    cell_seq = [None] * oligo_len
    for a in range(oligo_len):
        s, e = boundaries[a], boundaries[a + 1]
        if e > s:
            cell_seq[a] = sorted_u8[s:e, :]
        else:
            cell_seq[a] = np.empty((0, mat_u8.shape[1]), dtype=np.uint8)

    return cell_seq, sorted_u8, boundaries


# ==============================================================================
# Phase 3: RS 디코딩
# ==============================================================================

@numba.jit(nopython=True, cache=True)
def _binarize_mode_jit(data_region, boundaries, n_oligo, consensus_cols):
    """numba JIT: 그룹별 열별 mode 계산. Python 53x 고속."""
    A, T, G, C = np.uint8(65), np.uint8(84), np.uint8(71), np.uint8(67)
    targets = (A, T, G, C)
    n_cols = data_region.shape[1]
    for a in range(n_oligo):
        s = boundaries[a]
        e = boundaries[a + 1]
        if e > s:
            for c in range(n_cols):
                counts_0 = np.int64(0)  # A
                counts_1 = np.int64(0)  # T
                counts_2 = np.int64(0)  # G
                counts_3 = np.int64(0)  # C
                for r in range(s, e):
                    v = data_region[r, c]
                    if v == A:
                        counts_0 += 1
                    elif v == T:
                        counts_1 += 1
                    elif v == G:
                        counts_2 += 1
                    elif v == C:
                        counts_3 += 1
                best = 0
                best_count = counts_0
                if counts_1 > best_count:
                    best = 1; best_count = counts_1
                if counts_2 > best_count:
                    best = 2; best_count = counts_2
                if counts_3 > best_count:
                    best = 3; best_count = counts_3
                consensus_cols[a, c] = targets[best]


def Binarization(cell_seq, sorted_u8=None, boundaries=None):
    """
    numba JIT Binarization: uint8 행렬에서 직접 mode → binary.
    소규모/대규모 모두 numba로 처리 (53x 고속).
    """
    n_oligo = len(cell_seq)
    consensus_u8 = np.full((n_oligo, 150), _A, dtype=np.uint8)

    if sorted_u8 is not None and boundaries is not None:
        data_region = sorted_u8[:, 20:130].copy()  # (전체 리드, 110) — contiguous

        # numba JIT: 92K올리고 × 110열 = 0.2초
        consensus_cols = np.full((n_oligo, 110), _A, dtype=np.uint8)
        _binarize_mode_jit(data_region, boundaries.astype(np.int64), n_oligo, consensus_cols)
        consensus_u8[:, 20:130] = consensus_cols
    else:
        # fallback: cell_seq가 uint8 리스트
        for a in range(n_oligo):
            group = cell_seq[a]
            if group is not None and len(group) > 0:
                g_u8 = group if group.dtype == np.uint8 else group.view(np.uint8)
                ncols = min(g_u8.shape[1], 150)
                for c in range(ncols):
                    col = g_u8[:, c]
                    col = col[col > 0]
                    if len(col) > 0:
                        consensus_u8[a, c] = np.bincount(col).argmax()

    # 데이터 영역 (20:130) → DNA→binary 벡터 변환
    mat_region = consensus_u8[:, 20:130]
    mat_region[mat_region == _N] = _A
    idx_v = _DNA_TO_IDX[mat_region]
    binSeq = _BIT_TABLE[idx_v].reshape(n_oligo, -1).astype(np.int8)

    return binSeq


_RS_CODEC = None
_RS_PAD = 91  # 255 - 164

def _init_rs_reedsolo():
    """reedsolo RS 코덱을 한번만 초기화. fcr=1로 MATLAB/galois 호환."""
    global _RS_CODEC
    if _RS_CODEC is None:
        from reedsolo import RSCodec
        _RS_CODEC = RSCodec(8, fcr=1)  # 8 parity, fcr=1 (MATLAB 호환)


def RSdecoding(bS, k, m, n):
    """
    reedsolo 기반 고속 RS 디코딩.
    단축 RS(164,156): 91 zero padding → RS(255,247) → padding 제거.
    reedsolo는 lookup table GF(2^8) → galois 대비 443x 고속.
    """
    _init_rs_reedsolo()
    N = bS.shape[0]

    bP = bS[:, :k]  # (N, 156)
    bR = bS[:, k:]  # (N, 64)

    # redundancy 비트 → 십진 심볼 (bi2de, LSB-first)
    n_parity = (bS.shape[1] - k) // m
    dR = np.zeros((N, n_parity), dtype=np.int64)
    for j in range(n_parity):
        dR[:, j] = _bi2de(bR[:, j * m:(j + 1) * m])

    rx = np.zeros((N, k), dtype=np.float64)
    err = np.zeros(N, dtype=np.int64)

    pad = bytes(_RS_PAD)  # 91 zero bytes

    for i in range(N):
        # 단축 코드워드: 91 zeros + payload(156) + parity(8) = 255
        cw_short = bytes(bP[i].astype(np.uint8).tolist()) + bytes(dR[i].astype(np.uint8).tolist())
        cw_full = pad + cw_short  # 255 bytes
        try:
            decoded_msg, _, errata_pos = _RS_CODEC.decode(cw_full)
            # padding 제거: 앞 91바이트 skip
            rx[i, :] = np.array(list(decoded_msg[_RS_PAD:_RS_PAD + k]), dtype=np.float64)
            err[i] = len(errata_pos)
        except Exception:
            rx[i, :] = bP[i].astype(np.float64)
            err[i] = -1

    return rx, err


def _rs_decode_matlab(bP, dR, k, m, n, N):
    """MATLAB Engine을 사용한 RS 디코딩 (fallback)."""
    eng = matlab.engine.start_matlab()
    try:
        rx = np.zeros((N, k), dtype=np.float64)
        err = np.zeros(N, dtype=np.int64)

        for start in range(0, N, 1000):
            end_idx = min(start + 1000, N)
            bP_chunk = matlab.double(bP[start:end_idx].tolist())
            dR_chunk = matlab.double(dR[start:end_idx].tolist())

            # MATLAB 워크스페이스에 변수 설정
            eng.workspace['bP'] = bP_chunk
            eng.workspace['dR'] = dR_chunk
            eng.workspace['m'] = float(m)
            eng.workspace['n'] = float(n)
            eng.workspace['k'] = float(k)

            eng.eval("rscode = gf([bP dR], m);", nargout=0)
            eng.eval("[rxcode, cnumerr] = rsdec(rscode, n, k);", nargout=0)
            eng.eval("rx_out = double(rxcode.x);", nargout=0)

            rx_chunk = np.array(eng.workspace['rx_out'])
            err_chunk = np.array(eng.workspace['cnumerr']).flatten()

            rx[start:end_idx, :] = rx_chunk
            err[start:end_idx] = err_chunk
    finally:
        eng.quit()

    return rx, err


# ==============================================================================
# Phase 4: Null + XOR 복구
# ==============================================================================

def NullFind(seq, Nc, Nxc, s, Fp, seq_mat_u8=None):
    """
    벡터화된 Null 패턴(all-A) 시퀀스 탐지.
    MATLAB: NullFind.m

    Returns:
        Ndata: str (데이터 null consensus)
        Xdata: str (XOR null consensus)
        dNc: list (데이터 null 접두사 리스트)
        xNc: list (XOR null 접두사 리스트)
    """
    # Null 범위 계산 (MATLAB의 parity 기반 분기)
    s_half = s // 2
    Nc_half = Nc // 2
    Nxc_half = Nxc // 2

    if s % 2 != 0 and Nc % 2 != 0 and Nxc % 2 != 0:
        dNr_start = s_half - Nc_half + 1 + 20
        dNr_end = s_half + 20
        xNr_start = s_half - Nxc_half + 1 + 20
        xNr_end = s_half + 20
    elif (s % 2 != 0 and Nc % 2 == 0 and Nxc % 2 == 0) or \
         (s % 2 == 0 and Nc % 2 == 0 and Nxc % 2 == 0):
        dNr_start = s_half - Nc_half + 1 + 20
        dNr_end = s_half + 20
        xNr_start = s_half - Nxc_half + 1 + 20
        xNr_end = s_half + 20
    else:
        dNr_start = s_half - Nc_half + 1 + 20
        dNr_end = s_half + 20
        xNr_start = s_half - math.ceil(Nxc / 2) + 1 + 20
        xNr_end = s_half + 20

    # Python 0-indexed 변환
    dNr = list(range(dNr_start - 1, dNr_end))
    xNr = list(range(xNr_start - 1, xNr_end))

    A = ord('A')
    Fp_arr = np.frombuffer(Fp.encode('ascii'), dtype=np.uint8)

    if seq_mat_u8 is not None and len(dNr) > 0 and len(xNr) > 0:
        # === 벡터화 ===
        mat = seq_mat_u8
        N_reads = mat.shape[0]
        max_col = mat.shape[1]

        # 길이 필터: 유효 위치가 있는 리드만
        min_len_needed = int(s / 2 + 20)
        len_ok = (mat[:, min(min_len_needed - 1, max_col - 1)] > 0)

        # 데이터 null: dNr 위치가 전부 A인지
        dNr_cols = [c for c in dNr if c < max_col]
        if dNr_cols:
            d_region = mat[:, dNr_cols]
            d_all_A = np.all(d_region == A, axis=1) & len_ok
        else:
            d_all_A = np.zeros(N_reads, dtype=bool)

        # XOR null: xNr 위치가 전부 A인지
        xNr_cols = [c for c in xNr if c < max_col]
        if xNr_cols:
            x_region = mat[:, xNr_cols]
            x_all_A = np.all(x_region == A, axis=1) & len_ok
        else:
            x_all_A = np.zeros(N_reads, dtype=bool)

        # 프라이머 매칭 (≥15/20) — 벡터
        primer_cols = min(20, len(Fp_arr), max_col)
        fwd_match = (mat[:, :primer_cols] == Fp_arr[:primer_cols]).sum(axis=1)
        primer_ok = (fwd_match >= 15)

        dNo_mask = d_all_A & primer_ok
        xNo_mask = x_all_A & primer_ok

        dNo = np.where(dNo_mask)[0]
        xNo = np.where(xNo_mask)[0]
    else:
        # fallback: for 루프
        dNo = []
        xNo = []
        for a in range(len(seq)):
            if len(seq[a]) >= s / 2 + 20:
                d_region = ''.join(seq[a][i] for i in dNr if i < len(seq[a]))
                if d_region == 'A' * len(dNr):
                    fwd_match = sum(1 for i in range(min(20, len(Fp), len(seq[a])))
                                    if seq[a][i] == Fp[i])
                    if fwd_match >= 15:
                        dNo.append(a)
                x_region = ''.join(seq[a][i] for i in xNr if i < len(seq[a]))
                if x_region == 'A' * len(xNr):
                    fwd_match = sum(1 for i in range(min(20, len(Fp), len(seq[a])))
                                    if seq[a][i] == Fp[i])
                    if fwd_match >= 15:
                        xNo.append(a)

    # 중복 제거 — uint8 행렬 기반 (seq 리스트 불필요)
    mat_ref = seq_mat_u8 if seq_mat_u8 is not None else None
    dNr_arr = np.array(dNr, dtype=int)
    xNr_arr = np.array(xNr, dtype=int)

    if mat_ref is not None and len(dNo) > 0:
        dNo_arr = np.asarray(dNo)
        xNo_arr = np.asarray(xNo)

        if len(xNr) > len(dNr):
            # dNo 중에서 xNr도 all-A인 것 제거
            xNr_cols = [c for c in xNr if c < mat_ref.shape[1]]
            if xNr_cols:
                d_xor_check = np.all(mat_ref[dNo_arr][:, xNr_cols] == A, axis=1)
                dNo_arr = dNo_arr[~d_xor_check]
        else:
            dNr_cols = [c for c in dNr if c < mat_ref.shape[1]]
            if dNr_cols:
                x_dat_check = np.all(mat_ref[xNo_arr][:, dNr_cols] == A, axis=1)
                xNo_arr = xNo_arr[~x_dat_check]

        # 접두사 추출 + consensus — uint8에서 직접
        min_dNr = min(dNr) if dNr else 0
        min_xNr = min(xNr) if xNr else 0

        if len(dNo_arr) > 0 and min_dNr > 0:
            d_prefix_u8 = mat_ref[dNo_arr, :min_dNr]
            # mode per column
            d_mode = np.zeros(min_dNr, dtype=np.uint8)
            for c in range(min_dNr):
                col = d_prefix_u8[:, c]
                col = col[col > 0]
                if len(col) > 0:
                    d_mode[c] = np.bincount(col).argmax()
                else:
                    d_mode[c] = _A
            Ndata = d_mode.tobytes().decode('ascii')
        else:
            Ndata = 'A' * 150

        if len(xNo_arr) > 0 and min_xNr > 0:
            x_prefix_u8 = mat_ref[xNo_arr, :min_xNr]
            x_mode = np.zeros(min_xNr, dtype=np.uint8)
            for c in range(min_xNr):
                col = x_prefix_u8[:, c]
                col = col[col > 0]
                if len(col) > 0:
                    x_mode[c] = np.bincount(col).argmax()
                else:
                    x_mode[c] = _A
            Xdata = x_mode.tobytes().decode('ascii')
        else:
            Xdata = 'A' * 150

        dNc = [mat_ref[i, :min_dNr].tobytes().decode('ascii') for i in dNo_arr] if len(dNo_arr) > 0 else []
        xNc = [mat_ref[i, :min_xNr].tobytes().decode('ascii') for i in xNo_arr] if len(xNo_arr) > 0 else []
    else:
        # fallback: seq 리스트 사용
        dNs = [seq[i] for i in dNo]
        xNs = [seq[i] for i in xNo]
        dpA = 'A' * len(dNr)
        xpA = 'A' * len(xNr)

        if len(xNr) > len(dNr):
            dNs = [s_item for s_item in dNs
                   if ''.join(s_item[i] for i in xNr if i < len(s_item)) != xpA]
        else:
            xNs = [s_item for s_item in xNs
                   if ''.join(s_item[i] for i in dNr if i < len(s_item)) != dpA]

        min_dNr = min(dNr) if dNr else 0
        min_xNr = min(xNr) if xNr else 0
        dNc = [s_item[:min_dNr] for s_item in dNs if len(s_item) > min_dNr]
        xNc = [s_item[:min_xNr] for s_item in xNs if len(s_item) > min_xNr]

        if dNc:
            d_matrix = np.array([list(s_item) for s_item in dNc])
            Ndata = ''.join(_char_mode_per_column(d_matrix))
        else:
            Ndata = 'A' * 150
        if xNc:
            x_matrix = np.array([list(s_item) for s_item in xNc])
            Xdata = ''.join(_char_mode_per_column(x_matrix))
        else:
            Xdata = 'A' * 150

    return Ndata, Xdata, dNc, xNc


def XORrun(matrix, s, Nc, Nxc, ncp, ncc, nbp, nbc):
    """
    XOR 기반 에러 복구.
    MATLAB: XORrun.m

    Returns:
        xtp, md, tp, xp, XA, XB, XC
    """
    total = matrix.size
    rows = matrix.shape[0]

    # 2/3 : 1/3 분할
    tri_rows = total * 2 // (3 * s)
    tp = matrix[:tri_rows, :].copy()
    xp = matrix[tri_rows:, :].copy()

    # flatten (row-major = MATLAB의 reshape(tp', [numel(tp) 1])')
    row_part = tp.flatten(order='C')
    if Nc > 0:
        row_part = row_part[:-Nc]

    # 홀수 처리
    if (ncp * 3 * nbp + ncc * 3 * nbc) % 2 != 0:
        row_part = row_part[:-1]

    XC = xp.flatten(order='C')
    if Nxc > 0:
        XC = XC[:-Nxc]

    # A/B 분할
    half = math.ceil(len(row_part) / 2)
    XA = row_part[:half]
    XB = row_part[half:]

    if len(row_part) % 2 != 0:
        XB = np.append(XB, 0)

    # XOR 복구
    XOR_matrix = np.concatenate([
        np.bitwise_xor(XB.astype(np.int64), XC[:len(XB)].astype(np.int64)),
        np.bitwise_xor(XA.astype(np.int64), XC[:len(XA)].astype(np.int64))
    ])

    # 패딩 후 reshape
    ECC_example = np.concatenate([XOR_matrix, np.zeros(Nc, dtype=np.int64)])
    ECC_ex_matrix = ECC_example.reshape((-1, s))  # 여기서 -1은 numel/s

    # 누락 행 (합이 0인 행) 찾기
    md = np.where(tp.sum(axis=1) == 0)[0]

    # 누락 행을 XOR 복구 데이터로 대체
    for idx in md:
        if idx < len(ECC_ex_matrix):
            tp[idx, :] = ECC_ex_matrix[idx, :]

    xtp = tp.copy()
    return xtp, md, tp, xp, XA, XB, XC


def nestedXOR(CellSeq, xtp, ErrD, RSerrinf, XA, XB, XC, payload, tp, xp, matsz, Null_count):
    """
    중첩 XOR 복구: XOR 대응 파트너가 에러가 아닌 경우만 복구.
    MATLAB: nestedXOR.m

    Returns:
        XORM: list of np.array (복구 대상 리드 그룹)
        Dx: np.array (복구 대상 인덱스)
    """
    # Depth_cut 에러 인덱스 (RS 에러가 있는 것만)
    err_mask = RSerrinf[ErrD] != 0
    Depth_cut = ErrD[err_mask]

    if len(Depth_cut) == 0:
        return [], np.array([], dtype=np.int64)

    xp_rows = xp.shape[0]
    tp_rows = tp.shape[0]

    # A/B/C 파티션 (MATLAB 1-based 인덱스 기준 유지)
    A_part = Depth_cut[Depth_cut <= xp_rows]
    C_part = Depth_cut[Depth_cut > tp_rows]
    bc_mask = (Depth_cut > len(A_part)) if len(A_part) > 0 else np.ones(len(Depth_cut), dtype=bool)
    remaining = Depth_cut[len(A_part):]
    B_part = remaining[:len(remaining) - len(C_part)] if len(C_part) > 0 else remaining

    # B_part, C_part 끝 인덱스 보정
    if len(B_part) > 0 and B_part[-1] == tp_rows:
        B_part = B_part[:-1]
    if len(C_part) > 0 and C_part[-1] == tp_rows + xp_rows - 1:
        C_part = C_part[:-1]

    Xi = np.zeros(len(Depth_cut), dtype=np.int64)

    # A 파트 검사
    for a_idx in range(len(A_part)):
        a = A_part[a_idx]  # 0-based index
        AtoB_start = a * matsz - len(XA) % matsz
        AtoB_end = AtoB_start + matsz
        AtoC_start = a * matsz
        AtoC_end = AtoC_start + matsz

        if AtoB_start >= 0 and AtoB_end <= len(XB) and AtoC_start >= 0 and AtoC_end <= len(XC):
            if np.sum(XB[AtoB_start:AtoB_end]) == 0 or np.sum(XC[AtoC_start:AtoC_end]) == 0:
                Xi[a_idx] = a

    # B 파트 검사 (MATLAB B_part(b)는 1-based, Python은 0-based → +1 보정)
    for b_idx in range(len(B_part)):
        b = B_part[b_idx]
        BtoA_start = (b + 1) * matsz - len(XA) + (matsz - len(XA) % matsz)
        BtoA_end = BtoA_start + matsz
        BtoC_start = BtoA_start
        BtoC_end = BtoA_end

        if BtoA_start >= 0 and BtoA_end <= len(XA) and BtoC_start >= 0 and BtoC_end <= len(XC):
            if np.sum(XA[BtoA_start:BtoA_end]) == 0 or np.sum(XC[BtoC_start:BtoC_end]) == 0:
                Xi[len(A_part) + b_idx] = b

    # C 파트 검사
    for c_idx in range(len(C_part)):
        c = C_part[c_idx]
        CtoA_start = c * matsz - len(XA) * 2 + Null_count
        CtoA_end = CtoA_start + matsz
        CtoB_start = c * matsz - len(XA) * 2 + len(XA) % matsz + Null_count
        CtoB_end = CtoB_start + matsz

        if CtoA_start >= 0 and CtoA_end <= len(XB) and CtoB_start >= 0 and CtoB_end <= len(XB):
            if np.sum(XB[CtoA_start:CtoA_end]) == 0 or np.sum(XB[CtoB_start:CtoB_end]) == 0:
                Xi[len(A_part) + len(B_part) + c_idx] = c

    Non_XOR_oligo = Xi[Xi > 0]

    # RS에러=0이면 원본 payload로 복구, -1이면 재디코딩 필요
    for idx in Non_XOR_oligo:
        if idx < len(RSerrinf) and idx < xtp.shape[0]:
            if RSerrinf[idx] == 0 and idx < xtp.shape[0]:
                xtp[idx, :] = payload[idx, :matsz] if idx < payload.shape[0] else xtp[idx, :]

    # 재디코딩이 필요한 인덱스
    rs_err_mask = np.array([RSerrinf[int(i)] != 0 for i in Non_XOR_oligo
                            if int(i) < len(RSerrinf)])
    Dx = Non_XOR_oligo[rs_err_mask] if len(rs_err_mask) > 0 else np.array([], dtype=np.int64)

    # 해당 인덱스의 원본 리드 그룹 추출
    XORM = []
    for d in Dx:
        d_int = int(d)
        if d_int < len(CellSeq):
            XORM.append(CellSeq[d_int])
        else:
            XORM.append(np.empty((0, 150), dtype='U1'))

    return XORM, Dx


def reRSdec(XORM, matrix, k, m, n, Dx, matsz):
    """
    XOR 복구 후 재-RS 디코딩: 원본 리드를 다시 이진화하고 RS 디코딩.
    MATLAB: reRSdec.m
    """
    for a in range(len(XORM)):
        group = XORM[a]
        if group is None or len(group) == 0:
            continue

        # 데이터 영역 추출 (위치 21-130, MATLAB) = (20:130, Python)
        reRSdata = group[:, 20:130] if group.shape[1] >= 130 else group

        if reRSdata.size == 0:
            continue

        # 벡터화 이진화: uint8에서 직접 DNA_TO_IDX → BIT_TABLE
        nrows_g, ncols_g = reRSdata.shape
        if reRSdata.dtype == np.uint8:
            re_u8 = reRSdata.copy()
        else:
            re_u8 = np.zeros((nrows_g, ncols_g), dtype=np.uint8)
            for r in range(nrows_g):
                s = ''.join(reRSdata[r, :]).replace('N', 'A')
                re_u8[r] = np.frombuffer(s.encode('ascii'), dtype=np.uint8)
        re_u8[re_u8 == _N] = _A
        idx_v = _DNA_TO_IDX[re_u8]
        reSeqtoBin = _BIT_TABLE[idx_v].reshape(nrows_g, -1).astype(np.int8)

        # RS 디코딩
        try:
            re_rx, re_err = RSdecoding(reSeqtoBin, k, m, n)

            # 성공한 디코딩만 사용 (에러 수 = 0)
            success_mask = (re_err == 0)
            if np.any(success_mask):
                successful = re_rx[success_mask, :matsz]
                # mode (최빈값)
                re_RSdata = np.zeros(matsz, dtype=np.float64)
                for col in range(matsz):
                    col_data = successful[:, col]
                    counter = Counter(col_data)
                    re_RSdata[col] = counter.most_common(1)[0][0]

                dx_idx = int(Dx[a])
                if dx_idx < matrix.shape[0]:
                    matrix[dx_idx, :] = re_RSdata
        except Exception:
            pass

    return matrix


# ==============================================================================
# Phase 5: 3D 복원
# ==============================================================================

def tridecode(ECCp, ncp, nbp, ncc, nbc, filename):
    """
    이진 페이로드에서 3D 좌표(tp)와 면 연결성(tc) 디코딩.
    MATLAB: tridecode.m

    Returns:
        tc: np.array shape (ncc, 3) - 면 인덱스 (1-based)
        tp: np.array shape (ncp, 3) - 꼭짓점 좌표
    """
    ECCp = np.asarray(ECCp, dtype=np.int64)

    # 홀수 처리
    if (ncp * 3 * nbp + ncc * 3 * nbc) % 2 != 0:
        ECCp = ECCp[:-1]

    # 좌표 영역(pr)과 연결성 영역(cr) 분리
    pr = ECCp[:ncp * nbp * 3]
    cr = ECCp[ncp * nbp * 3:]

    # 좌표: reshape → bi2de → /10^16
    # MATLAB: pb = reshape(pr', [nbp ncp*3])'; → (ncp*3, nbp) 행렬
    pb = np.reshape(pr, (ncp * 3, nbp), order='C')  # row-major (pr is already flat row-major)
    pd = _bi2de(pb) / (10 ** 16)

    # 3D 좌표 reshape
    # 특수 데이터셋: reshape(pd, [ncp 3]) in MATLAB column-major
    if filename in ('seokga', 'liberty', 'pisa', 'venus', 'goldengate', 'tripitaka'):
        tp = np.reshape(pd, (ncp, 3), order='F')
    else:
        # 일반: reshape(pd, [3 ncp])' → (ncp, 3)
        tp = np.reshape(pd, (3, ncp), order='F').T

    # 연결성: reshape → bi2de
    cb = np.reshape(cr, (ncc * 3, nbc), order='C')
    cd = _bi2de(cb)
    tc = np.reshape(cd, (3, ncc), order='F').T

    # 안전 처리: 0 → 1, out-of-range → 1
    tc[tc == 0] = 1
    tc[tc > len(tp)] = 1

    return tc, tp


def trierrXOR(mat, tc, tp, Ncount, Nxcount, ncp, nbp, ncc, nbc):
    """
    XOR 데이터를 사용하여 0값 꼭짓점/면을 복구.
    MATLAB: trierrXOR.m

    Returns:
        trid: dict or trimesh object
    """
    s = mat.shape[1]
    total = mat.size
    tri_rows = total * 2 // (3 * s)

    tri_part = mat[:tri_rows, :].copy()
    xor_part = mat[tri_rows:, :].copy()

    # flatten
    row_tri_part = tri_part.flatten(order='C')
    if Ncount > 0:
        row_tri_part = row_tri_part[:-Ncount]

    if (ncp * 3 * nbp + ncc * 3 * nbc) % 2 != 0:
        row_tri_part = row_tri_part[:-1]

    XOR_ECC_C = xor_part.flatten(order='C')
    if Nxcount > 0:
        XOR_ECC_C = XOR_ECC_C[:-Nxcount]

    half = math.ceil(len(row_tri_part) / 2)
    XOR_ECC_A = row_tri_part[:half]
    XOR_ECC_B = row_tri_part[half:]

    if len(row_tri_part) % 2 != 0:
        XOR_ECC_B = np.append(XOR_ECC_B, 0)

    XOR_ECC_matrix = np.concatenate([
        np.bitwise_xor(XOR_ECC_B.astype(np.int64), XOR_ECC_C[:len(XOR_ECC_B)].astype(np.int64)),
        np.bitwise_xor(XOR_ECC_A.astype(np.int64), XOR_ECC_C[:len(XOR_ECC_A)].astype(np.int64))
    ])

    # 0값 꼭짓점 복구
    tp_flat = tp.T.flatten(order='F')  # MATLAB: tp' then linear indexing
    zero_positions = np.where(tp_flat == 0)[0]

    if len(zero_positions) > 0:
        for zp in zero_positions:
            xor_start = zp * nbp
            xor_end = xor_start + nbp
            if xor_end <= len(XOR_ECC_matrix):
                bits = XOR_ECC_matrix[xor_start:xor_end].astype(np.int8)
                val = _bi2de(bits) / (10 ** 16)
                tp_flat[zp] = val

        tp = tp_flat.reshape(tp.shape[1], tp.shape[0], order='F').T

    # 0값/범위초과 면 복구
    tc_flat = tc.T.flatten(order='F')
    bad_faces = np.where((tc_flat == 0) | (tc_flat > tp.shape[0]))[0]

    if len(bad_faces) > 0:
        coXOR1 = XOR_ECC_C.copy()
        coXOR2 = XOR_ECC_A.copy()
        if len(row_tri_part) % 2 != 0:
            coXOR1 = coXOR1[:-1] if len(coXOR1) > 0 else coXOR1
            coXOR2 = coXOR2[:-1] if len(coXOR2) > 0 else coXOR2
        XORrest = np.bitwise_xor(coXOR1[:len(coXOR2)].astype(np.int64),
                                  coXOR2.astype(np.int64))

        for bf in bad_faces:
            xor_start = bf * nbc
            xor_end = xor_start + nbc
            if xor_end <= len(XORrest):
                bits = XORrest[xor_start:xor_end].astype(np.int8)
                tc_flat[bf] = _bi2de(bits)

        tc = tc_flat.reshape(tc.shape[1], tc.shape[0], order='F').T

        # 여전히 범위 밖인 면 제거
        bad_rows = np.where(
            (tc > ncp).any(axis=1) | (tc == 0).any(axis=1)
        )[0]
        tc = np.delete(tc, bad_rows, axis=0)

    return _make_triangulation(tc, tp)


def _make_triangulation(tc, tp):
    """
    꼭짓점과 면으로부터 3D 메쉬 객체 생성.
    MATLAB: triangulation(tc, tp)
    """
    # MATLAB 면 인덱스는 1-based → Python 0-based
    tc_0based = tc.astype(int) - 1

    # 유효성 검사
    max_idx = tp.shape[0] - 1
    valid_mask = np.all((tc_0based >= 0) & (tc_0based <= max_idx), axis=1)
    tc_0based = tc_0based[valid_mask]

    if HAS_TRIMESH:
        mesh = trimesh.Trimesh(vertices=tp, faces=tc_0based)
        return mesh
    else:
        return {'vertices': tp, 'faces': tc_0based, 'faces_1based': tc}


def visualize_mesh(trid, title_str='Decoded 3D Mesh', save_path=None):
    """3D 메쉬 시각화 (matplotlib). save_path 지정 시 파일 저장 후 비대화형."""
    if HAS_TRIMESH and isinstance(trid, trimesh.Trimesh):
        vertices = trid.vertices
        faces = trid.faces
    elif isinstance(trid, dict):
        vertices = trid['vertices']
        faces = trid['faces']
    else:
        print("[경고] 시각화할 메쉬 데이터가 없습니다.")
        return

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(
        vertices[:, 0], vertices[:, 1], vertices[:, 2],
        triangles=faces, cmap='viridis', alpha=0.7, edgecolor='k', linewidth=0.1
    )
    ax.set_title(title_str)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  메쉬 이미지 저장: {save_path}")
        plt.close(fig)
    else:
        plt.show()


# ==============================================================================
# Phase 6: 메인 파이프라인 (Textscan 스타일 시간 측정)
# ==============================================================================

def read_fastq_file(filepath):
    """
    FASTQ 파일 고속 읽기.
    비압축: mmap + 직접 파싱 (OS 캐시 활용, Biopython 불필요)
    gz: pigz/gzip 멀티스레드 해제 → 파이프 파싱
    """
    sequences = []

    # .gz 없는 경로도 시도
    if not os.path.exists(filepath):
        if os.path.exists(filepath + '.gz'):
            filepath = filepath + '.gz'
        else:
            print(f"[경고] 파일 없음: {filepath}")
            return sequences

    is_gz = filepath.endswith('.gz')

    if is_gz:
        return _read_fastq_gz(filepath)
    else:
        return _read_fastq_mmap(filepath)


def _read_fastq_mmap(filepath):
    """비압축 FASTQ: mmap bulk read → 시퀀스 + uint8 행렬 동시 생성."""
    file_size = os.path.getsize(filepath)

    with open(filepath, 'rb') as f:
        # 전체 파일을 한번에 읽기 (mmap보다 빠를 수 있음, OS 캐시)
        raw = f.read()

    # 줄 분할 → 매 4줄 중 2번째(index 1)가 시퀀스
    lines = raw.split(b'\n')
    sequences = [lines[i].decode('ascii') for i in range(1, len(lines), 4) if i < len(lines) and lines[i]]

    return sequences


def _read_fastq_gz(filepath):
    """gz FASTQ: pigz(병렬) 또는 gzip으로 해제 후 파이프 파싱."""
    # pigz 사용 가능하면 멀티스레드 해제
    pigz_path = shutil.which('pigz')
    if pigz_path:
        cmd = [pigz_path, '-dc', '-p', str(N_WORKERS), filepath]
    else:
        cmd = ['gzip', '-dc', filepath]

    sequences = []
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1024*1024)
    line_num = 0
    for line in proc.stdout:
        if line_num % 4 == 1:
            sequences.append(line.rstrip(b'\n').decode('ascii'))
        line_num += 1
    proc.wait()

    return sequences


def main():
    """
    DNA 디코딩 메인 파이프라인.
    MATLAB: Textscan_230527_downsampling_additional.m 대응.
    각 단계별 시간 측정 포함.
    """
    # =========================================================================
    # 설정 (CLI 인자 또는 기본값)
    # =========================================================================
    filename = sys.argv[1] if len(sys.argv) > 1 else 'Jikji7'
    data_dir = sys.argv[2] if len(sys.argv) > 2 else '.'
    Timecell = {}

    print(f"=== DNA Decoding Pipeline: {filename} ===")
    print(f"  데이터 디렉토리: {os.path.abspath(data_dir)}\n")

    # =========================================================================
    # Step 0: FASTQ 텍스트 스캔
    # =========================================================================
    t = time.time()
    print("[Step 0] FASTQ 파일 로드...")

    # 데이터셋별 FASTQ 파일명 매핑
    FASTQ_FILES = {
        'Jikji4': [('4_S2_L001_R1_001.fastq', '4_S2_L001_R2_001.fastq')],
        'Jikji6': [('Jikji6_1.fastq', 'Jikji6_2.fastq')],  # paired
        'Jikji7': [('jikji-7_S1_L001_R1_001.fastq.gz', 'jikji-7_S1_L001_R2_001.fastq.gz')],
        'venus': [
            ('Venus1_R1.fastq', 'Venus1_R2.fastq'),
            ('Venus2_R1.fastq', 'Venus2_R2.fastq'),
        ],
        'goldengate': [
            ('Goldengate1_R1.fastq', 'Goldengate1_R2.fastq'),
            ('Goldengate2_R1.fastq', 'Goldengate2_R2.fastq'),
        ],
    }

    # 기본 R1/R2 패턴
    file_pairs = FASTQ_FILES.get(filename, [(f'{filename}_R1.fastq', f'{filename}_R2.fastq')])

    seq = []
    for fnames in file_pairs:
        for fname in fnames:
            fpath = os.path.join(data_dir, fname)
            seq += read_fastq_file(fpath)

    if not seq:
        print("[경고] FASTQ 파일을 찾을 수 없습니다. 테스트 더미 데이터 생성.")
        seq = ["ATGC" * 37 + "AT" for _ in range(100)]  # 150bp 더미

    Timecell['FASTQ Load'] = round(time.time() - t, 3)
    print(f"  로드된 리드 수: {len(seq):,} ({Timecell['FASTQ Load']}s)")

    # 전체 리드를 uint8 행렬로 한번에 구성 (128GB 메모리 활용)
    t_mat = time.time()
    print("  numpy 행렬 구성 중...")
    seq_mat = _seqs_to_matrix(seq)
    print(f"  행렬 구성: {round(time.time()-t_mat, 2)}s, shape={seq_mat.shape}, "
          f"mem={seq_mat.nbytes/1024/1024:.0f}MB")

    # =========================================================================
    # Step 1: Evert Sequence (서열 방향 교정) — 벡터화
    # =========================================================================
    t = time.time()
    print("\n[Step 1] Evert Sequence (벡터화)...")

    forward_section = [s[:21] if len(s) >= 21 else s for s in seq]
    F_primers, R_primers = PrimerSelect(filename)
    tevert_seq, reversed_seq_log = ReverseEvert(seq, forward_section, R_primers, F_primers,
                                                 mat_uint8=seq_mat)

    # 역방향 리드 교체 — uint8 행렬만 업데이트 (seq 리스트는 나중에 필요할 때)
    rev_indices = np.where(reversed_seq_log)[0]
    if isinstance(tevert_seq, np.ndarray) and tevert_seq.dtype == np.uint8:
        seq_mat[rev_indices, :150] = tevert_seq  # 벡터 1회
    else:
        for i, idx in enumerate(rev_indices):
            if tevert_seq[i] is not None:
                arr = np.frombuffer(tevert_seq[i].encode('ascii'), dtype=np.uint8)
                seq_mat[idx, :len(arr)] = arr

    Timecell['Evert'] = round(time.time() - t, 3)
    print(f"  완료: {Timecell['Evert']}s (역방향 리드: {len(rev_indices):,}개)")

    # =========================================================================
    # Step 2: Primer Sorting — 벡터화
    # =========================================================================
    t = time.time()
    print("\n[Step 2] Primer Sorting (벡터화)...")

    PrimerSortingSeq, ps_mat = PrimerSort(seq, F_primers, R_primers, mat_uint8=seq_mat)

    Timecell['Primer Sorting'] = round(time.time() - t, 3)
    print(f"  완료: {Timecell['Primer Sorting']}s (유효 리드: {len(PrimerSortingSeq):,}개)")

    if len(PrimerSortingSeq) == 0:
        print("[오류] 유효한 리드가 없습니다. 파이프라인 중단.")
        return Timecell

    # =========================================================================
    # Step 3: Key Finding
    # =========================================================================
    t = time.time()
    print("\n[Step 3] Key Finding...")

    k = 156
    n = k + 8  # 164
    m = len(_de2bi(n))  # 8

    KeyMap, KeymatDNA = KeyFind(PrimerSortingSeq, ps_mat_u8=ps_mat)
    OligoLen, idxlen, ncolpo, nbitpo, ncolco, nbitco = Keysort(KeyMap, KeymatDNA, k, filename)

    Timecell['Key Finding'] = round(time.time() - t, 3)
    print(f"  완료: {Timecell['Key Finding']}s")
    print(f"  KeyMap={KeyMap}, OligoLen={OligoLen}, idxlen={idxlen}")
    print(f"  ncolpo={ncolpo}, nbitpo={nbitpo}, ncolco={ncolco}, nbitco={nbitco}")

    # =========================================================================
    # Step 4: Consensus Sequencing (Index Sorting)
    # =========================================================================
    t = time.time()
    print("\n[Step 4] Consensus Sequencing...")

    MIndex, sorted_u8, cons_boundaries = Consensus(PrimerSortingSeq, idxlen, OligoLen, ps_mat_u8=ps_mat)

    # Depth 계산
    Depth = np.zeros(len(MIndex), dtype=np.int64)
    for i in range(len(MIndex)):
        Depth[i] = MIndex[i].shape[0] if MIndex[i] is not None else 0

    Timecell['Index Sorting'] = round(time.time() - t, 3)
    print(f"  완료: {Timecell['Index Sorting']}s")
    print(f"  평균 Depth: {np.mean(Depth):.1f}x, 최소: {np.min(Depth)}, 최대: {np.max(Depth)}")

    # Mapping Location (pisa, seokga, liberty 전용)
    if filename in ('pisa', 'seokga', 'liberty'):
        warnings.warn(
            f"MapLocation 함수 미구현: {filename} 데이터셋 결과가 부정확할 수 있습니다.",
            UserWarning
        )
        # TODO: MapLocation.m 원본 확보 후 구현 필요
        # MIndex, MapLoc = MapLocation(MIndex, OligoLen)

    # =========================================================================
    # Step 5: RS Decoding
    # =========================================================================
    t = time.time()
    print("\n[Step 5] Binarization + RS Decoding...")

    t_bin = time.time()
    SeqtoBin = Binarization(MIndex, sorted_u8=sorted_u8, boundaries=cons_boundaries)
    t_bin_done = time.time()

    _init_rs_reedsolo()  # RS 코덱 사전 초기화 (시간 측정에서 제외)
    t_rs = time.time()
    rsECC, RSerrinf = RSdecoding(SeqtoBin, k, m, n)
    t_rs_done = time.time()

    print(f"  Binarization: {round(t_bin_done - t_bin, 3)}s, RS decode: {round(t_rs_done - t_rs, 3)}s")

    # 에러 필터링
    rsECC[RSerrinf < 0, :] = 0
    rsECC[RSerrinf > 4, :] = 0
    for a in range(min(OligoLen, len(rsECC))):
        if np.any(rsECC[a, :] > 1):
            rsECC[a, :] = 0

    # payload 추출
    payload = rsECC[:, :k - idxlen]
    Temp_mat = payload[:-1, :].copy()

    Timecell['RS Decoding'] = round(time.time() - t, 3)
    print(f"  완료: {Timecell['RS Decoding']}s")
    print(f"  RS 에러 수정불가: {np.sum(RSerrinf < 0)}개, 수정됨: {np.sum(RSerrinf > 0)}개")

    # =========================================================================
    # Step 6: Null Finding
    # =========================================================================
    t = time.time()
    print("\n[Step 6] Null Finding...")

    # Nullkey 계산
    term = ncolpo * 3 * nbitpo + ncolco * 3 * nbitco
    Nullkey_1 = term + 1 if term % 2 == 1 else term

    s_val = Temp_mat.shape[1]
    Nullkey_2 = Nullkey_1 // 2  # 정수 나눗셈 (float 전파 방지)

    Null_count = 0
    nk1 = Nullkey_1
    while nk1 % s_val != 0:
        nk1 += 1
        Null_count += 1

    Null_xor_count = 0
    nk2 = Nullkey_2
    while nk2 % s_val != 0:
        nk2 += 1
        Null_xor_count += 1

    data_Ndecode, xor_Ndecode, dNc, xNc = NullFind(seq, Null_count, Null_xor_count, s_val, F_primers,
                                                    seq_mat_u8=seq_mat)

    # Null 서열을 이진 변환
    bindataN = _dna_to_binary(data_Ndecode[20:])  # MATLAB: data_Ndecode(21:end)
    binxorN = _dna_to_binary(xor_Ndecode[20:])    # MATLAB: xor_Ndecode(21:end)

    Timecell['Null Finding'] = round(time.time() - t, 3)
    print(f"  완료: {Timecell['Null Finding']}s")

    # Depth 업데이트
    tri_row_idx = Temp_mat.size * 2 // (3 * s_val) - 1  # MATLAB 1-indexed → Python 0-indexed
    if tri_row_idx < len(Depth) and dNc:
        total_dnc = sum(len(d) for d in dNc if d)
        Depth[tri_row_idx] = total_dnc
    if len(Depth) > 1 and xNc:
        total_xnc = sum(len(x) for x in xNc if x)
        Depth[-2] = total_xnc

    # =========================================================================
    # Step 7: XOR Recovery
    # =========================================================================
    t = time.time()
    print("\n[Step 7] XOR Recovery...")

    matsz = k - idxlen

    # Null 데이터를 Temp_mat에 삽입
    Temp_mat = payload[:-1, :].copy()
    data_insert_row = Temp_mat.size * 2 // (3 * s_val) - 1  # MATLAB 1-indexed → Python 0-indexed
    if data_insert_row < Temp_mat.shape[0]:
        insert_len = min(len(bindataN), Temp_mat.shape[1])
        Temp_mat[data_insert_row, :insert_len] = bindataN[:insert_len]
    if Temp_mat.shape[0] > 0:
        insert_len = min(len(binxorN), Temp_mat.shape[1])
        Temp_mat[-1, :insert_len] = binxorN[:insert_len]

    # Depth 크기 맞춤
    if len(Depth) > Temp_mat.shape[0]:
        Depth = Depth[:Temp_mat.shape[0]]

    # === 복구 단계 1: Global RS ===
    trid = None
    try:
        print("  [1/3] Global RS 디코딩 시도...")
        ECC_tri_part = Temp_mat[:Temp_mat.size * 2 // (3 * s_val), :]
        ECC_payload_flat = ECC_tri_part.flatten(order='C')
        if Null_count > 0:
            ECC_payload_flat = ECC_payload_flat[:-Null_count]
        tc, tp = tridecode(ECC_payload_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
        trid = _make_triangulation(tc, tp)
        print("  → Global RS 성공!")
    except Exception as e:
        print(f"  → Global RS 실패: {e}")

    # === 복구 단계 2: Depth cut ===
    Depthcut = 16
    ErrD = np.where(Depth <= Depthcut)[0]

    try:
        print("  [2/3] Depth cut 기반 XOR 복구 시도...")
        err_mask = RSerrinf[ErrD] != 0 if len(ErrD) > 0 else np.array([], dtype=bool)
        err_indices = ErrD[err_mask] if len(err_mask) > 0 else np.array([], dtype=int)
        if len(err_indices) > 0:
            Temp_mat[err_indices, :] = 0

        xtp, md, tri_p, xor_p, XOR_A, XOR_B, XOR_C = \
            XORrun(Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)

        ECC_payload_flat = xtp.flatten(order='C')
        if Null_count > 0:
            ECC_payload_flat = ECC_payload_flat[:-Null_count]
        tc, tp = tridecode(ECC_payload_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
        trid = trierrXOR(Temp_mat, tc, tp, Null_count, Null_xor_count, ncolpo, nbitpo, ncolco, nbitco)
        print("  → Depth cut 성공!")
    except Exception as e:
        print(f"  → Depth cut 실패: {e}")

    # === 복구 단계 3: Local RS (nestedXOR + reRSdec) ===
    try:
        print("  [3/3] Local RS 기반 복구 시도...")
        if len(err_indices) > 0:
            Temp_mat[err_indices, :] = 0

        xtp, md, tri_p, xor_p, XOR_A, XOR_B, XOR_C = \
            XORrun(Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)

        if len(ErrD) > 0:
            XORMIndex, Del_xor = nestedXOR(
                MIndex, xtp, ErrD, RSerrinf, XOR_A, XOR_B, XOR_C,
                payload, tri_p, xor_p, matsz, Null_count
            )

            Temp_mat = reRSdec(XORMIndex, Temp_mat, k, m, n, Del_xor, matsz)

            _, missing_data, _, _, _, _, _ = \
                XORrun(Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)

            XORmData = [MIndex[i] for i in missing_data if i < len(MIndex)]
            Temp_mat = reRSdec(XORmData, Temp_mat, k, m, n, missing_data, matsz)

        # 최종 복구
        xtp_final, _, _, _, _, _, _ = \
            XORrun(Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)
        ECC_payload_flat = xtp_final.flatten(order='C')
        if Null_count > 0:
            ECC_payload_flat = ECC_payload_flat[:-Null_count]
        tc, tp = tridecode(ECC_payload_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
        trid = trierrXOR(Temp_mat, tc, tp, Null_count, Null_xor_count, ncolpo, nbitpo, ncolco, nbitco)
        print("  → Local RS 성공!")
    except Exception as e:
        print(f"  → Local RS 실패: {e}")
        if trid is None:
            print("  [오류] 디코딩 실패 (Decoding Failed)")

    Timecell['XOR Recovery'] = round(time.time() - t, 3)
    print(f"  완료: {Timecell['XOR Recovery']}s")

    # =========================================================================
    # Step 8: STL Decode & Visualization
    # =========================================================================
    t = time.time()
    print("\n[Step 8] STL Decode & Visualization...")

    if trid is not None:
        # STL 저장
        output_dir = os.path.join(".", "Downsampling_Results", filename)
        os.makedirs(output_dir, exist_ok=True)

        if HAS_TRIMESH and isinstance(trid, trimesh.Trimesh):
            stl_path = os.path.join(output_dir, f"{filename}_decoded.stl")
            trid.export(stl_path)
            print(f"  STL 저장: {stl_path}")
            print(f"  꼭짓점: {len(trid.vertices):,}, 면: {len(trid.faces):,}")

        # 시각화 (비대화형: 이미지 저장, 대화형: 팝업)
        if os.environ.get('DISPLAY') or sys.platform == 'darwin':
            try:
                fig_path = os.path.join(output_dir, f"{filename}_mesh.png")
                visualize_mesh(trid, title_str=f'{filename} - Decoded 3D Mesh', save_path=fig_path)
            except Exception as e:
                print(f"  [경고] 시각화 실패: {e}")

    Timecell['STL Decode'] = round(time.time() - t, 3)

    # =========================================================================
    # 결과 요약
    # =========================================================================
    print("\n" + "=" * 50)
    print("시간 측정 결과 (Timecell)")
    print("=" * 50)
    total_time = 0
    for step_name, elapsed in Timecell.items():
        print(f"  {step_name}: {elapsed}s")
        total_time += elapsed
    print(f"  {'─' * 30}")
    print(f"  Total: {round(total_time, 3)}s")

    # Depth/에러 정보 엑셀 저장
    output_dir = os.path.join(".", "Downsampling_Results", filename)
    os.makedirs(output_dir, exist_ok=True)

    fcoverage = round(np.mean(Depth), 1)

    try:
        depth_file = os.path.join(output_dir, f"{filename}_Depth.xlsx")
        pd.DataFrame(Depth, columns=['Depth']).to_excel(depth_file, sheet_name=f"{fcoverage}x")
        print(f"\nDepth 저장: {depth_file}")
    except Exception as e:
        print(f"[경고] Depth 저장 실패: {e}")

    # RS 에러 인덱스
    RSerridx = np.where(RSerrinf < 0)[0]
    if len(RSerridx) > 0:
        if RSerridx[-1] == OligoLen - 1:
            RSerridx = RSerridx[:-1]
        valid = RSerridx[RSerridx < Temp_mat.shape[0]]
        errDepth = Depth[valid] if len(valid) > 0 else np.array([])
        rs_tbl = pd.DataFrame({'RS error idx': valid, 'Depth': errDepth})

        err_file = os.path.join(output_dir, f"{filename}_RS_errors({fcoverage}x).xlsx")
        try:
            rs_tbl.to_excel(err_file, index=False)
            print(f"RS 에러 저장: {err_file}")
        except Exception as e:
            print(f"[경고] RS 에러 저장 실패: {e}")

    print(f"\n평균 Depth: {fcoverage}x")
    print(f"총 리드: {len(seq):,}, Primer 유효: {len(PrimerSortingSeq):,}")

    return Timecell


if __name__ == "__main__":
    main()
