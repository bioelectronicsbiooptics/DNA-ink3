#!/usr/bin/env python3
"""
================================================================================
Figure 1 — iSCAN DNA Data Storage Decoding Pipeline
================================================================================
Role    : Core decoding pipeline for the iSCAN DNA-data-storage protocol.
          Python port of the MATLAB original (Single_Decoding_Time.m).
          (Final Production Version, 2026-05-07)

Pipeline: FASTQ -> ReverseEvert -> PrimerSort -> KeyFind -> Consensus
          -> Binarization -> RSdecoding -> NullFind -> XOR -> tridecode
================================================================================
"""

import os
import sys
import time
import math
import json
import warnings
from collections import Counter

import numpy as np
import pandas as pd
# Force non-interactive backend BEFORE importing pyplot — otherwise plt.show()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Module-level dict populated by KeyFind() each invocation. Lets main()
# capture KIS sub-timings without changing KeyFind's return signature.
KEYFIND_LAST_TIMING = {}

# Biopython (FASTQ parsing)
try:
    from Bio import SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("[warn] Biopython not installed. 'pip install biopython'")

# galois (Reed-Solomon)
try:
    import galois
    HAS_GALOIS = True
except ImportError:
    HAS_GALOIS = False
    print("[warn] galois not installed. 'pip install galois'")

# trimesh (3D mesh / STL)
try:
    import trimesh
    HAS_TRIMESH = True
except ImportError:
    HAS_TRIMESH = False
    print("[warn] trimesh not installed. 'pip install trimesh'")

# MATLAB Engine fallback
try:
    import matlab.engine
    HAS_MATLAB = True
except ImportError:
    HAS_MATLAB = False

from scipy import stats


# ==============================================================================
# PHASE 1: UTILITY FUNCTIONS
# ==============================================================================

def _seq_reverse_complement(seq_str):
    """
    Computes the reverse complement of a DNA sequence (MATLAB seqrcomplement equivalent).

    Args:
        seq_str (str): Input DNA sequence.
    Returns:
        str: Reverse complemented sequence (A<->T, G<->C, N->N).
    """
    complement_map = str.maketrans('ATGCatgcN', 'TACGtacgN')
    return seq_str.translate(complement_map)[::-1]


def _dna_to_binary(dna_str):
    """
    Converts a DNA string into a binary bit array.
    Mapping: A=00, T=01, G=10, C=11 (consistent with iSCAN protocol).

    Args:
        dna_str (str): Input DNA sequence.
    Returns:
        np.array: Binary representation as an int8 array.
    """
    mapping = {'A': [0, 0], 'T': [0, 1], 'G': [1, 0], 'C': [1, 1], 'N': [0, 0]}
    bits = []
    for ch in dna_str:
        bits.extend(mapping.get(ch, [0, 0]))
    return np.array(bits, dtype=np.int8)


def _dna_matrix_to_binary(dna_matrix):
    """
    Converts a 2D DNA character matrix into a binary matrix.
    Each base corresponds to 2 bits; the output width is doubled.

    Args:
        dna_matrix (np.array): 2D character array (rows=reads, cols=bases).
    Returns:
        np.array: 2D binary matrix (int8).
    """
    nrows, ncols = dna_matrix.shape
    result = np.zeros((nrows, ncols * 2), dtype=np.int8)
    # character -> index mapping: A=0, T=1, G=2, C=3
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
    Converts a binary bit array to a decimal integer (LSB-first).
    Compatible with MATLAB bi2de() functionality.

    Args:
        bits (array-like): 1D or 2D binary bits.
    Returns:
        int/np.array: Decimal equivalent(s).
    """
    bits = np.asarray(bits, dtype=np.int64)
    if bits.ndim == 1:
        return int(np.dot(bits, 2 ** np.arange(len(bits))))
    else:
        # 2D: convert each row
        return (bits * (2 ** np.arange(bits.shape[1]))).sum(axis=1).astype(np.int64)


def _de2bi(val, num_bits=None):
    """
    Converts a decimal integer to a binary bit array (LSB-first).
    Compatible with MATLAB de2bi() functionality.

    Args:
        val (int): Decimal value.
        num_bits (int, optional): Required bit length.
    Returns:
        np.array: Binary representation.
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
    Computes the mode (most frequent character) per column of a matrix.
    On ties, the lexicographically minimum character is returned (MATLAB mode behavior).

    Args:
        char_matrix (np.array): 2D character matrix.
    Returns:
        np.array: 1D array containing the mode of each column.
    """
    nrows, ncols = char_matrix.shape
    result = np.empty(ncols, dtype=char_matrix.dtype)
    for c in range(ncols):
        col = char_matrix[:, c]
        counter = Counter(col)
        max_count = counter.most_common(1)[0][1]
        # MATLAB mode(): returns minimum value on ties -> lexicographic min
        candidates = [k for k, v in counter.items() if v == max_count]
        result[c] = min(candidates)
    return result


# ==============================================================================
# PHASE 2: INITIAL PIPELINE FUNCTIONS
# ==============================================================================

def PrimerSelect(filename):
    """
    Retrieves the Forward and Reverse primers for a specified dataset.

    Args:
        filename (str): Name of the dataset (e.g., 'Jikji1').
    Returns:
        tuple: (Forward Primer, Reverse Primer)
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
        raise ValueError(f"Unknown dataset: {filename}")
    return PRIMER_DB[filename]


def ReverseEvert(seq, forward_section, r_primers, f_primers):
    """
    Detect reverse-orientation reads and apply reverse-complement conversion.
    MATLAB: ReverseEvert.m (Single_Decoding_Time.m v240108 style).

    First pass: flag reads where F20 OR F21 matches RC_R_primer with >= 18/20.
    Second pass: for flagged reads, EXACT 20/20 match decides slicing window.
                 (Two independent if statements; if both match, second overrides.)

    Args:
        seq: list of str (DNA sequences)
        forward_section: list of str (first 21 bp of each read)
        r_primers: str (reverse primer)
        f_primers: str (forward primer; kept for signature compatibility)

    Returns:
        tevert_seq: list of str (reverse-complemented sequences; None for unfit reads)
        reversed_seq_log: np.array of bool (whether each read was flagged)
    """
    # reverse complement of the reverse primer
    CompP = _seq_reverse_complement(r_primers)

    # First pass: flag any read with >= 18/20 match in either F20 or F21
    n = len(seq)
    reversed_seq_log = np.zeros(n, dtype=bool)

    for a in range(n):
        fs = forward_section[a]
        if len(fs) < 21:
            continue
        F20 = fs[0:20]
        F21 = fs[1:21]
        match_20 = sum(1 for i in range(min(20, len(CompP))) if F20[i] == CompP[i])
        match_21 = sum(1 for i in range(min(20, len(CompP))) if F21[i] == CompP[i])
        if match_20 >= 18 or match_21 >= 18:
            reversed_seq_log[a] = True

    # Second pass: only EXACT (20/20) matches drive the slicing window.
    # Two independent if-blocks: case 2 (length >= 151) overrides case 1 if both match.
    eSeq_indices = np.where(reversed_seq_log)[0]
    tevert_seq = [None] * len(eSeq_indices)

    for idx, a in enumerate(eSeq_indices):
        s = seq[a]
        seq_len = len(s)

        # case 1: strcmp(s[0:20], CompP) AND length >= 150
        if seq_len >= 20 and seq_len >= 150 and s[0:20] == CompP:
            tevert_seq[idx] = _seq_reverse_complement(s[0:150])
        # case 2: strcmp(s[1:21], CompP) AND length >= 151 (NOT elif — overrides)
        if seq_len >= 21 and seq_len >= 151 and s[1:21] == CompP:
            tevert_seq[idx] = _seq_reverse_complement(s[1:151])

    return tevert_seq, reversed_seq_log


def PrimerSort(seq, f_primers, r_primers):
    """
    Filter valid reads by primer matching.
    """
    valid = []
    for s in seq:
        chosen = None
        if len(s) >= 150:
            fwd_match = sum(1 for i in range(20) if s[i] == f_primers[i])
            rev_match = sum(1 for i in range(20) if s[130 + i] == r_primers[i])
            if fwd_match >= 18 and rev_match >= 18:
                chosen = s[0:150]
        if len(s) >= 151:
            fwd_match2 = sum(1 for i in range(20) if s[1 + i] == f_primers[i])
            rev_match2 = sum(1 for i in range(20) if s[131 + i] == r_primers[i])
            if fwd_match2 >= 18 and rev_match2 >= 18:
                chosen = s[1:151]  # overrides branch 1 if both matched
        if chosen is not None:
            valid.append(chosen)
    return valid


def KeyFind(primer_sorting_seq):
    """
    Returns:
        Key: int (key value, bi2de result)
        KmatDNA: str (key sequence consensus)
    """
    # build character matrix
    dna_matrix = np.array([list(s) for s in primer_sorting_seq])
    # payload: positions 21-98 (MATLAB 1-indexed) = 20:98 (Python 0-indexed)
    payload = dna_matrix[:, 20:98]  # shape: (N, 78)

    mDNA = ('TGCA' * 27 + 'TG')  # reference pattern (110 chars)
    Ks_list = []
    used_threshold = None
    # v8_10 sub-timers for the manuscript's KIS detection claim.
    kis_mask_seconds = 0.0   # cumulative across threshold attempts
    ref_match_seconds = 0.0  # cumulative across threshold attempts

    # Adaptive threshold: try 3, then 2, then 1 (matches 2320 changelog).
    # For each threshold, re-evaluate KIS mask AND reference-pattern matching
    # together — only break out of the threshold loop when reference hits exist.
    for threshold in (3, 2, 1):
        t0 = time.time()
        key_mask = np.zeros(len(payload), dtype=bool)
        for a in range(payload.shape[0]):
            row_str = ''.join(payload[a, :])
            count_5A = row_str.count('AAAAA')
            count_4A = row_str.count('AAAA')
            count_GAAAT = row_str.count('GAAAT')
            if count_5A >= 1:
                if (count_4A - count_5A + count_GAAAT) >= threshold:
                    key_mask[a] = True
        kis_mask_seconds += time.time() - t0

        matched_payloads = payload[key_mask, :]  # shape: (M, 78)
        if len(matched_payloads) == 0:
            continue

        # Reference-pattern matching at this threshold.
        # MATLAB: for b = 73 : -1 : 58 (1-indexed in 78bp payload)
        # Python: b = 72 down to 57 (0-indexed)
        # Boundary checks removed — b range guarantees indices stay valid.
        t1 = time.time()
        Ks_list = []
        for row in matched_payloads:
            temp = ''.join(row)
            matched = False
            for b in range(72, 56, -1):
                keyi = len(temp[b + 1:])  # length of temp(b+1:end)
                # poly-A check: temp[b-4:b+1] == 'AAAAA'
                if temp[b - 4:b + 1] == 'A' * 5:
                    # reference pattern match
                    if temp[b - keyi - 10:b - keyi] == mDNA[b - keyi - 10:b - keyi]:
                        matched = True
                        # MATLAB overwrites KeyD{a}; no break — semantic identical
            if matched:
                Ks_list.append(temp)
        ref_match_seconds += time.time() - t1

        if len(Ks_list) > 0:
            used_threshold = threshold
            print(f"  KIS adaptive: threshold={threshold}, KIS-mask "
                  f"{int(key_mask.sum())}, reference-match {len(Ks_list)}")
            print(f"  KIS sub-timing: KIS-mask scan {kis_mask_seconds:.3f}s, "
                  f"reference-match {ref_match_seconds:.3f}s")
            break

    # publish sub-timings for the caller via module-level dict.
    KEYFIND_LAST_TIMING.clear()
    KEYFIND_LAST_TIMING.update({
        'kis_mask_s': round(kis_mask_seconds, 3),
        'ref_match_s': round(ref_match_seconds, 3),
        'used_threshold': used_threshold,
        'n_reference_match': len(Ks_list),
    })

    if len(Ks_list) == 0:
        raise ValueError(
            "No valid Key sequence found (tried adaptive thresholds 3, 2, 1)."
        )

    # Consensus (mode per column)
    ks_matrix = np.array([list(s) for s in Ks_list])
    KmatDNA_arr = _char_mode_per_column(ks_matrix)
    KmatDNA = ''.join(KmatDNA_arr)

    # find the 'AAAA' delimiter and decode the prefix as Key.
    # MATLAB v240108: if KeyMaptemp contains 'N', skip this AAAA position and try
    # the next b (do NOT replace N with A and decode).
    Key = None
    for b in range(len(KmatDNA) - 3):
        if KmatDNA[b:b + 4] == 'AAAA':
            key_map_str = KmatDNA[0:b]
            if 'N' in key_map_str:
                # skip this AAAA position — try next b (matches MATLAB `continue`)
                continue
            # DNA -> binary -> decimal (bi2de, LSB-first)
            key_bits = _dna_to_binary(key_map_str)
            Key = _bi2de(key_bits)
            break

    if Key is None:
        raise ValueError("Could not extract Key value.")

    return Key, KmatDNA


def Keysort(keymap, KD, k, filename):
    """
    Extract mesh parameters (ncp, nbp, ncc, nbc) and oligolen, idxlen from the key.
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
        raise ValueError(f"Keymap digit count error: {strkey} ({len(strkey)} digits)")

    keysz = ncolposz + nbitposz + ncolcosz + nbitcosz

    # KD (Key DNA consensus) -> binary.
    # MATLAB v240108: ismember(KD,'ATGC') would error if N present (returns 0).
    # We raise explicitly to match MATLAB's hard-fail behavior on N.
    if 'N' in KD:
        raise ValueError(
            f"KeymatDNA contains 'N' — MATLAB ismember would fail here "
            f"(KD={KD!r}). Upstream KeyFind likely produced an unreliable consensus."
        )
    Keybin = _dna_to_binary(KD)

    # offset based on keymap bit length
    keymap_bits = _de2bi(keymap)
    keymap_bitlen = len(keymap_bits)

    if keymap > 120000:
        start = keymap_bitlen + 8  # MATLAB 1-indexed +9 -> Python 0-indexed +8
    else:
        start = keymap_bitlen + 9  # MATLAB 1-indexed +10 -> Python 0-indexed +9

    # extract binkloc (keysz+1 bits)
    binkloc = Keybin[start:start + keysz + 1]

    # extract each parameter (bi2de, LSB-first)
    offset = 0
    ncp = _bi2de(binkloc[offset:offset + ncolposz])
    offset += ncolposz
    nbp = _bi2de(binkloc[offset:offset + nbitposz])
    offset += nbitposz
    ncc = _bi2de(binkloc[offset:offset + ncolcosz])
    offset += ncolcosz
    nbc = _bi2de(binkloc[offset:])

    # compute nfile, idxlen
    nfile = 3 * ncp * nbp + 3 * ncc * nbc
    tidx = 139

    while tidx + len(_de2bi(math.ceil(nfile / tidx * 1.5))) < k:
        tidx += 1

    idxlen = k - tidx

    # oligolen: last idxlen bits of the Key binary
    keyidx = _bi2de(Keybin[-idxlen:])
    oligolen = int(keyidx)

    # Final: removed `oligolen += 1` for pisa/seokga/liberty.
    # MATLAB Single_Decoding_Time.m v240108 line 211-212 sets
    #   `OligoLen = keyidx;` with no +1 adjustment.
    # The +1 in earlier ports caused tridecode reshape mismatches by adding
    # exactly `s_val` (= k - idxlen) extra elements to ECC_payload after the
    # tri_part slice + Null_count subtraction.

    return oligolen, idxlen, ncp, nbp, ncc, nbc


def Consensus(dna_matrix_list, idxlen, oligo_len):
    """
    Group reads by index.
    MATLAB: Consensus.m

    Args:
        dna_matrix_list: list of str (150 bp sequences)
        idxlen: int (index bit length)
        oligo_len: int (total number of oligos)

    Returns:
        cell_seq: list of np.array (per-index read group as a character matrix)
    """
    # build character matrix
    dna_matrix = np.array([list(s) for s in dna_matrix_list])

    # extract the index region: MATLAB DNAmatrix(:, 98-ceil(idxlen/2)+1 : 98)
    idx_end = 98  # MATLAB 1-indexed column 98 = Python 0-indexed column 97, but slice goes up to [:98]
    idx_start = 98 - math.ceil(idxlen / 2)  # MATLAB: 98-ceil(idxlen/2)+1 (1-indexed) = Python: 98-ceil(idxlen/2) (0-indexed)
    IndexDNA = dna_matrix[:, idx_start:idx_end]

    # DNA -> binary
    nrows, ncols = IndexDNA.shape
    Indexbin = np.zeros((nrows, ncols * 2), dtype=np.int8)

    for r in range(nrows):
        row_str = ''.join(IndexDNA[r, :])
        if 'N' in row_str:
            row_str = 'A' * ncols
        Indexbin[r, :] = _dna_to_binary(row_str)

    # parity alignment: MATLAB rem(length(de2bi(OligoLen)),2)
    oligo_bits = _de2bi(oligo_len)
    if len(oligo_bits) % 2 != 0:
        if Indexbin.shape[1] % 2 == 0:
            Indexbin = Indexbin[:, 1:]  # drop the first column

    # bi2de (per row)
    Indexdec = _bi2de(Indexbin)

    # group by index (MATLAB CellSeq{a} = DNAmatrix(Indexdec==a, :), a=1:OligoLen)
    cell_seq = [None] * oligo_len
    for a in range(1, oligo_len + 1):
        mask = (Indexdec == a)
        if np.any(mask):
            cell_seq[a - 1] = dna_matrix[mask, :]
        else:
            cell_seq[a - 1] = np.empty((0, dna_matrix.shape[1]), dtype='U1')

    return cell_seq


# ==============================================================================
# Phase 3: RS decoding
# ==============================================================================

def Binarization(cell_seq):
    """
    Convert per-index group consensus to a binary matrix.
    MATLAB: Binarization.m

    Returns:
        binSeq: np.array shape (OligoLen, 220) of int8
    """
    n_oligo = len(cell_seq)

    # compute mode (most frequent base) per group
    data_mode = []
    for a in range(n_oligo):
        group = cell_seq[a]
        if group is not None and len(group) > 0:
            consensus = _char_mode_per_column(group)
            data_mode.append(consensus)
        else:
            # if empty, fill with all-A (150 bp)
            data_mode.append(np.array(list('A' * 150)))

    dna_data = np.array(data_mode)  # shape: (OligoLen, 150)

    # extract data region: positions 21-130 (MATLAB) = 20:130 (Python), 110 columns
    mat_data = dna_data[:, 20:130]

    # replace N -> A then convert to binary
    binSeq = np.zeros((mat_data.shape[0], mat_data.shape[1] * 2), dtype=np.int8)

    for r in range(mat_data.shape[0]):
        row_str = ''.join(mat_data[r, :])
        row_str = row_str.replace('N', 'A')
        binSeq[r, :] = _dna_to_binary(row_str)

    return binSeq


def RSdecoding(bS, k, m, n):
    """
    Reed-Solomon error-correction decoding.
    MATLAB: RSdecoding.m (gf + rsdec)

    Primary : galois library
    Fallback: MATLAB Engine

    Args:
        bS: np.array shape (N, 220) - binary sequence matrix
        k: int (156) - number of info symbols
        m: int (8) - bits per symbol
        n: int (164) - codeword length

    Returns:
        rx: np.array shape (N, k) - decoded data
        err: np.array shape (N,) - error count (-1 = uncorrectable)
    """
    N = bS.shape[0]

    # split into payload (first k columns) and redundancy (remaining columns)
    bP = bS[:, :k]  # shape: (N, 156)
    bR = bS[:, k:]  # shape: (N, 64) = 8 symbols * 8 bits

    # convert redundancy bits into decimal symbols in 8-bit chunks (bi2de, LSB-first)
    n_parity = (bS.shape[1] - k) // m  # = 8
    dR = np.zeros((N, n_parity), dtype=np.int64)
    for j in range(n_parity):
        dR[:, j] = _bi2de(bR[:, j * m:(j + 1) * m])

    if HAS_GALOIS:
        return _rs_decode_galois(bP, dR, k, m, n, N)
    elif HAS_MATLAB:
        return _rs_decode_matlab(bP, dR, k, m, n, N)
    else:
        raise RuntimeError("RS decoding unavailable: requires galois or MATLAB Engine.")


def _rs_decode_galois(bP, dR, k, m, n, N):
    """
    Shortened RS decoding via the galois library.
    RS(164,156) over GF(2^8) = shortened RS(255,247).
    91-symbol zero padding -> RS(255,247) decode -> remove padding
    """
    GF = galois.GF(2 ** m)

    # full RS codec: RS(255, 247) over GF(2^8)
    n_full = 2 ** m - 1  # 255
    k_full = n_full - (n - k)  # 255 - 8 = 247
    pad_len = n_full - n  # 255 - 164 = 91

    rs_full = galois.ReedSolomon(n_full, k_full, field=GF)

    rx = np.zeros((N, k), dtype=np.float64)
    err = np.zeros(N, dtype=np.int64)

    # chunked processing (1000 rows at a time, MATLAB-compatible)
    chunk_size = 1000

    for start in range(0, N, chunk_size):
        end = min(start + chunk_size, N)
        chunk_bP = bP[start:end]
        chunk_dR = dR[start:end]
        chunk_n = end - start

        # build codeword: payload(0/1) + parity symbols(0-255)
        codeword_short = np.hstack([
            chunk_bP.astype(np.int64),
            chunk_dR.astype(np.int64)
        ])  # shape: (chunk_n, 164)

        # 91-symbol zero padding -> 255 symbols
        padding = np.zeros((chunk_n, pad_len), dtype=np.int64)
        codeword_full = np.hstack([padding, codeword_short])  # (chunk_n, 255)

        # decode per row
        for i in range(chunk_n):
            try:
                gf_cw = GF(codeword_full[i:i + 1, :])
                decoded_full, n_errors = rs_full.decode(gf_cw, errors=True)
                # remove padding: skip first pad_len(91) symbols
                decoded_short = decoded_full[0, pad_len:].view(np.ndarray).astype(np.float64)
                rx[start + i, :] = decoded_short[:k]
                err[start + i] = int(n_errors[0])
            except Exception:
                # uncorrectable error
                rx[start + i, :] = bP[start + i].astype(np.float64)
                err[start + i] = -1

    return rx, err


def _rs_decode_matlab(bP, dR, k, m, n, N):
    """RS decoding via MATLAB Engine (fallback)."""
    eng = matlab.engine.start_matlab()
    try:
        rx = np.zeros((N, k), dtype=np.float64)
        err = np.zeros(N, dtype=np.int64)

        for start in range(0, N, 1000):
            end_idx = min(start + 1000, N)
            bP_chunk = matlab.double(bP[start:end_idx].tolist())
            dR_chunk = matlab.double(dR[start:end_idx].tolist())

            # set variables in the MATLAB workspace
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
# Phase 4: Null + XOR recovery
# ==============================================================================

def NullFind(seq, Nc, Nxc, s, Fp):
    """
    Find Null-pattern (all-A) sequences and return data/XOR null sequences.
    MATLAB: NullFind.m

    Returns:
        Ndata: str (data null consensus)
        Xdata: str (XOR null consensus)
        dNc: list (data null prefix list)
        xNc: list (XOR null prefix list)
    """
    # Compute null range (parity-based branching from MATLAB v240108).
    # MATLAB note: when s is odd, MATLAB's `s/2` is fractional. For odd s + odd Nc
    # the fractional halves cancel, leaving an integer. Python's s_half = s // 2
    # plus Nc_half = Nc // 2 gives the same integer result because:
    #   (floor(s/2)+0.5) - (floor(Nc/2)+0.5) = floor(s/2) - floor(Nc/2)
    s_half = s // 2  # floor(s/2)
    Nc_half = Nc // 2
    Nxc_half = Nxc // 2

    if s % 2 != 0 and Nc % 2 != 0 and Nxc % 2 != 0:
        # branch 1: all-odd. MATLAB uses fractional s/2 - Nc/2 which simplifies
        # to floor(s/2) - floor(Nc/2). Python equivalent below.
        dNr_start = s_half - Nc_half + 1 + 20  # MATLAB 1-indexed + 20
        dNr_end = s_half + 20
        xNr_start = s_half - Nxc_half + 1 + 20
        xNr_end = s_half + 20
    elif (s % 2 != 0 and Nc % 2 == 0 and Nxc % 2 == 0) or \
         (s % 2 == 0 and Nc % 2 == 0 and Nxc % 2 == 0):
        # branch 2: Nc even AND Nxc even. floor(s/2) - Nc/2 with both halves integer.
        dNr_start = s_half - Nc_half + 1 + 20
        dNr_end = s_half + 20
        xNr_start = s_half - Nxc_half + 1 + 20
        xNr_end = s_half + 20
    else:
        # branch 3: at least one of (Nc, Nxc) is odd.
        # MATLAB v240108 uses Null_count/2 directly which yields a fractional index
        # when Nc is odd — producing a MATLAB error. Python applies ceil(Nc/2) to
        # mirror MATLAB's xor-side handling and stay defined.
        dNr_start = s_half - math.ceil(Nc / 2) + 1 + 20
        dNr_end = s_half + 20
        xNr_start = s_half - math.ceil(Nxc / 2) + 1 + 20
        xNr_end = s_half + 20

    # convert to Python 0-indexed
    dNr = list(range(dNr_start - 1, dNr_end))  # MATLAB 1-indexed -> Python 0-indexed
    xNr = list(range(xNr_start - 1, xNr_end))

    dpA = 'A' * len(dNr)
    xpA = 'A' * len(xNr)

    # detect null patterns in each read
    dNo = []
    xNo = []
    for a in range(len(seq)):
        if len(seq[a]) >= s / 2 + 20:  # float compare to match MATLAB behavior (handles odd s)
            # check data null
            d_region = ''.join(seq[a][i] for i in dNr if i < len(seq[a]))
            if d_region == dpA:
                # check primer match (>=15/20)
                fwd_match = sum(1 for i in range(min(20, len(Fp), len(seq[a])))
                                if seq[a][i] == Fp[i])
                if fwd_match >= 15:
                    dNo.append(a)

            # check XOR null
            x_region = ''.join(seq[a][i] for i in xNr if i < len(seq[a]))
            if x_region == xpA:
                fwd_match = sum(1 for i in range(min(20, len(Fp), len(seq[a])))
                                if seq[a][i] == Fp[i])
                if fwd_match >= 15:
                    xNo.append(a)

    # deduplicate: if xNr is longer than dNr, drop data-null entries that also match xor-null
    dNs = [seq[i] for i in dNo]
    xNs = [seq[i] for i in xNo]

    if len(xNr) > len(dNr):
        dNs_clean = []
        for s_item in dNs:
            x_region = ''.join(s_item[i] for i in xNr if i < len(s_item))
            if x_region != xpA:
                dNs_clean.append(s_item)
        dNs = dNs_clean
    else:
        xNs_clean = []
        for s_item in xNs:
            d_region = ''.join(s_item[i] for i in dNr if i < len(s_item))
            if d_region != dpA:
                xNs_clean.append(s_item)
        xNs = xNs_clean

    # extract prefixes
    min_dNr = min(dNr) if dNr else 0
    min_xNr = min(xNr) if xNr else 0

    dNc = [s_item[:min_dNr] for s_item in dNs if len(s_item) > min_dNr]
    xNc = [s_item[:min_xNr] for s_item in xNs if len(s_item) > min_xNr]

    # Consensus (mode)
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
    XOR-based error recovery.
    MATLAB: XORrun.m

    Returns:
        xtp, md, tp, xp, XA, XB, XC
    """
    total = matrix.size
    rows = matrix.shape[0]

    # split 2/3 : 1/3
    tri_rows = total * 2 // (3 * s)
    tp = matrix[:tri_rows, :].copy()
    xp = matrix[tri_rows:, :].copy()

    # flatten (row-major = MATLAB reshape(tp', [numel(tp) 1])')
    row_part = tp.flatten(order='C')
    if Nc > 0:
        row_part = row_part[:-Nc]

    # handle odd length
    if (ncp * 3 * nbp + ncc * 3 * nbc) % 2 != 0:
        row_part = row_part[:-1]

    XC = xp.flatten(order='C')
    if Nxc > 0:
        XC = XC[:-Nxc]

    # split A/B
    half = math.ceil(len(row_part) / 2)
    XA = row_part[:half]
    XB = row_part[half:]

    if len(row_part) % 2 != 0:
        XB = np.append(XB, 0)

    # XOR recovery
    XOR_matrix = np.concatenate([
        np.bitwise_xor(XB.astype(np.int64), XC[:len(XB)].astype(np.int64)),
        np.bitwise_xor(XA.astype(np.int64), XC[:len(XA)].astype(np.int64))
    ])

    # pad then reshape
    ECC_example = np.concatenate([XOR_matrix, np.zeros(Nc, dtype=np.int64)])
    ECC_ex_matrix = ECC_example.reshape((-1, s))  # -1 = numel/s

    # find missing rows (rows whose sum is 0)
    md = np.where(tp.sum(axis=1) == 0)[0]

    # replace missing rows with XOR-recovered data
    for idx in md:
        if idx < len(ECC_ex_matrix):
            tp[idx, :] = ECC_ex_matrix[idx, :]

    xtp = tp.copy()
    return xtp, md, tp, xp, XA, XB, XC


def nestedXOR(CellSeq, xtp, ErrD, RSerrinf, XA, XB, XC, payload, tp, xp, matsz, Null_count):
    """
    Nested XOR recovery: only recover when the corresponding XOR partner is error-free.
    MATLAB: nestedXOR.m

    Returns:
        XORM: list of np.array (read groups targeted for recovery)
        Dx: np.array (indices targeted for recovery)
    """
    # Depth_cut error indices (only those that still have RS errors)
    err_mask = RSerrinf[ErrD] != 0
    Depth_cut = ErrD[err_mask]

    if len(Depth_cut) == 0:
        return [], np.array([], dtype=np.int64)

    xp_rows = xp.shape[0]
    tp_rows = tp.shape[0]

    # A/B/C partitions (kept against MATLAB 1-based indexing convention)
    A_part = Depth_cut[Depth_cut <= xp_rows]
    C_part = Depth_cut[Depth_cut > tp_rows]
    bc_mask = (Depth_cut > len(A_part)) if len(A_part) > 0 else np.ones(len(Depth_cut), dtype=bool)
    remaining = Depth_cut[len(A_part):]
    B_part = remaining[:len(remaining) - len(C_part)] if len(C_part) > 0 else remaining

    # adjust trailing indices for B_part, C_part
    if len(B_part) > 0 and B_part[-1] == tp_rows:
        B_part = B_part[:-1]
    if len(C_part) > 0 and C_part[-1] == tp_rows + xp_rows - 1:
        C_part = C_part[:-1]

    Xi = np.zeros(len(Depth_cut), dtype=np.int64)

    # A part check
    for a_idx in range(len(A_part)):
        a = A_part[a_idx]  # 0-based index
        AtoB_start = a * matsz - len(XA) % matsz
        AtoB_end = AtoB_start + matsz
        AtoC_start = a * matsz
        AtoC_end = AtoC_start + matsz

        if AtoB_start >= 0 and AtoB_end <= len(XB) and AtoC_start >= 0 and AtoC_end <= len(XC):
            if np.sum(XB[AtoB_start:AtoB_end]) == 0 or np.sum(XC[AtoC_start:AtoC_end]) == 0:
                Xi[a_idx] = a

    # B part check (MATLAB B_part(b) is 1-based, Python is 0-based -> +1 adjustment)
    for b_idx in range(len(B_part)):
        b = B_part[b_idx]
        BtoA_start = (b + 1) * matsz - len(XA) + (matsz - len(XA) % matsz)
        BtoA_end = BtoA_start + matsz
        BtoC_start = BtoA_start
        BtoC_end = BtoA_end

        if BtoA_start >= 0 and BtoA_end <= len(XA) and BtoC_start >= 0 and BtoC_end <= len(XC):
            if np.sum(XA[BtoA_start:BtoA_end]) == 0 or np.sum(XC[BtoC_start:BtoC_end]) == 0:
                Xi[len(A_part) + b_idx] = b

    # C part check
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

    # if RS error == 0, restore from original payload; if -1, needs re-decoding
    for idx in Non_XOR_oligo:
        if idx < len(RSerrinf) and idx < xtp.shape[0]:
            if RSerrinf[idx] == 0 and idx < xtp.shape[0]:
                xtp[idx, :] = payload[idx, :matsz] if idx < payload.shape[0] else xtp[idx, :]

    # indices that still need re-decoding
    rs_err_mask = np.array([RSerrinf[int(i)] != 0 for i in Non_XOR_oligo
                            if int(i) < len(RSerrinf)])
    Dx = Non_XOR_oligo[rs_err_mask] if len(rs_err_mask) > 0 else np.array([], dtype=np.int64)

    # extract original read groups at those indices
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
    Post-XOR re-RS decoding: re-binarize original reads and run RS decoding.
    MATLAB: reRSdec.m
    """
    for a in range(len(XORM)):
        group = XORM[a]
        if group is None or len(group) == 0:
            continue

        # extract data region (positions 21-130 MATLAB) = (20:130 Python)
        reRSdata = group[:, 20:130] if group.shape[1] >= 130 else group

        if reRSdata.size == 0:
            continue

        # binarize
        nrows_g, ncols_g = reRSdata.shape
        reSeqtoBin = np.zeros((nrows_g, ncols_g * 2), dtype=np.int8)

        for r in range(nrows_g):
            row_str = ''.join(reRSdata[r, :])
            row_str = row_str.replace('N', 'A')
            reSeqtoBin[r, :] = _dna_to_binary(row_str)

        # RS decode
        try:
            re_rx, re_err = RSdecoding(reSeqtoBin, k, m, n)

            # use only successful decodings (error count == 0)
            success_mask = (re_err == 0)
            if np.any(success_mask):
                successful = re_rx[success_mask, :matsz]
                # mode (most frequent value)
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
# Phase 5: 3D reconstruction
# ==============================================================================

def tridecode(ECCp, ncp, nbp, ncc, nbc, filename):
    """
    Decode 3D coordinates (tp) and face connectivity (tc) from the binary payload.
    MATLAB: tridecode.m

    Returns:
        tc: np.array shape (ncc, 3) - face indices (1-based)
        tp: np.array shape (ncp, 3) - vertex coordinates
    """
    ECCp = np.asarray(ECCp, dtype=np.int64)

    # handle odd length
    if (ncp * 3 * nbp + ncc * 3 * nbc) % 2 != 0:
        ECCp = ECCp[:-1]

    # split into coordinate region (pr) and connectivity region (cr)
    pr = ECCp[:ncp * nbp * 3]
    cr = ECCp[ncp * nbp * 3:]

    # coordinates: reshape -> bi2de -> /10^16
    # MATLAB: pb = reshape(pr', [nbp ncp*3])'; -> (ncp*3, nbp) matrix
    pb = np.reshape(pr, (ncp * 3, nbp), order='C')  # row-major (pr is already flat row-major)
    pd = _bi2de(pb) / (10 ** 16)

    # reshape into 3D coordinates
    # special datasets: reshape(pd, [ncp 3]) in MATLAB column-major
    if filename in ('seokga', 'liberty', 'pisa', 'venus', 'goldengate', 'tripitaka'):
        tp = np.reshape(pd, (ncp, 3), order='F')
    else:
        # default: reshape(pd, [3 ncp])' -> (ncp, 3)
        tp = np.reshape(pd, (3, ncp), order='F').T

    # connectivity: reshape -> bi2de
    cb = np.reshape(cr, (ncc * 3, nbc), order='C')
    cd = _bi2de(cb)
    tc = np.reshape(cd, (3, ncc), order='F').T

    # safe defaults: 0 -> 1, out-of-range -> 1
    tc[tc == 0] = 1
    tc[tc > len(tp)] = 1

    return tc, tp


def trierrXOR(mat, tc, tp, Ncount, Nxcount, ncp, nbp, ncc, nbc):
    """
    Recover zero-valued vertices/faces using XOR data.
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

    # recover zero-valued vertices
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

    # recover zero-valued / out-of-range faces
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

        # drop faces still out of range
        bad_rows = np.where(
            (tc > ncp).any(axis=1) | (tc == 0).any(axis=1)
        )[0]
        tc = np.delete(tc, bad_rows, axis=0)

    return _make_triangulation(tc, tp)


def _make_triangulation(tc, tp):
    """
    Build a 3D mesh object from vertices and faces.
    MATLAB: triangulation(tc, tp)
    """
    # MATLAB face indices are 1-based -> Python 0-based
    tc_0based = tc.astype(int) - 1

    # validity check
    max_idx = tp.shape[0] - 1
    valid_mask = np.all((tc_0based >= 0) & (tc_0based <= max_idx), axis=1)
    tc_0based = tc_0based[valid_mask]

    if HAS_TRIMESH:
        mesh = trimesh.Trimesh(vertices=tp, faces=tc_0based)
        return mesh
    else:
        return {'vertices': tp, 'faces': tc_0based, 'faces_1based': tc}


def visualize_mesh(trid, title_str='Decoded 3D Mesh', save_png_path=None):
    """3D mesh visualization (matplotlib). v8_10: PNG save only — no plt.show()
    (Agg backend; plt.show() blocks for minutes in headless runs)."""
    if HAS_TRIMESH and isinstance(trid, trimesh.Trimesh):
        vertices = trid.vertices
        faces = trid.faces
    elif isinstance(trid, dict):
        vertices = trid['vertices']
        faces = trid['faces']
    else:
        print("[warn] no mesh data to visualize.")
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
    if save_png_path is not None:
        try:
            plt.savefig(save_png_path, dpi=120, bbox_inches='tight')
            print(f"  PNG saved: {save_png_path}")
        except Exception as e:
            print(f"  [warn] PNG save failed: {e}")
    plt.close(fig)


# ==============================================================================
# PHASE 6: MAIN PIPELINE
# ==============================================================================

def read_fastq_file(filepath):
    """Read a FASTQ file (Biopython preferred, plain-text fallback). Supports .gz."""
    import gzip
    sequences = []
    if not os.path.exists(filepath):
        print(f"[warn] file not found: {filepath}")
        return sequences

    is_gz = filepath.endswith('.gz')

    if HAS_BIOPYTHON:
        try:
            handle = gzip.open(filepath, 'rt') if is_gz else open(filepath, 'r')
            with handle as fh:
                for record in SeqIO.parse(fh, "fastq"):
                    sequences.append(str(record.seq))
            return sequences
        except Exception as e:
            print(f"[warn] Biopython parse failed, falling back to text mode: {e}")

    # Plain-text fallback (MATLAB textscan style: line 2,6,10,... = sequence)
    opener = gzip.open if is_gz else open
    with opener(filepath, 'rt') as f:
        lines = f.readlines()
    sequences = [lines[i].strip() for i in range(1, len(lines), 4)]
    return sequences


def main():
    """
    DNA decoding main pipeline.
    MATLAB equivalent: Textscan_230527_downsampling_additional.m
    Includes per-step timing.
    """
    # =========================================================================
    # Run config — filename can be overridden via CLI arg.
    # Usage: python3 dna_decode_baseline_260506_1020.py <dataset>
    # =========================================================================
    filename = sys.argv[1] if len(sys.argv) > 1 else 'Jikji6'
    Timecell = {}

    print(f"=== DNA Decoding Pipeline: {filename} ===\n")
    print("[Step 0] Loading FASTQ files...")

    # =========================================================================
    # Dataset -> paired FASTQ resolution (if/elif/else style).
    # Paths are resolved against '3D culture heritage seq data/' subfolders.
    # J1..J7 are short aliases for Jikji1..Jikji7.
    # read_fastq_file() also accepts .gz inputs transparently.
    # =========================================================================
    DATA_ROOT = '/Users/young/Desktop/DNA ink3/3D culture heritage seq data'

    if filename in ('Jikji1', 'Jikji2', 'Jikji3', 'Jikji4', 'Jikji5', 'Jikji6', 'Jikji7'):
        folder = 'Jikji'
        r1_path = os.path.join(DATA_ROOT, folder, f'{filename}_1.fastq')
        r2_path = os.path.join(DATA_ROOT, folder, f'{filename}_2.fastq')
    elif filename in ('J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7'):
        # short alias: J{N} -> Jikji{N}
        folder = 'Jikji'
        n = filename[1:]
        r1_path = os.path.join(DATA_ROOT, folder, f'Jikji{n}_1.fastq')
        r2_path = os.path.join(DATA_ROOT, folder, f'Jikji{n}_2.fastq')
    elif filename == 'tripitaka':
        # The original tripitaka_{1,2}.fastq{.gz} files were tar-wrapped
        # (POSIX 'PaxHeader/...' prefix); they were extracted to a clean
        # 'Tripitaka_extracted/' subfolder once at setup time.
        folder = 'Tripitaka_extracted'
        r1_path = os.path.join(DATA_ROOT, folder, 'tripitaka_1.fastq')
        r2_path = os.path.join(DATA_ROOT, folder, 'tripitaka_2.fastq')
    elif filename == 'liberty':
        folder = 'Statue of Liberty'
        r1_path = os.path.join(DATA_ROOT, folder, 'liberty_1.fastq')
        r2_path = os.path.join(DATA_ROOT, folder, 'liberty_2.fastq')
    elif filename == 'pisa':
        folder = 'Tower of Pisa'
        r1_path = os.path.join(DATA_ROOT, folder, 'pisa_1.fastq')
        r2_path = os.path.join(DATA_ROOT, folder, 'pisa_2.fastq')
    elif filename == 'seokga':
        folder = 'seokgatap'
        r1_path = os.path.join(DATA_ROOT, folder, 'seokga_1.fastq')
        r2_path = os.path.join(DATA_ROOT, folder, 'seokga_2.fastq')
    elif filename == 'venus':
        folder = 'Venus Statue'
        r1_path = os.path.join(DATA_ROOT, folder, 'V_1.fastq')
        r2_path = os.path.join(DATA_ROOT, folder, 'V_2.fastq')
    elif filename == 'goldengate':
        folder = 'Golden gate'
        r1_path = os.path.join(DATA_ROOT, folder, 'goldengate_1.fastq')
        r2_path = os.path.join(DATA_ROOT, folder, 'goldengate_2.fastq')
    else:
        # generic fallback: assume the file lives directly under DATA_ROOT
        # following the {name}_1.fastq / {name}_2.fastq convention.
        r1_path = os.path.join(DATA_ROOT, f'{filename}_1.fastq')
        r2_path = os.path.join(DATA_ROOT, f'{filename}_2.fastq')

    seq = read_fastq_file(r1_path) + read_fastq_file(r2_path)

    print(f"  loaded reads: {len(seq):,}")

    # =========================================================================
    # Step 1: Evert Sequence (orient correction)
    # =========================================================================
    t = time.time()
    print("\n[Step 1] Evert Sequence...")

    forward_section = [s[:21] if len(s) >= 21 else s for s in seq]
    F_primers, R_primers = PrimerSelect(filename)
    tevert_seq, reversed_seq_log = ReverseEvert(seq, forward_section, R_primers, F_primers)

    # replace reverse-orientation reads
    rev_indices = np.where(reversed_seq_log)[0]
    for i, idx in enumerate(rev_indices):
        if tevert_seq[i] is not None:
            seq[idx] = tevert_seq[i]

    Timecell['Evert'] = round(time.time() - t, 3)
    print(f"  done: {Timecell['Evert']}s (reverse reads: {len(rev_indices):,})")

    # =========================================================================
    # Step 2: Primer Sorting
    # =========================================================================
    t = time.time()
    print("\n[Step 2] Primer Sorting...")

    PrimerSortingSeq = PrimerSort(seq, F_primers, R_primers)

    Timecell['Primer Sorting'] = round(time.time() - t, 3)
    print(f"  done: {Timecell['Primer Sorting']}s (valid reads: {len(PrimerSortingSeq):,})")

    if len(PrimerSortingSeq) == 0:
        print("[error] no valid reads. Pipeline aborted.")
        return Timecell

    # =========================================================================
    # Step 3: Key Finding
    # =========================================================================
    t = time.time()
    print("\n[Step 3] Key Finding...")

    k = 156
    n = k + 8  # 164
    m = len(_de2bi(n))  # 8

    KeyMap, KeymatDNA = KeyFind(PrimerSortingSeq)
    OligoLen, idxlen, ncolpo, nbitpo, ncolco, nbitco = Keysort(KeyMap, KeymatDNA, k, filename)

    Timecell['Key Finding'] = round(time.time() - t, 3)
    print(f"  done: {Timecell['Key Finding']}s")
    print(f"  KeyMap={KeyMap}, OligoLen={OligoLen}, idxlen={idxlen}")
    print(f"  ncolpo={ncolpo}, nbitpo={nbitpo}, ncolco={ncolco}, nbitco={nbitco}")

    # =========================================================================
    # Step 4: Consensus Sequencing (Index Sorting)
    # =========================================================================
    t = time.time()
    print("\n[Step 4] Consensus Sequencing...")

    MIndex = Consensus(PrimerSortingSeq, idxlen, OligoLen)

    # compute Depth
    Depth = np.zeros(len(MIndex), dtype=np.int64)
    for i in range(len(MIndex)):
        Depth[i] = MIndex[i].shape[0] if MIndex[i] is not None else 0

    Timecell['Index Sorting'] = round(time.time() - t, 3)
    print(f"  done: {Timecell['Index Sorting']}s")
    print(f"  mean Depth: {np.mean(Depth):.1f}x, min: {np.min(Depth)}, max: {np.max(Depth)}")

    # Mapping Location (only for pisa, seokga, liberty).
    # MATLAB Single_Decoding_Time.m v240108 lines 252-278.
    # For each read in MIndex[-1] (the last index group), scan payload positions
    # 21..98 (MATLAB) for three non-overlapping 'AAAA' patterns. Keep only reads
    # with exactly 3 AAAA hits as the true KeyIndex; drop the rest.
    if filename in ('pisa', 'seokga', 'liberty'):
        if MIndex and MIndex[-1] is not None and len(MIndex[-1]) > 0:
            key_map = MIndex[-1]
            n_reads = key_map.shape[0]
            KeyIndex_count = np.zeros(n_reads, dtype=np.int64)

            for a in range(n_reads):
                # 1st AAAA: b in MATLAB 21..98 -> Python 0-idx 20..97
                b_found = -1
                for b in range(20, 98):
                    if (key_map[a, b] == 'A' and key_map[a, b + 1] == 'A'
                            and key_map[a, b + 2] == 'A' and key_map[a, b + 3] == 'A'):
                        KeyIndex_count[a] = 1
                        b_found = b
                        break
                if b_found < 0:
                    continue

                # 2nd AAAA: c starts at b+3 (MATLAB) = Python b_found+3
                c_found = -1
                for c in range(b_found + 3, 98):
                    if (key_map[a, c] == 'A' and key_map[a, c + 1] == 'A'
                            and key_map[a, c + 2] == 'A' and key_map[a, c + 3] == 'A'):
                        KeyIndex_count[a] += 1
                        c_found = c
                        break
                if c_found < 0:
                    continue

                # 3rd AAAA: d starts at c+3
                for d in range(c_found + 3, 98):
                    if (key_map[a, d] == 'A' and key_map[a, d + 1] == 'A'
                            and key_map[a, d + 2] == 'A' and key_map[a, d + 3] == 'A'):
                        KeyIndex_count[a] += 1
                        break

            # keep only reads with exactly 3 AAAA hits
            keep_mask = (KeyIndex_count == 3)
            MIndex[-1] = key_map[keep_mask, :]
            print(f"  MapLocation: kept {int(keep_mask.sum())}/{n_reads} key-index reads")

    # =========================================================================
    # Step 5: RS Decoding
    # =========================================================================
    t = time.time()
    print("\n[Step 5] RS Decoding...")

    SeqtoBin = Binarization(MIndex)
    rsECC, RSerrinf = RSdecoding(SeqtoBin, k, m, n)

    # error filtering
    rsECC[RSerrinf < 0, :] = 0
    rsECC[RSerrinf > 4, :] = 0
    for a in range(min(OligoLen, len(rsECC))):
        if np.any(rsECC[a, :] > 1):
            rsECC[a, :] = 0

    # extract payload
    payload = rsECC[:, :k - idxlen]
    Temp_mat = payload[:-1, :].copy()

    Timecell['RS Decoding'] = round(time.time() - t, 3)
    print(f"  done: {Timecell['RS Decoding']}s")
    print(f"  RS uncorrectable: {np.sum(RSerrinf < 0)}, corrected: {np.sum(RSerrinf > 0)}")

    # =========================================================================
    # Step 6: Null Finding
    # =========================================================================
    t = time.time()
    print("\n[Step 6] Null Finding...")

    # compute Nullkey
    term = ncolpo * 3 * nbitpo + ncolco * 3 * nbitco
    Nullkey_1 = term + 1 if term % 2 == 1 else term

    s_val = Temp_mat.shape[1]
    Nullkey_2 = Nullkey_1 // 2  # integer division (avoid float propagation)

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

    data_Ndecode, xor_Ndecode, dNc, xNc = NullFind(seq, Null_count, Null_xor_count, s_val, F_primers)

    # convert Null sequences to binary
    bindataN = _dna_to_binary(data_Ndecode[20:])  # MATLAB: data_Ndecode(21:end)
    binxorN = _dna_to_binary(xor_Ndecode[20:])    # MATLAB: xor_Ndecode(21:end)

    Timecell['Null Finding'] = round(time.time() - t, 3)
    print(f"  done: {Timecell['Null Finding']}s")

    # update Depth
    tri_row_idx = Temp_mat.size * 2 // (3 * s_val) - 1  # MATLAB 1-indexed -> Python 0-indexed
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

    # insert Null data into Temp_mat
    Temp_mat = payload[:-1, :].copy()
    data_insert_row = Temp_mat.size * 2 // (3 * s_val) - 1  # MATLAB 1-indexed -> Python 0-indexed
    if data_insert_row < Temp_mat.shape[0]:
        insert_len = min(len(bindataN), Temp_mat.shape[1])
        Temp_mat[data_insert_row, :insert_len] = bindataN[:insert_len]
    if Temp_mat.shape[0] > 0:
        insert_len = min(len(binxorN), Temp_mat.shape[1])
        Temp_mat[-1, :insert_len] = binxorN[:insert_len]

    # align Depth size
    if len(Depth) > Temp_mat.shape[0]:
        Depth = Depth[:Temp_mat.shape[0]]

    # === Recovery stage 1: Global RS ===
    trid = None
    try:
        print("  [1/3] Trying Global RS decoding...")
        ECC_tri_part = Temp_mat[:Temp_mat.size * 2 // (3 * s_val), :]
        ECC_payload_flat = ECC_tri_part.flatten(order='C')
        if Null_count > 0:
            ECC_payload_flat = ECC_payload_flat[:-Null_count]
        tc, tp = tridecode(ECC_payload_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
        trid = _make_triangulation(tc, tp)
        print("  -> Global RS succeeded!")
    except Exception as e:
        print(f"  -> Global RS failed: {e}")

    # === Recovery stage 2: Depth cut ===
    Depthcut = 16
    ErrD = np.where(Depth <= Depthcut)[0]

    try:
        print("  [2/3] Trying Depth-cut-based XOR recovery...")
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
        print("  -> Depth cut succeeded!")
    except Exception as e:
        print(f"  -> Depth cut failed: {e}")

    # === Recovery stage 3: Local RS (nestedXOR + reRSdec) ===
    try:
        print("  [3/3] Trying Local RS recovery...")
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

        # final recovery
        xtp_final, _, _, _, _, _, _ = \
            XORrun(Temp_mat, s_val, Null_count, Null_xor_count, ncolpo, ncolco, nbitpo, nbitco)
        ECC_payload_flat = xtp_final.flatten(order='C')
        if Null_count > 0:
            ECC_payload_flat = ECC_payload_flat[:-Null_count]
        tc, tp = tridecode(ECC_payload_flat, ncolpo, nbitpo, ncolco, nbitco, filename)
        trid = trierrXOR(Temp_mat, tc, tp, Null_count, Null_xor_count, ncolpo, nbitpo, ncolco, nbitco)
        print("  -> Local RS succeeded!")
    except Exception as e:
        print(f"  -> Local RS failed: {e}")
        if trid is None:
            print("  [error] Decoding Failed")

    Timecell['XOR Recovery'] = round(time.time() - t, 3)
    print(f"  done: {Timecell['XOR Recovery']}s")

    # =========================================================================
    # Step 8: STL Decode & Visualization
    # =========================================================================
    t = time.time()
    print("\n[Step 8] STL Decode & Visualization...")

    if trid is not None:
        # save STL
        output_dir = os.path.join(".", "Downsampling_Results", filename)
        os.makedirs(output_dir, exist_ok=True)

        if HAS_TRIMESH and isinstance(trid, trimesh.Trimesh):
            stl_path = os.path.join(output_dir, f"{filename}_decoded.stl")
            trid.export(stl_path)
            print(f"  STL saved: {stl_path}")
            print(f"  vertices: {len(trid.vertices):,}, faces: {len(trid.faces):,}")

        # visualization (PNG only, no GUI)
        try:
            png_path = os.path.join(output_dir, f"{filename}_decoded.png")
            visualize_mesh(trid, title_str=f'{filename} - Decoded 3D Mesh',
                           save_png_path=png_path)
        except Exception as e:
            print(f"  [warn] visualization failed: {e}")

    Timecell['STL Decode'] = round(time.time() - t, 3)

    # =========================================================================
    # Result summary
    # =========================================================================
    print("\n" + "=" * 50)
    print("Per-step timing (Timecell)")
    print("=" * 50)
    total_time = 0
    for step_name, elapsed in Timecell.items():
        print(f"  {step_name}: {elapsed}s")
        total_time += elapsed
    print(f"  {'-' * 30}")
    print(f"  Total: {round(total_time, 3)}s")

    # save Depth / error info to Excel
    output_dir = os.path.join(".", "Downsampling_Results", filename)
    os.makedirs(output_dir, exist_ok=True)

    fcoverage = round(np.mean(Depth), 1)

    try:
        depth_file = os.path.join(output_dir, f"{filename}_Depth.xlsx")
        pd.DataFrame(Depth, columns=['Depth']).to_excel(depth_file, sheet_name=f"{fcoverage}x")
        print(f"\nDepth saved: {depth_file}")
    except Exception as e:
        print(f"[warn] Depth save failed: {e}")

    # RS error indices
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
            print(f"RS errors saved: {err_file}")
        except Exception as e:
            print(f"[warn] RS error save failed: {e}")

    print(f"\nMean Depth: {fcoverage}x")
    print(f"Total reads: {len(seq):,}, Primer valid: {len(PrimerSortingSeq):,}")

    # v8_10: dump per-dataset metrics JSON for the Figure 2 Excel aggregator.
    metrics = {
        'dataset': filename,
        'version': 'Final-260507',
        'fastq_size_bytes': (
            (os.path.getsize(r1_path) if os.path.exists(r1_path) else 0) +
            (os.path.getsize(r2_path) if os.path.exists(r2_path) else 0)
        ),
        'reads_loaded': len(seq),
        'reads_valid': len(PrimerSortingSeq),
        'mean_depth': float(fcoverage),
        'oligo_len': int(OligoLen),
        'idxlen': int(idxlen),
        'ncolpo': int(ncolpo),
        'nbitpo': int(nbitpo),
        'ncolco': int(ncolco),
        'nbitco': int(nbitco),
        'vertices': int(trid.vertices.shape[0]) if (trid is not None and HAS_TRIMESH and isinstance(trid, trimesh.Trimesh)) else None,
        'faces': int(trid.faces.shape[0]) if (trid is not None and HAS_TRIMESH and isinstance(trid, trimesh.Trimesh)) else None,
        'timecell': {k: float(v) if v is not None else None for k, v in Timecell.items()},
        'keyfind_subtiming': dict(KEYFIND_LAST_TIMING),
    }
    metrics_dir = os.path.join('.', 'Downsampling_Results', filename)
    os.makedirs(metrics_dir, exist_ok=True)
    metrics_file = os.path.join(metrics_dir, f'{filename}_metrics_Final.json')
    try:
        with open(metrics_file, 'w', encoding='utf-8') as f:
            json.dump(metrics, f, indent=2, ensure_ascii=False)
        print(f"Metrics JSON saved: {metrics_file}")
    except Exception as e:
        print(f"[warn] metrics JSON save failed: {e}")

    return Timecell


if __name__ == "__main__":
    main()
