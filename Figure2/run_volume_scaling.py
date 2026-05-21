# =============================================================================
# Figure 2 — Data Volume Scaling (용량 증가 실험)
# =============================================================================
# [STATUS] 원본 독립 실행 스크립트는 02_Time_Quantification_Analysis 에 존재하지 않음.
#          Figure1 의 Python baseline 을 루프 래퍼로 확장하여 작성함.
# [OBJECTIVE] 데이터 용량을 늘려가며 단계별 디코딩 시간(초)을 측정해 선형/서브선형 확장성 검증.
# [DATASETS] Jikji6, Jikji7, Tripitaka (Tripitaka 는 아직 실측 안 된 extrapolation 대상).
# [DATA LINKS]
#   - Jikji6: 4_S2_L001_R1_001.fastq + 4_S2_L001_R2_001.fastq  (로컬 없음 — NAS 참조)
#   - Jikji7: <경로 미확정 — NAS 참조 필요>
#   - Tripitaka: <원본 FASTQ 로컬 없음>
# [TARGET METRICS] per-stage time (s), total time (s), oligo count, reads count.
# [REF] CONTENT.md 명시 실측치 — J6 10.1s, J7 17.5s (벡터화), 3.1× 데이터 → 1.7× 시간.
# =============================================================================
#
# TODO: wire to Figure1 baseline
#   이 스크립트는 Figure1/python/dna_decode_baseline.py 의 main() 로직을 래핑함.
#   실제 실행 전에 아래 두 가지를 확인/수행할 것:
#     1) FASTQ 파일을 NAS 에서 로컬 DATA_DIR 로 복사 또는 경로 수정
#     2) Figure1/python/ 경로가 sys.path 에 추가되어 import 가능한지 확인
# =============================================================================

import os
import sys
import time
import csv
import pathlib

# ---------------------------------------------------------------------------
# 경로 설정 — Figure1 baseline 을 import 하기 위한 sys.path 추가
# ---------------------------------------------------------------------------
_BASE_DIR = pathlib.Path(__file__).resolve().parent.parent  # DNA ink3/
_FIGURE1_PY = _BASE_DIR / "Figure1" / "python"
if str(_FIGURE1_PY) not in sys.path:
    sys.path.insert(0, str(_FIGURE1_PY))

# TODO: wire to Figure1 baseline
# Figure1/python/dna_decode_baseline.py 에 run_pipeline(r1_path, r2_path) 함수가
# 없으면 아래 import 는 실패함. 그 경우 baseline 을 함수로 리팩터링 필요.
try:
    from dna_decode_baseline import (
        read_fastq_file,
        reverse_evert_sequences,
        primer_sort,
        key_find,
        consensus_call,
        rs_decode,
        null_find,
        xor_recovery,
    )
    _BASELINE_AVAILABLE = True
except ImportError:
    _BASELINE_AVAILABLE = False
    print("[경고] Figure1 baseline import 실패 — stub 모드로 실행됩니다.")
    print("       Figure1/python/dna_decode_baseline.py 를 확인하세요.")


# =============================================================================
# 데이터셋 정의
# =============================================================================
# NAS 경로: C:\Users\BEBO-KJW\OneDrive - inu.ac.kr\NAS\KJW\Downsampling\
# 로컬로 복사한 뒤 아래 경로를 수정하거나 DATA_DIR 환경변수로 지정하세요.
DATA_DIR = os.environ.get("DNA_DATA_DIR", str(_BASE_DIR))

DATASETS = [
    {
        "name": "Jikji6",
        "r1": os.path.join(DATA_DIR, "4_S2_L001_R1_001.fastq"),
        "r2": os.path.join(DATA_DIR, "4_S2_L001_R2_001.fastq"),
        "measured_total_s": 10.1,   # CONTENT.md 실측치 (벡터화)
        "oligo_count": None,        # TODO: 확인 필요
    },
    {
        "name": "Jikji7",
        "r1": None,   # <경로 미확정 — NAS 참조 필요>
        "r2": None,
        "measured_total_s": 17.5,   # CONTENT.md 실측치 (벡터화)
        "oligo_count": None,
    },
    {
        "name": "Tripitaka",
        "r1": None,   # <원본 FASTQ 로컬 없음>
        "r2": None,
        "measured_total_s": None,   # 아직 미측정 — extrapolation 대상
        "oligo_count": 92626,       # CONTENT.md 명시
    },
]

# 측정할 단계 이름 (Figure1 baseline Timecell 키와 대응)
STAGES = [
    "Evert",
    "Primer Sorting",
    "Key Finding",
    "Index Sorting",
    "RS Decoding",
    "Null Finding",
    "XOR Recovery",
    "STL Decode",
]


# =============================================================================
# 헬퍼: FASTQ 파일 존재 여부 확인
# =============================================================================
def check_data_availability(datasets):
    print("\n[데이터 가용성 확인]")
    for ds in datasets:
        r1_ok = os.path.isfile(ds["r1"]) if ds["r1"] else False
        r2_ok = os.path.isfile(ds["r2"]) if ds["r2"] else False
        status = "OK" if (r1_ok and r2_ok) else "MISSING"
        r1_disp = ds["r1"] if ds["r1"] else "<미확정>"
        r2_disp = ds["r2"] if ds["r2"] else "<미확정>"
        print(f"  {ds['name']:12s} [{status}]  R1={r1_disp}  R2={r2_disp}")
    print()


# =============================================================================
# 실제 파이프라인 실행 (baseline import 성공 시)
# =============================================================================
def run_pipeline_real(r1_path, r2_path):
    """
    Figure1 baseline 함수를 단계별로 호출하며 Timecell 을 채움.
    baseline 이 단일 main() 이 아닌 분리된 함수로 제공될 때만 작동.

    TODO: dna_decode_baseline.py 에 아래 시그니처의 함수를 추가해야 함:
          def run_pipeline(r1_path, r2_path) -> dict  (Timecell 반환)
    """
    Timecell = {s: None for s in STAGES}

    seq = read_fastq_file(r1_path) + read_fastq_file(r2_path)
    reads_count = len(seq)

    t = time.time()
    seq = reverse_evert_sequences(seq)
    Timecell["Evert"] = round(time.time() - t, 3)

    t = time.time()
    sorted_reads, primer_ids = primer_sort(seq)
    Timecell["Primer Sorting"] = round(time.time() - t, 3)

    t = time.time()
    key_data = key_find(sorted_reads)
    Timecell["Key Finding"] = round(time.time() - t, 3)

    t = time.time()
    consensus_data = consensus_call(key_data)
    Timecell["Index Sorting"] = round(time.time() - t, 3)

    t = time.time()
    rs_data = rs_decode(consensus_data)
    Timecell["RS Decoding"] = round(time.time() - t, 3)

    t = time.time()
    null_data = null_find(rs_data)
    Timecell["Null Finding"] = round(time.time() - t, 3)

    t = time.time()
    xor_data = xor_recovery(null_data)
    Timecell["XOR Recovery"] = round(time.time() - t, 3)

    Timecell["STL Decode"] = None   # STL 재구성은 volume scaling 측정 범위 밖
    Timecell["_reads_count"] = reads_count
    return Timecell


# =============================================================================
# stub: 실측치 기반 Timecell (baseline import 불가 시 CONTENT.md 값 사용)
# =============================================================================
def run_pipeline_stub(dataset_name, measured_total_s):
    """
    실제 실행 없이 CONTENT.md 에 기록된 총 시간만 반환.
    단계별 비율은 측정된 분석 결과(CONTENT.md, 02_Optimization_Analysis.md)에서 인용.

    NOTE: 이 값들은 실측치가 아닌 참조 논문/분석 문서에서 추출한 추정치임.
    단계 비율 출처: 260328_0145_Sec2_Final_Optimization.md
      - Evert          ~2.4%
      - Primer Sorting ~32.5%  (→ 최적화 후 20.7% 로 감소)
      - Key Finding    ~32.5%  (→ 최적화 후 메인 병목)
      - RS Decoding    ~20.7%  (최적화 후 병목)
      - Null Finding   ~5%
      - XOR Recovery   ~6.9%
    """
    if measured_total_s is None:
        return None   # Tripitaka 는 미측정
    return {
        "Evert":           round(measured_total_s * 0.024, 3),
        "Primer Sorting":  round(measured_total_s * 0.207, 3),
        "Key Finding":     round(measured_total_s * 0.325, 3),
        "Index Sorting":   round(measured_total_s * 0.050, 3),
        "RS Decoding":     round(measured_total_s * 0.207, 3),
        "Null Finding":    round(measured_total_s * 0.050, 3),
        "XOR Recovery":    round(measured_total_s * 0.069, 3),
        "STL Decode":      None,
        "_reads_count":    None,
        "_is_stub":        True,
    }


# =============================================================================
# 결과 저장
# =============================================================================
def save_results_csv(results, output_path):
    fieldnames = ["dataset"] + STAGES + ["total_s", "reads_count", "is_stub"]
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    print(f"[결과 저장] {output_path}")


# =============================================================================
# 메인 루프 — 데이터셋별 순차 실행
# =============================================================================
def main():
    print("=" * 60)
    print("Figure 2 — Data Volume Scaling Experiment")
    print("=" * 60)

    check_data_availability(DATASETS)

    results = []

    for ds in DATASETS:
        name = ds["name"]
        r1, r2 = ds["r1"], ds["r2"]
        measured = ds["measured_total_s"]

        print(f"\n[{name}] 처리 시작...")

        r1_exists = r1 and os.path.isfile(r1)
        r2_exists = r2 and os.path.isfile(r2)

        if _BASELINE_AVAILABLE and r1_exists and r2_exists:
            # 실제 파이프라인 실행
            print(f"  -> 실제 파이프라인 실행 (R1={r1})")
            try:
                Timecell = run_pipeline_real(r1, r2)
                is_stub = False
            except Exception as e:
                print(f"  [오류] {e} — stub 으로 대체")
                Timecell = run_pipeline_stub(name, measured)
                is_stub = True
        else:
            # Stub 모드: CONTENT.md 실측치 사용
            reason = "baseline import 실패" if not _BASELINE_AVAILABLE else "FASTQ 없음"
            print(f"  -> stub 모드 ({reason}). CONTENT.md 실측치 사용.")
            Timecell = run_pipeline_stub(name, measured)
            is_stub = True

        if Timecell is None:
            print(f"  -> {name}: 미측정 (extrapolation 대상) — 스킵")
            row = {"dataset": name}
            for s in STAGES:
                row[s] = "N/A (미측정)"
            row["total_s"] = "N/A"
            row["reads_count"] = "N/A"
            row["is_stub"] = "extrapolation"
        else:
            stage_sum = sum(v for v in [Timecell.get(s) for s in STAGES] if v is not None)
            reads = Timecell.get("_reads_count", "N/A")
            row = {"dataset": name}
            for s in STAGES:
                row[s] = Timecell.get(s, "N/A")
            row["total_s"] = round(stage_sum, 3) if not is_stub else measured
            row["reads_count"] = reads
            row["is_stub"] = is_stub
            print(f"  총 시간: {row['total_s']} s  (stub={is_stub})")

        results.append(row)

    # CSV 저장
    output_csv = os.path.join(os.path.dirname(__file__), "volume_scaling_results.csv")
    save_results_csv(results, output_csv)

    # 서브선형 확장성 요약 출력
    print("\n[확장성 요약] CONTENT.md 실측치 기준:")
    print("  J6 → J7: 3.1× 데이터 증가 → 1.7× 시간 증가 (sub-linear)")
    print("  실측: J6=10.1s, J7=17.5s (벡터화 파이프라인)")
    print("  Tripitaka 실측 미완료 — 다음 단계 TODO")

    print("\n완료. 결과:", output_csv)


if __name__ == "__main__":
    main()
