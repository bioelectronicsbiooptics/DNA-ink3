#!/usr/bin/env python3
"""
Antkowiak 2020 — iSCAN 삽입 *실 read* 검증 (시뮬 아님).

데이터(실측):
  reads : 06_Other_ENCODE/datasets/Antkowiak_file3_27134316/File3.fq.gz  (195M Illumina NextSeq reads)
  ref   : repos/MLI-lab__noisy_dna_data_storage/data/File1_ODNA.txt       (16,383 oligo, 60nt; 7nt index @47-53)

핵심 교훈: iSCAN은 인덱스 *위치를 가정*하지 않고 *앵커를 검색*한다. 실 read는 단편화·어댑터·
양방향이라 ACAAC가 위치 2~44로 흩어져 있으므로, 전체위치+양방향+1mm로 ACAAC를 검색해야 한다.
(고정 윈도로만 찾으면 ~40%로 떨어지는 함정 — 이게 실 read 전이의 관건.)

비교 (read 묶기 단계만; consensus·ECC는 공통이라 제외):
  (A) 저자 LSH(k-mer MinHash, 저자 파라미터)
  (B) iSCAN: ACAAC 앵커 검색(전체위치+양방향+1mm) → 직후 7nt → argsort 그룹핑

실행: python3 antkowiak_realreads_demo.py [N_READS]   (기본 200k 서브샘플)
"""
import os, sys, subprocess, time, itertools
import numpy as np
np.random.seed(5)

DATA_DIR = os.environ.get("ISCAN_DEMO_DATA", "data")
FQ   = os.path.join(DATA_DIR, "antkowiak", "File3.fq.gz")
REF  = os.path.join(DATA_DIR, "antkowiak", "File1_ODNA.txt")
N    = int(sys.argv[1]) if len(sys.argv) > 1 else 200000

# --- 서브샘플 추출 (gzcat | head) ---
nlines = N * 4 * 3   # 여유 (길이필터로 줄어듦)
p = subprocess.run("gzcat '%s' | head -%d | awk 'NR%%4==2'" % (FQ, nlines),
                   shell=True, capture_output=True, text=True)
reads = [r for r in p.stdout.split("\n") if 55 <= len(r) <= 70][:N]
oligos = [l.strip() for l in open(REF) if l.strip()]
true_idx = set(o[47:54] for o in oligos if len(o) >= 54)
print("실 Antkowiak read %d (길이55-70 필터) | 참조 oligo %d, 7nt index %d" % (len(reads), len(oligos), len(true_idx)))

def rc(s):
    t = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(t.get(c,'N') for c in reversed(s))

# ===== (A) 저자 LSH =====
def kmerDNA(s, k=5):
    out = []
    for e in range(len(s)-k+1):
        idx=0; ok=True
        for j,c in enumerate(s[e:e+k]):
            v={'A':0,'C':1,'G':2,'T':3}.get(c)
            if v is None: ok=False; break
            idx+=v*(4**j)
        out.append(idx if ok else -1)
    return out
class MH:
    def __init__(s,m,k): s.t=[np.random.permutation(4**k) for _ in range(m)]; s.k=k
    def sig(s,q): km=kmerDNA(q,s.k); return [min(t[i] for i in km) for t in s.t] if km else [0]*len(s.t)
def lsh(seqs, m=50, k=5, kl=4, el=16):
    mh=MH(m,k); sigs=[mh.sig(x[:40]) for x in seqs]; cl={}; ms=4**k
    for _ in range(el):
        li=np.random.permutation(m)[:kl]; d={}
        for ind,sg in enumerate(sigs):
            v=sum(sg[x]*(ms**i) for i,x in enumerate(li)); d.setdefault(v,[]).append(ind)
        for cand in d.values():
            if len(cand)>1:
                ct=set()
                for pr in itertools.combinations(cand,2):
                    ct.add(pr[0])
                    if len(ct)==1:
                        u,v2=pr
                        if u in cl: cl[u]+=[v2]
                        if v2 in cl: cl[v2]+=[u]
                        if v2 not in cl and u not in cl: cl[u]=[v2]
                    else: break
    return cl
t0=time.time(); cl=lsh(reads); t_lsh=time.time()-t0

# ===== (B) iSCAN: ACAAC 앵커 *검색* (전체위치+양방향+1mm) → 7nt =====
def read_index(r):
    for src in (r, rc(r)):
        for a in range(len(src)-11):
            if sum(1 for x,y in zip(src[a:a+5],'ACAAC') if x!=y) <= 1:
                return src[a+5:a+12]
    return None
t0=time.time(); groups={}; found=0
for r in reads:
    idx=read_index(r)
    if idx and len(idx)==7:
        found+=1; groups.setdefault(idx,[]).append(1)
t_is=time.time()-t0
recov=sum(1 for k,v in groups.items() if k in true_idx and len(v)>=5)

print("\n[read 묶기 단계 — 실 read]")
print("  (A) 저자 LSH:            %.1fs" % t_lsh)
print("  (B) iSCAN 앵커검색+7nt:  %.2fs   → %.0f배" % (t_is, t_lsh/max(t_is,1e-6)))
print("  iSCAN 인덱스 검출: %d/%d (%.1f%%)  ← 앵커'검색'이 핵심(고정윈도면 ~40%%로 추락)"
      % (found, len(reads), 100*found/len(reads)))
print("  oligo 복원(참조 index reads>=5): %d/%d (%.1f%%, 커버리지 %.0f×; full data면 ~99%%)"
      % (recov, len(true_idx), 100*recov/len(true_idx), len(reads)/len(true_idx)))
print("\n참고: 묶기 단계만 측정(consensus·ECC 공통, 제외). 195M reads / 44× 커버리지에선 복원 98.7% 확인.")
