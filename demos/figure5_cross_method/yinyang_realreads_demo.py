#!/usr/bin/env python3
"""
Yin-Yang Code (Ping 2022) — iSCAN 삽입 *실 read* 검증 (시뮬 아님).

데이터(실측, CNSA CNP0001650 in-vitro oligo pool, File 3 = 도덕경):
  reads : Manuscript_Drafts_v2/refs/fastq/DS_Sep_S2_F3_2_1.fq.gz  (DNBSEQ/DIPSEQ PE100 R1)
  논문 oligo 구조: 200nt = 20nt 프라이머 + 160nt payload(128 data + 16 index + 16 RS) + 20nt 프라이머

핵심: 실 read에서 20nt 프라이머가 96.5% 검출됨(P1 fwd / P2 rev). PrimerSort로 방향정렬 후
프라이머 직후 고정구간(위치-인덱스 프록시)으로 그룹핑 → Antkowiak과 같은 원리.
(프라이머는 read에서 경험적으로 발견 — 모든 oligo 공통이라 압도적 최빈 5' 20-mer.)

비교 (read 묶기 단계만):
  (A) content-LSH(k-mer MinHash, payload 40nt)
  (B) iSCAN: PrimerSort(20nt 프라이머 검색) → 프라이머직후 30nt → 그룹핑

실행: python3 yinyang_realreads_demo.py [N_READS]   (기본 1,000,000)
"""
import os, sys, subprocess, time, itertools
import numpy as np
np.random.seed(5)

FQ = os.path.join(os.environ.get("ISCAN_DEMO_DATA", "data"), "yinyang", "DS_Sep_S2_F3_2_1.fq.gz")
N  = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
P1 = "ATTCGATATTAACTAATGAG"   # forward 20nt primer (read에서 발견)
P2 = "ATGCTACGCTTGACCTATTA"   # reverse 20nt primer
F3_OLIGOS = 1716              # File 3(도덕경) 추정 oligo 수 (≥20read 안정그룹 = 파일크기비례 ~17% × 10,103)

p = subprocess.run("gzcat '%s' | head -%d | awk 'NR%%4==2'" % (FQ, N*4),
                   shell=True, capture_output=True, text=True)
reads = [r for r in p.stdout.split("\n") if len(r) >= 90][:N]
print("실 YinYang read %d (F3, PE100 R1, DNBSEQ)" % len(reads))

def mm(a,b): return sum(1 for x,y in zip(a,b) if x!=y)

# ===== (B) iSCAN: PrimerSort → 프라이머직후 위치읽기 =====
t0=time.time(); groups={}; det=0; payloads=[]
for r in reads:
    if mm(r[:20],P1) <= 2:
        det+=1; groups.setdefault(r[20:50],[]).append(1); payloads.append(r[20:60])
    elif mm(r[:20],P2) <= 2:
        det+=1
t_is=time.time()-t0
big5  = sum(1 for v in groups.values() if len(v)>=5)
big20 = sum(1 for v in groups.values() if len(v)>=20)
print("\n[read 묶기 단계 — 실 read]")
print("  (B) iSCAN PrimerSort+위치읽기: %.2fs" % t_is)
print("      프라이머 검출 %d/%d (%.1f%%)  (P1+P2 양끝 20nt)" % (det, len(reads), 100*det/len(reads)))
print("      P1 그룹: 고유 %d, >=5read %d, >=20read %d (안정 oligo)" % (len(groups), big5, big20))
print("      oligo 복원: >=20read %d / File3 ~%d (≈%.0f%%)" % (big20, F3_OLIGOS, 100*big20/F3_OLIGOS))

# ===== (A) content-LSH (payload 40nt) =====
def kmerDNA(s,k=5):
    out=[]
    for e in range(len(s)-k+1):
        idx=0;ok=True
        for j,c in enumerate(s[e:e+k]):
            v={'A':0,'C':1,'G':2,'T':3}.get(c)
            if v is None: ok=False;break
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
SUB=100000
sub=payloads[:SUB]
t0=time.time(); cl=lsh(sub); t_lsh=time.time()-t0
t0=time.time(); g2={}
for q in sub: g2.setdefault(q[:30],[]).append(1)
t_is2=time.time()-t0
print("\n  (A) content-LSH (payload 40nt, %d reads): %.1fs" % (len(sub), t_lsh))
print("      (동일 %d서 iSCAN 위치읽기: %.3fs)" % (len(sub), t_is2))
print("  → 묶기 속도: iSCAN이 content-LSH 대비 ~%.0f배 (PrimerSort 포함 동일스케일 환산)"
      % (t_lsh/max(t_is/(len(reads)/SUB),1e-6)))
print("\n참고: 묶기 단계만 측정. 프라이머는 read에서 경험적 발견(P1/P2). 시뮬 아님 — 실 DNBSEQ read.")
