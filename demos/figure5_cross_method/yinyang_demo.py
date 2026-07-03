#!/usr/bin/env python3
"""
Yin-Yang Code (YYC, Ping 2022) — iSCAN 삽입 재현 데모 (read 묶기 단계)

정직성 전제:
  - YYC의 인덱스는 'DNA 위치필드'가 아니라 디코딩된 binary의 선두 비트다
    (yyc/pipeline.py:130 index_operator.divide_all → 선두 index 비트 분리).
    support_spacing=0(기본)이므로 선두 DNA 염기가 선두 index 비트에 1:1 대응
    → 선두 K nt가 '위치-인덱스 프록시'로 작동(아래서 유일성 실측 검증).
  - YYC는 codec 내부에 read 묶기 단계가 없다(논문은 외부 SOAPdenovo 어셈블 사용,
    본 환경에서 구동 불가). 따라서 baseline은 '일반 content 클러스터링'을 대표하는
    LSH(k-mer MinHash, Antkowiak 저자 코드와 동일 알고리즘)로 둔다. → illustrative.

데이터(실측 oligo, 시뮬 read):
  encoded.dna : 3,251 실 YYC 출력 서열 (133 nt)  ← 진짜 YinYang 인코딩 결과
  reads       : 위 oligo의 노이즈 카피(sub/ins/del) 시뮬

실행: python3 yinyang_demo.py
"""
import numpy as np, time, itertools, random, os
random.seed(11); np.random.seed(11)

YY = os.path.join(os.environ.get("ISCAN_DEMO_DATA", "data"), "yinyang", "encoded.dna")
oligos = [l.strip() for l in open(YY) if l.strip()]
print("실 YYC oligo: %d개, 길이 %s" % (len(oligos), sorted(set(len(o) for o in oligos))))

# ---- 위치-인덱스 프록시(선두 K nt) 유일성 실측: 최소 K 자동선택 ----
N_OLIGO = 500; COPIES = 30
sub = oligos[:N_OLIGO]
K = None
for k in (12,16,20,24,28,32):
    pref = [o[:k] for o in sub]
    if len(set(pref)) == len(sub):
        K = k; break
if K is None:
    K = 40
print("선두 K=%d nt → %d/%d oligo 유일 (위치-인덱스 프록시 성립)"
      % (K, len({o[:K] for o in sub}), N_OLIGO))
true_pref = {o[:K]: i for i,o in enumerate(sub)}

# ---- 시뮬 reads (sub 2% / ins 0.5% / del 0.5% = DNBSEQ/Illumina 급 저에러) ----
B='ACGT'; ps,pi,pd = 0.02,0.005,0.005
reads=[]; src=[]
for i,o in enumerate(sub):
    for _ in range(COPIES):
        r=[]
        for ch in o:
            x=random.random()
            if x<pd: continue
            elif x<pd+pi: r.append(random.choice(B)); r.append(ch)
            elif x<pd+pi+ps: r.append(random.choice([b for b in B if b!=ch]))
            else: r.append(ch)
        reads.append(''.join(r)); src.append(i)
print("시뮬 reads: %d (oligo %d x %d copies, ~3%% error)" % (len(reads),N_OLIGO,COPIES))

# ================= (A) 일반 content 클러스터링 — LSH(k-mer MinHash, 저자 알고리즘) =================
def kmerDNA(s,k=3):
    out=[]
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
    def sig(s,q): km=kmerDNA(q,s.k); return [min(t[i] for i in km) for t in s.t]
def sp(sigs,m,kl,el,ms):
    p=set()
    for _ in range(el):
        li=np.random.permutation(m)[:kl]; ls=[]
        for sg in sigs:
            v=0
            for i,x in enumerate(li): v+=sg[x]*(ms**i)
            ls.append(v)
        d={}
        for ind,sg in enumerate(ls): d.setdefault(sg,[]).append(ind)
        for cand in d.values():
            ct=set()
            if len(cand)>1:
                for pr in itertools.combinations(cand,2):
                    ct.add(pr[0])
                    if len(ct)==1: p.add(pr)
                    else: break
        yield p; p=set()
def lsh(seqs,m,k,kl,el):
    mh=MH(m,k); sigs=[mh.sig(x[:40]) for x in seqs]; cl={}
    for out in sp(sigs,m,kl,el,4**k):
        for u,v in sorted(out):
            if u in cl: cl[u]+=[v]
            if v in cl: cl[v]+=[u]
            if v not in cl and u not in cl: cl[u]=[v]
    return cl
t0=time.time(); cl=lsh(reads,50,5,4,16); t_lsh=time.time()-t0
print("\n(A) 일반 content 클러스터링 LSH(k-mer MinHash): %.2fs, %d clusters" % (t_lsh,len(cl)))

# ================= (B) iSCAN 위치-인덱스(선두 K nt) 그룹핑 =================
t0=time.time()
groups={}; correct=0
for ri,r in enumerate(reads):
    pref=r[:K]                       # 선두 K nt = 위치-인덱스 프록시
    groups.setdefault(pref,[]).append(ri)
    if true_pref.get(pref,-1)==src[ri]: correct+=1
t_iscan=time.time()-t0
print("(B) iSCAN 위치-인덱스(선두 %dnt) 그룹핑: %.4fs, %d groups" % (K,t_iscan,len(groups)))
print("    정확 oligo 배정(per-read, exact prefix): %d/%d (%.1f%%)"
      % (correct,len(reads),100*correct/len(reads)))

# oligo별 복원 지표 (진짜 prefix 그룹에 reads>=5)
recov=0; pure=0
for oi,o in enumerate(sub):
    key=o[:K]; g=groups.get(key,[])
    if len(g)>=5:
        recov+=1
        same=sum(1 for ri in g if src[ri]==oi)
        if same/len(g)>=0.8: pure+=1
print("    복원가능 oligo(reads>=5): %d/%d (%.1f%%), 그중 purity>=80%%: %d/%d (%.1f%%)"
      % (recov,N_OLIGO,100*recov/N_OLIGO,pure,N_OLIGO,100*pure/N_OLIGO))

print("\n=== YinYang read 묶기 단계 비교 ===")
print("  일반 LSH 클러스터링: %.2fs" % t_lsh)
print("  iSCAN 위치읽기:      %.4fs" % t_iscan)
print("  → 속도: iSCAN이 LSH 대비 %.0f배" % (t_lsh/t_iscan))
print("  주의: YYC는 in-codec 묶기 부재 → baseline은 일반 content 클러스터(LSH) 대용.")
print("        실제 YinYang 논문의 외부 어셈블(SOAPdenovo)은 더 무거우므로 실격차는 더 큼.")
