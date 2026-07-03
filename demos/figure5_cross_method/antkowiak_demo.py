#!/usr/bin/env python3
"""
Antkowiak 2020 (Nat Commun) — iSCAN 삽입 재현 데모 (가장 강한 케이스)

근거(실측):
  - 실 레퍼런스 oligo 16,383개 (File1_ODNA.txt), 전부 60 nt
  - 보존모티프 ACAAC @ pos42-46 = 100% 보존(16,383/16,383)
  - 7nt 인덱스 @ pos47-53 = 16,383개 완전 유일 (= 논문 N=16,383, S=120 bits와 정합)
  - read 묶기 baseline = 저자(MLI-lab) LSH(k-mer MinHash) 코드 그대로

두 단계 비교:
  [구간1] read 묶기:  (A) 저자 LSH 클러스터링  vs  (B) iSCAN(ACAAC앵커→7nt→argsort)
  [구간2] consensus:  각 방식의 '자기 그룹(>=5 reads)'으로 MUSCLE MSA → 공정 end-to-end

실행: python3 antkowiak_demo.py
"""
import numpy as np, time, itertools, random, subprocess, collections, os
random.seed(5); np.random.seed(5)

DATA_DIR = os.environ.get("ISCAN_DEMO_DATA", "data")
ODNA = os.path.join(DATA_DIR, "antkowiak", "File1_ODNA.txt")
MUSCLE = os.environ.get("MUSCLE", "muscle")

oligos=[l.strip() for l in open(ODNA) if l.strip()]
# 데이터 사실 자가검증
anc100 = all(o[42:47]=='ACAAC' for o in oligos if len(o)>=47)
uniq7  = len({o[47:54] for o in oligos if len(o)>=54})
print("실 oligo: %d개, 60nt | ACAAC@42-46 100%%보존=%s | 7nt@47-53 유일값=%d (4^7=16384)"
      % (len(oligos), anc100, uniq7))

N_OLIGO=500; COPIES=30
sub=oligos[:N_OLIGO]
true_idx={o[47:54]:i for i,o in enumerate(sub)}
B='ACGT'; ps,pi,pd=0.02,0.005,0.005   # ~3% (sub2%+ins0.5%+del0.5%)
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
print("시뮬 reads: %d (oligo %d x %d copies, ~3%% error)\n" % (len(reads),N_OLIGO,COPIES))

# ---- (A) 저자 LSH(k-mer MinHash) ----
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
# LSH 파라미터 = 저자 노트북 cell11 실제값과 동일: m,k=50,5 / k_lsh=4 / ell_lsh=int(1/0.5**4)=16
t0=time.time(); cl=lsh(reads,50,5,4,16); t_lsh=time.time()-t0
lsh_groups=[sorted(set(m)|{c}) for c,m in cl.items()]
lsh_big=[g for g in lsh_groups if len(g)>=5]

# ---- (B) iSCAN: ACAAC앵커 → 7nt → argsort 그룹핑 ----
def fa(r):
    best=None                              # 정확매치(mm==0) 우선, 없으면 첫 mm==1
    for a in range(38,47):
        seg=r[a:a+5]
        if len(seg)<5: continue
        mm=sum(1 for x,y in zip(seg,'ACAAC') if x!=y)
        if mm==0: return a
        if mm==1 and best is None: best=a
    return best
t0=time.time(); groups={}; found=0; correct=0
for ri,r in enumerate(reads):
    a=fa(r)
    if a is None: continue
    idx=r[a+5:a+12]
    if len(idx)==7:
        found+=1; groups.setdefault(idx,[]).append(ri)
        if true_idx.get(idx,-1)==src[ri]: correct+=1
t_iscan=time.time()-t0
iscan_big=[g for g in groups.values() if len(g)>=5]

print("[구간1 read묶기]")
print("  (A) 저자 LSH:  %.2fs, 클러스터 %d개(>=5read %d개)" % (t_lsh,len(lsh_groups),len(lsh_big)))
print("  (B) iSCAN:     %.4fs, 그룹 %d개(>=5read %d개)  앵커검출 %.1f%%, per-read정확 %.1f%%"
      % (t_iscan,len(groups),len(iscan_big),100*found/len(reads),100*correct/len(reads)))
print("  → 묶기 속도: iSCAN이 LSH 대비 %.0f배" % (t_lsh/t_iscan))

# oligo 복원 지표
recov=pure=0
for oi,o in enumerate(sub):
    g=groups.get(o[47:54],[])
    if len(g)>=5:
        recov+=1
        if sum(1 for ri in g if src[ri]==oi)/len(g)>=0.8: pure+=1
print("  → oligo 복원(reads>=5): %d/%d (%.1f%%), purity>=80%%: %d/%d (%.1f%%)"
      % (recov,N_OLIGO,100*recov/N_OLIGO,pure,N_OLIGO,100*pure/N_OLIGO))

# ---- [구간2] 공정 consensus: 각 방식의 자기 >=5 그룹으로 MUSCLE MSA ----
def consensus_one(ids):
    with open('/tmp/_ac.fa','w') as f:
        for j,ri in enumerate(ids[:25]): f.write(">s%d\n%s\n"%(j,reads[ri]))
    subprocess.run([MUSCLE,'-align','/tmp/_ac.fa','-output','/tmp/_aco.fa'],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
def time_cons(grouplist, cap):
    gl=grouplist[:cap]
    t0=time.time()
    for g in gl: consensus_one(g)
    return time.time()-t0, len(gl)
# 공정 비교는 과분할 처리 가정에 민감 → 두 시나리오 모두 측정해 '범위'로 보고
# 시나리오 X: 각 방식 자기 >=5그룹 '전부' consensus (LSH 과분할이 그대로 비용)
tcX_lsh,nX_lsh=time_cons(lsh_big,10**9)
tcX_is, nX_is =time_cons(iscan_big,10**9)
# 시나리오 Y: consensus를 진짜 oligo수(500)로 cap (다운스트림 병합/dedup이 과분할 흡수)
tcY_lsh,nY_lsh=time_cons(lsh_big,N_OLIGO)
tcY_is, nY_is =time_cons(iscan_big,N_OLIGO)
print("\n[구간2 consensus] 각 방식 자기 >=5그룹으로 MUSCLE MSA")
print("  시나리오X(과분할 전부 consensus):  LSH %.2fs(%d MSA) | iSCAN %.2fs(%d MSA)"
      % (tcX_lsh,nX_lsh,tcX_is,nX_is))
print("  시나리오Y(병합후 500 cap):         LSH %.2fs(%d MSA) | iSCAN %.2fs(%d MSA)"
      % (tcY_lsh,nY_lsh,tcY_is,nY_is))
print("\n[end-to-end 구간1+2]")
print("  시나리오X: LSH %.2fs vs iSCAN %.2fs → %.1f배"
      % (t_lsh+tcX_lsh, t_iscan+tcX_is, (t_lsh+tcX_lsh)/(t_iscan+tcX_is)))
print("  시나리오Y: LSH %.2fs vs iSCAN %.2fs → %.1f배"
      % (t_lsh+tcY_lsh, t_iscan+tcY_is, (t_lsh+tcY_lsh)/(t_iscan+tcY_is)))
print("\n  ★견고한 핵심수치 = [구간1] 묶기 속도 %.0f배 (재현 안정)." % (t_lsh/t_iscan))
print("   end-to-end는 LSH 과분할(>=5그룹 %d개 vs 진짜 %d) 처리 가정에 따라 X~Y로 변동."
      % (len(lsh_big),N_OLIGO))
print("   본 데모는 LSH만(저자 파라미터 동일); 저자 풀 파이프라인은 LSH 뒤 SSW 클러스터")
print("   정제(local_pairwise_align_ssw, match>=35)가 있어 과분할이 줄어 시나리오Y에 가까움.")
