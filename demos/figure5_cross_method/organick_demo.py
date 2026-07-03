#!/usr/bin/env python3
"""
Organick/Rashtchian (Microsoft) — iSCAN 삽입 재현 데모 (read 묶기 단계)

데이터(실측 nanopore reads + 정답라벨, microsoft/clustered-nanopore-reads-dataset):
  Centers.txt  : 10,000 센터 서열(111 nt, 무작위)  ← '정답 oligo'
  Clusters.txt : 269,709 실 ONT MinION reads, 클러스터(=====구분)별로 묶임.
                 클러스터 i의 reads = 센터 i 의 노이즈 카피 → 정답 라벨 내장.
  (Organick 2018 random-access 라인의 클러스터링 벤치마크; 클러스터링 알고리즘 = Rashtchian NIPS2017.)

비교 (구간1 read 묶기):
  (A) content 클러스터링 baseline = LSH(k-mer MinHash, 저자 파라미터 m,k=50,5/k_lsh=4/ell_lsh=16)
  (B) iSCAN 위치-인덱스(선두 K nt) — 2변형:
      B1 exact   : 선두 K nt 그대로 그룹 (저에러용; O(N))
      B2 nearest : 선두 K nt를 센터 선두에 인델허용 최근접 (고에러 nanopore용; Yazdi식)

정직성: 이 벤치마크의 센터는 '무작위 111-mer'(에러정정 인덱스 미설계)라, nanopore 고에러에서
exact-prefix는 도전받음 → Organick이 🟡'에러내성 인덱스 필요' 군임을 실측으로 보여줌.

실행: python3 organick_demo.py
"""
import os
import numpy as np, time, itertools, random
random.seed(13); np.random.seed(13)

BASE = os.path.join(os.environ.get("ISCAN_DEMO_DATA", "data"), "organick")
N_CLUST=500          # 앞 500 클러스터(≈Antkowiak/YinYang 규모로 맞춤)
K=20                 # 위치-인덱스(선두) 윈도

centers=[l.strip() for l in open(BASE+"/Centers.txt")][:N_CLUST]
# Clusters.txt 파싱: '='로 시작하는 줄이 클러스터 구분자
reads=[]; src=[]
ci=-1
for l in open(BASE+"/Clusters.txt"):
    s=l.strip()
    if not s: continue
    if s.startswith("="):
        ci+=1
        if ci>=N_CLUST: break
        continue
    if 0<=ci<N_CLUST:
        reads.append(s); src.append(ci)
print("센터 %d개(111nt) | 실 nanopore reads %d개(앞 %d 클러스터, 평균 %.1f read/클러스터)"
      % (len(centers), len(reads), N_CLUST, len(reads)/N_CLUST))

# 위치-인덱스(선두 K) 유일성
pref_centers=[c[:K] for c in centers]
print("선두 %dnt 센터-프록시 유일성: %d/%d unique" % (K,len(set(pref_centers)),len(centers)))
true_pref={c[:K]:i for i,c in enumerate(centers)}

# 실 read 에러수준 가늠(선두 K가 센터와 정확히 일치하는 비율)
exact_pref_clean=sum(1 for r,c in zip(reads,src) if r[:K]==centers[c][:K])
print("read 선두 %dnt가 센터와 정확일치: %d/%d (%.1f%%) ← nanopore 에러/인델 수준 가늠"
      % (K,exact_pref_clean,len(reads),100*exact_pref_clean/len(reads)))

# ================= (A) LSH content 클러스터링 (저자 파라미터) =================
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
print("\n(A) LSH content 클러스터링: %.2fs, %d clusters" % (t_lsh,len(cl)))

# ================= (B1) iSCAN exact-prefix =================
t0=time.time(); g_exact={}; correct_e=0
for ri,r in enumerate(reads):
    p=r[:K]; g_exact.setdefault(p,[]).append(ri)
    if true_pref.get(p,-1)==src[ri]: correct_e+=1
t_exact=time.time()-t0

# ================= (B2) iSCAN nearest (인델허용, 센터 선두에 최근접) =================
def hamming(a,b): return sum(1 for x,y in zip(a,b) if x!=y)+abs(len(a)-len(b))
pc=[c[:K] for c in centers]
t0=time.time(); g_near={}; correct_n=0
for ri,r in enumerate(reads):
    pr=r[:K]; best=-1; bd=10**9
    for ti,tp in enumerate(pc):
        for off in (0,1,2,-1,-2):
            d=hamming(pr[off:off+K-2], tp[:K-2]) if off>=0 else hamming(pr[:K-2], tp[-off:-off+K-2])
            if d<bd: bd=d; best=ti
    g_near.setdefault(best,[]).append(ri)
    if best==src[ri]: correct_n+=1
t_near=time.time()-t0

def recov(groupmap, keyfn):
    rc=pu=0
    for oi,c in enumerate(centers):
        g=groupmap.get(keyfn(oi,c),[])
        if len(g)>=5:
            rc+=1
            if sum(1 for ri in g if src[ri]==oi)/len(g)>=0.8: pu+=1
    return rc,pu
re_e,pu_e=recov(g_exact, lambda oi,c:c[:K])
re_n,pu_n=recov(g_near,  lambda oi,c:oi)

print("(B1) iSCAN exact-prefix: %.4fs, %d groups, per-read정확 %.1f%%, oligo복원(>=5) %d/%d, purity %d/%d"
      % (t_exact,len(g_exact),100*correct_e/len(reads),re_e,N_CLUST,pu_e,N_CLUST))
print("(B2) iSCAN nearest(인델허용): %.2fs, %d groups, per-read정확 %.1f%%, oligo복원(>=5) %d/%d, purity %d/%d"
      % (t_near,len(g_near),100*correct_n/len(reads),re_n,N_CLUST,pu_n,N_CLUST))

print("\n=== Organick read 묶기 단계 비교 (실 nanopore reads) ===")
print("  (A) LSH content 클러스터링: %.2fs" % t_lsh)
print("  (B1) iSCAN exact-prefix:    %.4fs  → %.0f배 (단 고에러로 정확도↓)" % (t_exact,t_lsh/t_exact))
print("  (B2) iSCAN nearest:         %.2fs  → %.1f배 (정확도 회복)" % (t_near,t_lsh/t_near))
print("  해석: 무작위 센터(에러정정 인덱스 미설계)+nanopore 고에러 → exact는 빠르나 정확도 낮음,")
print("        nearest(인델허용)로 회복. Organick은 🟡'에러내성 인덱스 필요' 군임을 실측 확인.")
