"""Operator local markets (basis bands), hash-locked cross-planet DvP, light-time expiries, aggregate receipts."""
import json, math, numpy as np, mpx, scenarios as S
HUB="Earth"; SIG=0.30; Z=2.62; FEE=0.0010; out={}
def path(t,d,step=0.005,n=8,B=1.5):
    k=min(n,max(0,int(math.floor(t/B)))); return (1+d*step)**k
# 1. basis: operator's local price = hub price as of one light-time ago
for d,lab in ((+1,"rising"),(-1,"falling")):
    for x in ("Mars","Jupiter","Saturn","Uranus","Neptune"):
        pk=0.0
        for t in np.arange(0,36,0.25):
            l=mpx.best_route(HUB,x,max(0,t-4.5))["delay_h"]; b=abs(path(t,d)-path(t-l,d))/path(t,d); pk=max(pk,b)
        l0=mpx.best_route(HUB,x,0.0)["delay_h"]; band=Z*SIG*math.sqrt(l0/8760)+FEE
        out[f"basis_{lab}_{x}"]=dict(peak=pk,band=band,inside=pk<=band)
        print(lab,x,f"peak basis {pk*100:.2f}% band {band*100:.2f}%",pk<=band)
# 2. hash-locked DvP: Earth client pays 100 credits for an asset held on Neptune's ledger
t0=0.0; s=S.Sim(t0)
a1=s.send("LOCK notice (Earth credits locked, hash H)",HUB,"Neptune",t0)
a2=s.send("LOCK notice (Neptune asset locked, hash H)","Neptune",HUB,a1)
a3=s.send("CLAIM asset with secret","Earth","Neptune",a2)          # secret revealed at Neptune
a4=s.send("CLAIM credits with secret","Neptune",HUB,a3)
TB=a3-t0+1.0                     # Neptune lock must outlast the claim's arrival, plus margin
TA=(a4-a1)+1.0                   # Earth lock must outlast Neptune's claim, measured from Neptune's lock
out["dvp"]=dict(lockA=t0,lockB_seen=a1,lockB_known=a2,claimB=a3,claimA=a4,complete=a4-t0,TB_min=TB,TA_rel=TA,orig=sum(m["orig"] for m in s.msgs),launch=sum(m["launch"] for m in s.msgs))
print("DvP",{k:round(v,3) for k,v in out["dvp"].items()})
# 3. expiry from light-time: E = min(30 d, 4 R_e), R_e = 2 T0 + 24 h
for x in mpx.SETTLEMENTS:
    if x==HUB: continue
    rt=mpx.best_route(x,HUB,0.0); Re=2*mpx.T0(rt["route"],0.0)+24; out[f"expiry_{x}"]=min(720.0,4*Re)
print("expiry h",{x:round(out[f'expiry_{x}'],1) for x in mpx.SETTLEMENTS if x!=HUB})
# 4. aggregate receipts: one custody signal per hop per batch instead of one per copy
T=json.load(open("out/tiers.json")); INNER=["Mercury","Venus","Mars","Ceres","Jupiter"]; OUTER=["Saturn","Uranus","Neptune"]
before=sum(T[x]["k_mean"]*2*16 for x in INNER)+sum(T[x]["k_mean"]*2*4 for x in OUTER)
after=sum(1*2*16 for x in INNER)+sum(1*2*4 for x in OUTER)
out["receipts"]=dict(before=before,after=after,saved=before-after)
print("order-copy receipt launches per day",round(before),"->",after,"saved",round(before-after))
json.dump(out,open("out/dvp_basis.json","w"),indent=1,default=float)
