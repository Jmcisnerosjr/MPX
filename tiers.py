"""Adaptive copies, two-speed batch schedule, quota over 200 years. Writes out/tiers.json."""
import numpy as np, json, math, mpx
HUB="Earth"; S=[s for s in mpx.SETTLEMENTS if s!=HUB]
INNER=["Mercury","Venus","Mars","Ceres","Jupiter"]; OUTER=["Saturn","Uranus","Neptune"]
T=200*365.25*24; g=np.arange(0,T,24.0); TARGET=0.005
def best(a,b,t):
    bd=np.full(len(t),np.inf); bf=np.zeros(len(t))
    for rt in mpx.candidate_routes(a,b):
        ready=t.copy(); ok=np.ones(len(t),bool); D=np.zeros(len(t))
        for k in range(len(rt)-1):
            te=ready+mpx.SEC; st=mpx.link_state(rt[k],rt[k+1],te,maintenance=False)
            ok&=st["sun_ok"]; D+=st["d"]; ready=st["ta"]+(mpx.SEC if rt[k+1] in mpx.RELAYS else 0)
        d=ready-t; m=ok&(d<bd); bd=np.where(m,d,bd); bf=np.where(m,np.exp(-0.02*D),bf)
    return bd,bf
K={};F={};DL={};FL={}
for x in S:
    di,fi=best(x,HUB,g); do,fo=best(HUB,x,g); f=np.minimum(fi,fo)
    K[x]=np.ceil(np.log(TARGET)/np.log(1-f)).astype(int); F[x]=f; DL[x]=(di,do)
    FL[x]=mpx.light_time(x,HUB,g)[0]-g
out={}
for x in S: out[x]=dict(k_mean=float(K[x].mean()),k_max=int(K[x].max()),k_min=int(K[x].min()),k0=int(K[x][0]),floor_max=float(FL[x].max()))
for x in S: print(f"{x:8s} adaptive copies mean {K[x].mean():.2f} min {K[x].min()} max {K[x].max()} at h0 {K[x][0]} floor max {FL[x].max():.3f} h")
def daily(Bf,full_every,adaptive=True):
    nf=24/Bf; nfull=24/(Bf*full_every)
    kk={x:(K[x] if adaptive else np.full(len(g),K[x].max())) for x in S}
    inner=sum((kk[x]+1)*nf for x in INNER); outer=sum((kk[x]+1)*nfull for x in OUTER)
    return inner+outer
res={}
for Bf in (1.0,1.5,2.0):
    for fe in (3,4,6):
        if (Bf*fe)%1: continue
        a=daily(Bf,fe); w=daily(Bf,fe,False)
        res[f"{Bf}_{fe}"]=dict(Bf=Bf,full=Bf*fe,adapt_max=float(a.max()),adapt_mean=float(a.mean()),fixed=float(w.max()))
        print(f"fast every {Bf} h, full every {Bf*fe} h: adaptive worst day {a.max():.0f} mean {a.mean():.0f} | fixed worst-case {w.max():.0f}")
single=sum((K[x]+1)*8 for x in S); print("single 3h adaptive worst",single.max(),"mean",single.mean(),"fixed",sum((K[x].max()+1)*8 for x in S))
out["schedules"]=res; out["single3"]=dict(adapt_max=float(single.max()),adapt_mean=float(single.mean()))
# inner-tier clearing lag: max floor among inner
lag_in=np.max([FL[x] for x in INNER],axis=0); lag_all=np.max([FL[x] for x in S],axis=0)
out["lag_inner_max"]=float(lag_in.max()); out["lag_all_max"]=float(lag_all.max())
print("clearing lag inner max",lag_in.max(),"all max",lag_all.max())
json.dump(out,open("out/tiers.json","w"),indent=1)
np.savez_compressed("out/tiers_arrays.npz",g=g,**{f"K_{x}":K[x] for x in S},**{f"F_{x}":F[x] for x in S},**{f"FL_{x}":FL[x] for x in S},**{f"DI_{x}":DL[x][0] for x in S},**{f"DO_{x}":DL[x][1] for x in S})
