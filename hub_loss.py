"""Hub choice by loss: for each hub, over 200 yr (daily samples), best-delay open route from every
settlement; report mean/worst route first-try success and route abandonment P(any hop loses all 4)."""
import numpy as np, os, json
import mpx
OUT=os.path.join(mpx.HERE,"out"); SC=os.path.join(OUT,"scan")
cl=json.load(open(os.path.join(OUT,"e4_closures.json"))); T=cl["T_end_h"]
grid=np.arange(0.0,T+1,1.0); idx=np.arange(0,len(grid),24)
Ld={}
for (s,r) in mpx.DIRECTED:
    z=np.load(os.path.join(SC,f"{s}__{r}.npz".replace(" ","_")))
    Ld[(s,r)]=(z["clr"][idx].astype(float),z["tau"][idx].astype(float),z["d"][idx].astype(float))
g=grid[idx]
def rt(route):
    ready=g.copy(); ok=np.ones(len(g),bool); ps=[]
    for k in range(len(route)-1):
        c,t,d=Ld[(route[k],route[k+1])]; te=ready+mpx.SEC
        ok&=np.interp(te,g,c)>=0.1; ta=te+np.interp(te,g,t); ps.append(mpx.loss_backbone(np.interp(te,g,d)))
        ready=ta+mpx.SEC if route[k+1] in mpx.RELAYS else ta
    ps=np.array(ps); return ready-g, ok, np.prod(1-ps,0), 1-np.prod(1-ps**4,0)
rows=[]
for H in mpx.SETTLEMENTS:
    fs,ab=[],[]
    for X in mpx.SETTLEMENTS:
        if X==H: continue
        for a,b in ((X,H),(H,X)):
            bd=np.full(len(g),np.inf); bf=np.zeros(len(g)); bab=np.ones(len(g))
            for r in mpx.candidate_routes(a,b):
                d,ok,f,A=rt(r); m=ok&(d<bd); bd=np.where(m,d,bd); bf=np.where(m,f,bf); bab=np.where(m,A,bab)
            fs.append(bf); ab.append(bab)
    fs=np.array(fs); ab=np.array(ab)
    rows.append((H,fs.mean(),fs.min(),ab.mean(),ab.max()))
rows.sort(key=lambda r:-r[1])
print("| Hub | Mean first-try route success | Worst first-try | Mean P(route abandoned) | Worst P(route abandoned) |")
print("|---|---|---|---|---|")
for H,a,b,c,d in rows: print(f"| {H} | {a:.4f} | {b:.4f} | {c:.2e} | {d:.2e} |")
