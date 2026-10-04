import json, numpy as np, importlib, e5_runs as E
LQ=[0.15,0.20,0.25,0.30,0.35,0.40]; TR=[0.05,0.10,0.15,0.20]; out={}
for d in ("rising","falling"):
    M=np.zeros((len(TR),len(LQ)))
    for i,tr in enumerate(TR):
        for j,lq in enumerate(LQ):
            r=E.run(tr if d=="rising" else -tr,"t3",True,n=600,liqf=lq,seed=11)
            M[i,j]=min(1.0,r["P_default"])
    out[d]=M.tolist(); print(d,np.round(M,3).tolist(),flush=True)
json.dump(dict(LQ=LQ,TR=TR,**out),open("out/phase.json","w"))
