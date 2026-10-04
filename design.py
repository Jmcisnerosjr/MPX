import numpy as np, json, math, mpx
HUB="Earth"; S=[s for s in mpx.SETTLEMENTS if s!=HUB]
T=200*365.25*24; g=np.arange(0,T,24.0)   # daily grid, 200 yr
def route_arr(rt,t0):
    ready=t0.copy(); ok=np.ones(len(t0),bool); D=np.zeros(len(t0))
    for k in range(len(rt)-1):
        te=ready+mpx.SEC; st=mpx.link_state(rt[k],rt[k+1],te,maintenance=False)
        ok&=st["sun_ok"]; D+=st["d"]; ready=st["ta"]+(mpx.SEC if rt[k+1] in mpx.RELAYS else 0)
    return ready-t0, ok, np.exp(-0.02*D)
res={}
for x in S:
    out={}
    for a,b,tag in ((x,HUB,"in"),(HUB,x,"out")):
        bd=np.full(len(g),np.inf); bf=np.zeros(len(g))
        for rt in mpx.candidate_routes(a,b):
            d,ok,f=route_arr(rt,g); m=ok&(d<bd); bd=np.where(m,d,bd); bf=np.where(m,f,bf)
        fl=mpx.LT_H*np.linalg.norm(mpx.pos(b,g+bd)-mpx.pos(a,g),axis=-1)  # straight-line floor (approx, receiver at arrival)
        # exact floor: solve light time a->b directly
        fl=mpx.light_time(a,b,g)[0]-g
        out[tag]=dict(dmax=float(bd.max()),fmin=float(bf.min()),fmean=float(bf.mean()),excess_max=float((bd-fl).max()),floor_max=float(fl.max()))
    fmin=min(out["in"]["fmin"],out["out"]["fmin"])
    k=math.ceil(math.log(0.005)/math.log(1-fmin))
    res[x]=dict(**{f"{t}_{kk}":v for t in out for kk,v in out[t].items()},k=k,ontime=1-(1-fmin)**k)
    r=res[x];print(f"{x:8s} worst first-try {fmin:.4f} copies {k} on-time {r['ontime']:.4f} | worst in {r['in_dmax']:.3f}h out {r['out_dmax']:.3f}h RT {r['in_dmax']+r['out_dmax']:.3f}h | max send-early {r['in_excess_max']*60:.1f} min floor max {r['in_floor_max']:.3f}h")
json.dump(res,open("out/design.json","w"),indent=1)
p=sum(res[x]["k"]+1 for x in S);print("originations per batch p =",p,"max batches/day",600/p,"min interval min",1440/(600/p))
