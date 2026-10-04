import json, copy, mpx, scenarios as S
def with_maint(extra=(), drop_bn=False):
    saved=copy.deepcopy(mpx.MAINT)
    if drop_bn: mpx.MAINT.pop(frozenset(("Relay B","Neptune")),None)
    for e,w in extra: mpx.MAINT.setdefault(frozenset(e),[]).append(w)
    return saved
def incident_links(t0,t1): return [((x,"Relay A"),(t0,t1)) for x in mpx.SETTLEMENTS]+[(("Relay A","Relay B"),(t0,t1))]
def mincoll(s):
    m={};best={}
    for t,acct,d,_,bal in sorted(s.events,key=lambda e:e[0]):
        if acct.startswith("margin:"): best[acct]=min(best.get(acct,1e9),bal)
    return best
def pack(s,r): 
    r=dict(r); r["min_coll"]=mincoll(s); r["msgs"]=[{k:m[k] for k in("kind","a","b","t","arr","copies","orig","launch","route")} for m in s.msgs]; return r
out={}
# S2: falling run (Neptune long) + forced loss at Relay A hours 0-6
for tag,drop,win in (("S2_maint",False,(0,6)),("S2_nomaint",True,(0,6)),("S2_maint_12h",False,(0,12))):
    saved=with_maint(incident_links(*win),drop); s,r=S.price_run(0.0,-1); out[tag]=pack(s,r); mpx.MAINT.clear(); mpx.MAINT.update(saved)
    print(tag,"complete",round(r["complete"],2),"calls",r["calls"],"peak shortfall",round(r["peak_exposure"],3),"min coll",{k:round(v,3) for k,v in out[tag]["min_coll"].items()})
s,r=S.price_run(0.0,-1); out["S2_base"]=pack(s,r)
# S3: value move to the hub from best, median, worst settlement
for x in ("Ceres","Jupiter","Neptune"):
    s,r=S.value_move(0.0,src=x,dst="Earth"); out[f"S3_{x}"]=pack(s,r); print("S3",x,"credit",round(r["credit"],3),"final",round(r["complete"],3))
# E5: all four scenarios at later epochs, balances reset
for lab,t in (("+1y",8766.0),("+10y",87660.0),("+100y",876600.0),("difficult",650208.0)):
    for name,fn in (("value",lambda t=t:S.value_move(t)),("long",lambda t=t:S.long_lived(t)),("rise",lambda t=t:S.price_run(t,+1)),("fall",lambda t=t:S.price_run(t,-1))):
        s,r=fn(); out[f"E5_{lab}_{name}"]=pack(s,r)
        print("E5",lab,name,"complete",round(r["complete"],2),"orig",r["orig"],"peak enc",round(r["peak_enc"],2),"asset-h",round(r["asset_hours"],1),"conserved",r["conserved"],"repins",r.get("repins","-"),"calls",r.get("calls","-"),"shortfall",round(r.get("peak_exposure",0),2))
json.dump(out,open("out/runs_all.json","w"),default=float)
