"""E5 Tiers 2-3: rising and falling price runs through the full margin system. Writes out/e5_runs.json."""
import numpy as np, json, math
A=np.load("out/tiers_arrays.npz"); D=json.load(open("out/design.json"))
Bf=1.5; G=1.0; SIG=0.30; Z=2.62; KS=6.0; MM=0.8; LIQ=0.30
# book: settlement, signed notional (M credits), worst round trip (h)
BOOK=[("Bank A","Earth",+600,0.0),("Bank B","Mars",-400,None),("Operator C","Jupiter",-350,None),("Operator D","Neptune",+150,None)]
def RT(w): return 0.0 if w=="Earth" else D[w]["in_dmax"]+D[w]["out_dmax"]
MEM=[]
for n,w,N,_ in BOOK:
    M=RT(w)+G+Bf; f=SIG*math.sqrt(M/8760); MEM.append(dict(n=n,w=w,N=N,M=M,IM=Z*f*abs(N),S=(KS-Z)*f*abs(N)))
Sv=[m["S"] for m in MEM]; FUND=sum(sorted(Sv)[-2:]); SKIN=.25*FUND
for m,s in zip(MEM,Sv): m["share"]=FUND*s/sum(Sv)
# on-time probability per leg at the difficult epoch (year 74.17), adaptive copies
gi=int(round(650208/24)); 
def q(w):
    if w=="Earth": return 1.0
    f=float(A["F_"+w][gi]); k=int(A["K_"+w][gi]); return 1-(1-f)**k, k
QK={m["w"]:(q(m["w"]) if m["w"]!="Earth" else (1.0,0)) for m in MEM}
def run(trend, tail, incident, n=4000, H=72.0, seed=1, liqf=None):
    global LIQ
    if liqf is not None: LIQ=liqf
    rng=np.random.default_rng(seed); steps=int(H/Bf); dt=Bf/8760
    if tail=="gauss": eps=rng.standard_normal((n,steps))
    else: eps=rng.standard_t(3,(n,steps))/math.sqrt(3)
    mu=math.log(1+trend)/steps
    lr=mu+SIG*math.sqrt(dt)*eps-0.5*SIG**2*dt
    if tail=="t3": lr[:,20]+=math.log(1+math.copysign(0.10,trend))   # a 10% jump with the trend at hour 30
    P=np.exp(np.cumsum(lr,axis=1)); P=np.hstack([np.ones((n,1)),P])
    stats=dict(cause_missed=0,cause_liq=0,calls=0,calls_missed=0,defaults=0,hedges=0,tier=np.zeros(8),loss_beyond_own=[],orig_calls=[])
    for p in range(n):
        coll=[m["IM"] for m in MEM]; liq=[LIQ*abs(m["N"]) for m in MEM]; alive=[True]*4; pend=[None]*4; orig=0; wf=0.0; deepest=0
        for s in range(1,steps+1):
            t=s*Bf; px=P[p,s]; prev=P[p,s-1]
            for i,m in enumerate(MEM):
                if not alive[i]: continue
                coll[i]+=m["N"]*(px-prev)          # variation margin, settled inside the hub
                iso = incident and m["w"]=="Neptune" and 30<=t<102
                if pend[i] is not None and t>=pend[i][0]:   # call deadline passed
                    ok=pend[i][1]
                    if not ok and pend[i][2]==0 and not iso:      # undelivered call: one resend round before default
                        qq,k=QK[m["w"]]; orig+=2*k if m["w"]!="Earth" else 0
                        pend[i]=(t+RT(m["w"])+G, rng.random()<qq and rng.random()<qq, 1); continue
                    if not ok: stats["cause_missed"]+=0 if iso else 1
                    if ok:
                        need=m["IM"]-coll[i]
                        if need<=liq[i]: coll[i]+=need; liq[i]-=need; pend[i]=None
                        else: ok=False; stats["cause_liq"]+=1
                    if not ok:
                        alive[i]=False
                        loss=max(0.0,-(coll[i]+m["N"]*(P[p,min(s+1,steps)]-px)))   # close out within one fast batch
                        if iso: stats["hedges"]+=1
                        else: stats["defaults"]+=1
                        # waterfall beyond own margin: own fund share, skin, survivors, assessments
                        rem=loss; tiers=[m["share"],SKIN,FUND-m["share"],FUND-m["share"]]
                        for k_,cap in enumerate(tiers):
                            if rem<=0: break
                            use=min(rem,cap); rem-=use; deepest=max(deepest,k_+3)
                        if rem>0: deepest=7
                        stats["loss_beyond_own"].append(loss); continue
                if pend[i] is None and coll[i]<MM*m["IM"]:
                    stats["calls"]+=1; qq,k=QK[m["w"]]
                    if m["w"]!="Earth": orig+=2*k
                    delivered = (not iso) and rng.random()<qq and rng.random()<qq
                    if not delivered: stats["calls_missed"]+=1
                    pend[i]=(t+RT(m["w"])+G, delivered, 0)
        stats["tier"][deepest]+=1; stats["orig_calls"].append(orig/(H/24))
    out=dict(cause_missed=stats["cause_missed"]/n,cause_liq=stats["cause_liq"]/n,calls_per_run=stats["calls"]/n,missed=stats["calls_missed"]/max(1,stats["calls"]),P_default=stats["defaults"]/n,P_hedge=stats["hedges"]/n,
             mean_loss_beyond=float(np.mean(stats["loss_beyond_own"])) if stats["loss_beyond_own"] else 0.0,
             max_loss_beyond=float(np.max(stats["loss_beyond_own"])) if stats["loss_beyond_own"] else 0.0,
             tier_hist=(stats["tier"]/n).tolist(),orig_calls_day_mean=float(np.mean(stats["orig_calls"])),orig_calls_day_max=float(np.max(stats["orig_calls"])))
    return out
R={"members":MEM,"fund":FUND,"skin":SKIN,"qk":{k:list(v) if isinstance(v,tuple) else v for k,v in QK.items()}}
if __name__=="__main__":
  for tier,tail,inc,mag,lq in (("T2","gauss",False,0.05,0.30),("T3","t3",True,0.15,0.30),("T3_L20","t3",True,0.15,0.20),("T3_L40","t3",True,0.15,0.40)):
     for d in ("rising","falling"):
         tr=mag if d=="rising" else -mag
         r=run(tr,tail,inc,liqf=lq); R[f"{tier}_{d}"]=r
         print(tier,d,f"trend {tr:+.0%}",json.dumps({k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items() if k!='tier_hist'}),"deepest tier hist",[round(x,4) for x in r["tier_hist"]])
  for m in MEM: print(m)
  print("fund",FUND,"skin",SKIN,"QK",QK)
  json.dump(R,open("out/e5_runs.json","w"),indent=1,default=float)
