import numpy as np, json, math, mpx
fh=mpx.fmt_h; HUB="Earth"; D=json.load(open("out/design.json")); B=3.0; G=1.0
def floor(a,b,t): return float(mpx.light_time(a,b,t)[0])-t
def first_try(rt): return float(np.prod([1-h["p_loss"] for h in rt["hops"]]))
out={}
# ---------------- S1: Neptune order through batch with cutoff Tb=6 h
Tb=6.0
fl=floor("Neptune",HUB,Tb); Dl=Tb+fl
rt0=mpx.best_route("Neptune",HUB,Tb-5,wait=False)   # route chosen near send time
delay=rt0["delay_h"]; k=D["Neptune"]["k"]
ts=Dl-delay-k*mpx.SEC
rt=mpx.route_timing(rt0["route"],ts,wait=False)
arr_last=rt["arrival"]+ (k-1)*mpx.SEC
maxfl=max(floor(x,HUB,Tb) for x in D)   # clearing waits for the farthest admitted floor
clear=Tb+maxfl
back=mpx.best_route(HUB,"Neptune",clear,wait=True)
ft=first_try(rt); ftb=first_try(back)
out["S1"]=dict(Tb=Tb,floor=fl,deadline=Dl,route=rt["route"],delay=delay,send=ts,arr=arr_last,clear=clear,back_route=back["route"],back_arr=back["arrival"],ft=ft,ontime=1-(1-ft)**k,ftb=ftb,k=k)
print("S1",json.dumps(out["S1"],default=str))
# per-hop rows
out["S1_hops"]=[(h["link"],h["te"],h["ta"],h["d"],h["p_loss"]) for h in rt["hops"]]+[(h["link"],h["te"],h["ta"],h["d"],h["p_loss"]) for h in back["hops"]]
# ---------------- E3: value move Earth -> Neptune at t=10 h
t0=10.0
v=mpx.best_route(HUB,"Neptune",t0,wait=True); r=mpx.best_route("Neptune",HUB,v["arrival"]+mpx.SEC,wait=True)
T0v=mpx.T0(v["route"],t0); Re=2*T0v+24
out["E3"]=dict(t0=t0,route=v["route"],arr=v["arrival"],held=sum(h["held_h"] for h in v["hops"]),ack_route=r["route"],ack_arr=r["arrival"],ft=first_try(v),abandon=mpx.route_abandon_bound([h["p_loss"] for h in v["hops"]])[0],Re=Re,
  hops=[(h["link"],h["te"],h["ta"],h["d"],h["p_loss"],h["held_h"]) for h in v["hops"]],ack_hops=[(h["link"],h["te"],h["ta"],h["d"],h["p_loss"]) for h in r["hops"]])
print("E3",json.dumps({k:v for k,v in out["E3"].items() if k not in("hops","ack_hops")},default=str))
# ---------------- S2: Neptune links during 0-30 h (maintenance B-Neptune 2-26) + forced loss at Relay A 0-6 h
ts_=np.arange(0,30,0.25); rows=[]
for t in ts_:
    a=bool(mpx.link_state("Neptune","Relay A",t)["ok"]); b=bool(mpx.link_state("Neptune","Relay B",t)["ok"])
    rows.append((t,a,b))
openA=[t for t,a,b in rows if a]; openB=[t for t,a,b in rows if b]
firstB=mpx.next_launch("Neptune","Relay B",0.0)
out["S2"]=dict(B_first_launch=firstB,A_open_all=all(a for _,a,_ in rows))
print("S2 Neptune->B first launch after maintenance:",fh(firstB),"A open every sample 0-30h:",out["S2"]["A_open_all"])
# 72 h isolation: missed batches, queue
iso=72.0; out["S2"]["missed_batches"]=iso/B
# ---------------- MPOR / book at B=3, G=1
sig=.30; z=2.33; ks=4
book=[("Bank A","Earth",0,600),("Bank B","Mars",None,400),("Operator C","Jupiter",None,250),("Operator D","Neptune",None,150)]
S=[];rows=[]
for n,w,rt_,N in book:
    RT=0 if w==HUB else D[w]["in_dmax"]+D[w]["out_dmax"]
    M=RT+G+B;f=sig*math.sqrt(M/8760);im=z*f*N;s=(ks-z)*f*N;S.append(s);rows.append([n,w,RT,M,N,im,s,f])
F=sum(sorted(S)[-2:]);sh=[F*s/sum(S) for s in S]
for r_,x in zip(rows,sh):r_.append(x)
out["book"]=rows;out["F"]=F;out["skin"]=.25*F
for r_ in rows:print(f"{r_[0]:11s} RT {r_[2]:.3f} MPOR {r_[3]:.3f} IM {r_[5]:.3f} ({100*r_[5]/r_[4]:.2f}%) S {r_[6]:.3f} share {r_[8]:.3f}")
print("F",F,"skin",.25*F,"IMtot",sum(r_[5] for r_ in rows),"escrow D",sh[3],"total",sum(r_[5] for r_ in rows)+F+sh[3]+.25*F+1)
# isolation: margin needed for 72h + MPOR vs posted
fD=rows[3][7]; need=z*sig*math.sqrt((72+rows[3][3])/8760)*150
out["S2"]["iso_need"]=need; out["S2"]["iso_posted"]=rows[3][5]
print("S2 72h isolation: 99% margin needed",need,"posted",rows[3][5])
# ---------------- E5 price runs
rng=np.random.default_rng(42); n=2_000_000
for name,M in [("Neptune",rows[3][3]),("Earth",rows[0][3])]:
    s=sig*math.sqrt(M/8760)
    g_=rng.standard_normal(n)*s
    nu=3; t_=rng.standard_t(nu,n)*s/math.sqrt(nu/(nu-2))
    for lab,x in (("gauss",g_),("t3",t_)):
        exc=np.mean(x<-z*s)  # long position loss beyond IM (one-sided)
        beyond_stress=np.mean(x<-ks*s)
        out.setdefault("E5",{})[f"{name}_{lab}"]=dict(M=M,exc=float(exc),beyond_stress=float(beyond_stress),q999=float(-np.quantile(x,.001)/s))
        print("E5",name,lab,"MPOR",round(M,2),"P(loss>IM)",round(exc,5),"P(loss>4sd)",round(beyond_stress,6),"99.9% loss in sd",round(-np.quantile(x,.001)/s,2))
json.dump(out,open("out/traces.json","w"),default=str,indent=1)
