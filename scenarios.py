"""Deterministic MPX protocol engine on mpx.py. No-loss scenario traces (S1), incident runs (S2),
access timelines (S3), ledgers and conservation (E3), later epochs (E5). Usage: run(epoch_h, ...)."""
import numpy as np, json, math, copy, mpx
HUB="Earth"; SEC=mpx.SEC; BF=1.5; FULL=6.0; G=1.0; Z=2.62; MM=0.8; LIQ=0.30; SIG=0.30
INNER={"Mercury","Venus","Earth","Mars","Ceres","Jupiter"}
D=json.load(open("out/design.json"))
def kcopies(x,t):
    if x==HUB: return 1
    rt_in=mpx.best_route(x,HUB,t,wait=True); rt_out=mpx.best_route(HUB,x,t,wait=True)
    f=min(np.prod([1-h["p_loss"] for h in rt_in["hops"]]),np.prod([1-h["p_loss"] for h in rt_out["hops"]]))
    return int(math.ceil(math.log(0.005)/math.log(1-f)))
def RT(x): return 0.0 if x==HUB else D[x]["in_dmax"]+D[x]["out_dmax"]
def IMf(x,N): return Z*SIG*math.sqrt((RT(x)+G+BF)/8760)*abs(N)
class Sim:
    def __init__(s,t0): s.t0=t0; s.msgs=[]; s.led={}; s.events=[]; s.enc=[]  # enc: (t, account, amount_encumbered)
    def send(s,kind,a,b,t,copies=1):
        if a==b: arr=t+SEC; s.msgs.append(dict(kind=kind,a=a,b=b,t=t,te=t,arr=arr,copies=copies,orig=0,launch=0,route="local",f=1.0,pk=1.0)); return arr
        rt=mpx.best_route(a,b,t,wait=True); n=len(rt["hops"]); f=float(np.prod([1-h["p_loss"] for h in rt["hops"]]))
        s.msgs.append(dict(kind=kind,a=a,b=b,t=t,te=rt["hops"][0]["te"],arr=rt["arrival"],copies=copies,orig=copies,launch=2*n*copies,route="-".join(r.replace("Relay ","R") for r in rt["route"]),held=sum(h["held_h"] for h in rt["hops"]),f=f,pk=1-(1-f)**copies))
        return rt["arrival"]
    def book(s,t,acct,delta,note): s.led[acct]=s.led.get(acct,0)+delta; s.events.append((t,acct,delta,note,s.led[acct]))
def summarize(s,t_end):
    # encumbrance series from events on accounts flagged as encumbered (hold:, margin:, transit:)
    ts=sorted(set([e[0] for e in s.events]+[t_end])); bal={}; series=[]
    ev=sorted(s.events,key=lambda e:e[0]); i=0
    for t in ts:
        while i<len(ev) and ev[i][0]<=t: bal[ev[i][1]]=bal.get(ev[i][1],0)+ev[i][2]; i+=1
        series.append((t,sum(v for k,v in bal.items() if k.split(":")[0] in("hold","margin","transit"))))
    peak=max(v for _,v in series); ah=sum(series[j][1]*(series[j+1][0]-series[j][0]) for j in range(len(series)-1))
    # asset-hours by asset use: held/in-flight credits vs margin collateral
    bycat={}
    for cat in ("hold","transit","margin"):
        bal={};i2=0;ser=[]
        for t in ts:
            while i2<len(ev) and ev[i2][0]<=t: bal[ev[i2][1]]=bal.get(ev[i2][1],0)+ev[i2][2]; i2+=1
            ser.append((t,sum(v for k,v in bal.items() if k.split(":")[0]==cat)))
        bycat[cat]=sum(ser[j][1]*(ser[j+1][0]-ser[j][0]) for j in range(len(ser)-1))
    pall=float(np.prod([m["pk"] for m in s.msgs]))
    return dict(ah_by=bycat,p_conditional=pall,peak_enc=peak,asset_hours=ah,orig=sum(m["orig"] for m in s.msgs),launch=sum(m["launch"] for m in s.msgs),nmsg=len(s.msgs))
def conserved(s):
    tot={}; 
    for t,acct,d,_,_ in s.events: tot[round(t,9)]=tot.get(round(t,9),0)+d
    return abs(sum(tot.values()))<1e-9 and all(abs(v)<1e-9 for v in tot.values())
# ---------------------------------------------------------------- scenarios
def value_move(t0,src="Neptune",dst="Mars",amt=100.0):
    s=Sim(t0); s.book(t0,f"cash:{src}",1000,"funding");s.book(t0,f"cash:{dst}",1000,"funding");s.book(t0,"funding:ext",-2000,"funding")
    s.book(t0,f"cash:{src}",-amt,"debit T1");s.book(t0,"flight:T1",amt,"T1 in transit")
    ta=s.send("TRANSFER T1",src,dst,t0)
    s.book(ta,"flight:T1",-amt,"T1 credited");s.book(ta,f"cash:{dst}",amt,"credit T1")
    tb=s.send("CREDITED T1",dst,src,ta)
    s.book(t0,f"hold:{src}",amt,"source hold (memo)");s.book(t0,f"memo:{src}",-amt,"hold contra");s.book(tb,f"hold:{src}",-amt,"final");s.book(tb,f"memo:{src}",amt,"final")
    r=summarize(s,tb); r.update(complete=tb-t0,credit=ta-t0,exposure=0.0,capital=amt,conserved=conserved(s)); return s,r
def long_lived(t0,x="Mercury",N=100.0,days=90):
    s=Sim(t0); im=IMf(x,N); te=t0+24*days
    s.book(t0,f"cash:{x}",1000,"funding");s.book(t0,"funding:ext",-1000,"funding")
    s.book(t0,f"cash:{x}",-im,"post IM");s.book(t0,"transit:IM",im,"IM in transit")
    ta=s.send("MARGIN DEPOSIT",x,HUB,t0,kcopies(x,t0)); s.book(ta,"transit:IM",-im,"");s.book(ta,f"margin:{x}",im,"IM at hub")
    s.send("CREDITED",HUB,x,ta)
    t_open=BF*math.ceil((ta-t0+0.2)/BF)+t0; s.send("ORDER open",x,HUB,t_open-1.0,kcopies(x,t_open-1)); s.send("FILL",HUB,x,t_open+0.9)
    # sessions: re-pin before predicted closures on the pinned route (handshake = 1 origination each way)
    repins=0; t=t_open; pinned=mpx.best_route(x,HUB,t)["route"]
    for d in range(days):
        tt=t0+24*(d+1); s.send("STATEMENT",HUB,x,tt)
        cur=mpx.best_route(x,HUB,tt+24,wait=False)
        if cur is None or tuple(cur["route"])!=tuple(pinned):
            s.send("SYN re-pin",x,HUB,tt); repins+=1; pinned=mpx.best_route(x,HUB,tt+24)["route"]
    tp=s.send("EXPIRY PAYOUT + IM release",HUB,x,te+0.9)
    s.book(te+0.9,f"margin:{x}",-im,"release");s.book(te+0.9,"transit:IM2",im,"");s.book(tp,"transit:IM2",-im,"");s.book(tp,f"cash:{x}",im,"IM returned")
    tc=s.send("CREDITED",x,HUB,tp)
    r=summarize(s,tc); r.update(complete=tc-t0,repins=repins,IM=im,exposure=0.0,capital=im,conserved=conserved(s)); return s,r
def price_run(t0,direction=+1,step=0.005,nsteps=8,N=150.0,long_at="Neptune",short_at="Mars",H=36.0,incident=None):
    s=Sim(t0); M={long_at:+N,short_at:-N}; im={x:IMf(x,N) for x in M}; liq={x:LIQ*N for x in M}
    for x in M:
        s.book(t0,f"cash:{x}",1000,"funding");s.book(t0,"funding:ext",-1000,"funding")
        s.book(t0,f"cash:{x}",-im[x],"IM posted earlier");s.book(t0,f"margin:{x}",im[x],"IM at hub")
    for x in M: s.send("ORDER open",x,HUB,t0-3.0,kcopies(x,t0-3.0)); s.send("FILL",HUB,x,t0+0.01)
    coll=dict(im); p=1.0; pend={x:None for x in M}; calls=0; shortfall={x:0.0 for x in M}; under_h={x:0.0 for x in M}; peak_exp=0.0
    t=t0; steps=int(H/BF)
    for k in range(1,steps+1):
        t=t0+k*BF; newp=p*(1+direction*step) if k<=nsteps else p
        for x in M:
            vm=M[x]*(newp-p); coll[x]+=vm; s.book(t,f"margin:{x}",vm,"VM"); s.book(t,"vm:pool",-vm,"VM")
            gap=max(0.0,im[x]-coll[x]); shortfall[x]=max(shortfall[x],gap); peak_exp=max(peak_exp,gap)
            if gap>0: under_h[x]+=BF
            if pend[x] is None and coll[x]<MM*im[x]:
                need=im[x]-coll[x]; c=kcopies(x,t); ta=s.send("MARGIN CALL",HUB,x,t,c); calls+=1
                s.book(ta,f"cash:{x}",-need,"top-up debit");s.book(ta,f"transit:TU{x}{k}",need,"")
                tb=s.send("TOP-UP",x,HUB,ta,c); s.book(tb,f"transit:TU{x}{k}",-need,"");s.book(tb,f"margin:{x}",need,"top-up")
                s.send("CREDITED",HUB,x,tb); pend[x]=(tb,need); liq[x]-=need
            if pend[x] is not None and t>=pend[x][0]: coll[x]+=pend[x][1]; pend[x]=None
        p=newp
    for x in M: s.book(t,f"margin:{x}",-coll[x] if False else 0,"")
    end=max(m["arr"] for m in s.msgs); r=summarize(s,end); r.update(liq_used={x:LIQ*N-liq[x] for x in M},calls=calls,peak_exposure=peak_exp,shortfall=shortfall,under_h=under_h,IM=im,liq_left=liq,capital=sum(im.values()),conserved=conserved(s),complete=end-t0)
    return s,r
def bind_step(direction,lo=0.0005,hi=0.05):
    """Smallest per-batch step at which a member's top-ups exhaust its liquidity (the binding constraint)."""
    for _ in range(18):
        m=(lo+hi)/2; _,r=price_run(0.0,direction,step=m)
        if min(r["liq_left"].values())<0: hi=m
        else: lo=m
    return hi
if __name__=="__main__":
    out={}
    for name,fn in (("value",lambda:value_move(0.0)),("long",lambda:long_lived(0.0)),("rise",lambda:price_run(0.0,+1)),("fall",lambda:price_run(0.0,-1))):
        s,r=fn(); out[name]=dict(r=r,msgs=s.msgs)
        print(name,{k:(round(v,3) if isinstance(v,float) else v) for k,v in r.items() if not isinstance(v,dict)})
    json.dump(out,open("out/s1.json","w"),indent=1,default=float)
