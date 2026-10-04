import json, numpy as np, scenarios as S
KN={"TRANSFER T1":("request from client; destination state unknown","debit 100 into hold; send TRANSFER"),
 "CREDITED T1":("first copy of T1; ID not tombstoned","credit 100; send CREDITED"),
 "MARGIN DEPOSIT":("margin requirement 1.93","send deposit as copies"),
 "CREDITED":("payment arrived","credit account; confirm"),
 "ORDER open":("own book; hub price one light-time old","send order as copies"),
 "FILL":("batch cleared","report fill"),
 "STATEMENT":("positions, margin","daily statement"),
 "SYN re-pin":("predicted closure on pinned link","re-pin session (1 origination)"),
 "EXPIRY PAYOUT + IM release":("expiry batch price","pay out; release margin"),
 "MARGIN CALL":("collateral below 80\\% of IM","call as copies"),
 "TOP-UP":("call received","debit cash; send top-up as copies")}
def actor(x): return "Earth hub" if x=="Earth" else ("Operator D" if x=="Neptune" else f"{x} ledger")
def bal(s,t,keys):
    b={}
    for tt,a,d,_,_ in s.events:
        if tt<=t+1e-9: b[a]=b.get(a,0)+d
    return "; ".join(f"{k.split(':')[1]} {k.split(':')[0]} {b.get(k,0):.2f}" for k in keys if abs(b.get(k,0))>1e-9 or k.startswith("cash"))
def esc(x): return x.replace("&","\\&").replace("%","\\%")
def rows(s,keys,merge_statements=False):
    R=[];stm=[m for m in s.msgs if m["kind"]=="STATEMENT"]
    for m in s.msgs:
        if merge_statements and m["kind"]=="STATEMENT":
            if m is not stm[0]: continue
            R.append((f"{m['t']:.2f}--{stm[-1]['t']:.2f}","Earth hub","positions, margin","daily statement $\\times$"+str(len(stm)),f"backbone, 1 each",f"{m['arr']-m['t']:.3f} h later","unchanged"));continue
        k,a=KN.get(m["kind"],("",m["kind"].lower()))
        pk=f" (p {m['pk']:.3f})" if m["route"]!="local" else ""
        F=lambda x: f"{x:.4f}" if x<100 else f"{x:.2f}"
        pkt="local" if m["route"]=="local" else "backbone "+m["route"].replace("-"," $\\to$ ")+", $\\times$"+str(m["copies"])+", launch "+format(m["te"],".4f")
        R.append((F(m['t']),actor(m["a"]),k,a,pkt,f"{F(m['arr'])}{pk}",bal(s,m["arr"],keys)))
    return R
def table(cap,R):
    head=r"""\setlength{\tabcolsep}{3pt}\begin{longtable}{@{}p{0.07\linewidth}p{0.09\linewidth}p{0.15\linewidth}p{0.14\linewidth}p{0.16\linewidth}p{0.1\linewidth}p{0.17\linewidth}@{}}
\caption{"""+cap+r"""}\\
\toprule
Time (h) & Actor & Knows & Action & Packet & Arrival & State after \\
\midrule
\endhead
"""
    return head+"\n".join(" & ".join(esc(c) if i not in (0,4,5) else c for i,c in enumerate(r))+r" \\" for r in R)+"\n\\bottomrule\n\\end{longtable}\n"
out=[]
s,r=S.value_move(0.0); out.append(table("S1 value move, conditional no-loss trace: Neptune client sends 100 credits to a Mars client. Probabilities are first-try success per message.",rows(s,["cash:Neptune","hold:Neptune","flight:T1","cash:Mars"])))
s,r=S.long_lived(0.0); out.append(table("S1 long-lived obligation, conditional trace: Mercury client holds a 90-day futures position, notional 100.",rows(s,["cash:Mercury","margin:Mercury","transit:IM"],True)))
for d,lab in ((+1,"rising"),(-1,"falling")):
    s,r=S.price_run(0.0,d); out.append(table(f"S1 {lab} price run, conditional trace: Neptune long 150 against Mars short 150; price {'+' if d>0 else '-'}0.5\\% per fast batch for 12 h.",rows(s,["cash:Neptune","margin:Neptune","cash:Mars","margin:Mars"])))
for x in ("Ceres","Jupiter","Neptune"):
    s,r=S.value_move(0.0,src=x,dst="Earth"); out.append(table(f"S3 value move from {x} to the Earth hub, conditional trace from reset.",rows(s,[f"cash:{x}",f"hold:{x}","cash:Earth"])))
open("out/traces.tex","w").write("\n".join(out)); print("ok",sum(t.count("\\\\") for t in out))
