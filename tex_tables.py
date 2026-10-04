import numpy as np, json, math, mpx
fh=lambda h: mpx.fmt_h(h).replace("h","\\,h\\,").replace("m","\\,m\\,").replace("s","\\,s")
def F(h):
    H=int(h);m=(h-H)*60;M=int(m);s=(m-M)*60;return f"{H:02d}:{M:02d}:{s:04.1f}"
HUB="Earth";D=json.load(open("out/design.json"));L=[]
def fl(a,b,t):return float(mpx.light_time(a,b,t)[0])-t
# ---------- S1
Tb=6.0;floor=fl("Neptune",HUB,Tb);dl=Tb+floor;k=D["Neptune"]["k"];safety=60/3600
rt0=mpx.best_route("Neptune",HUB,Tb-5,wait=False);ts=dl-rt0["delay_h"]-safety-(k-1)*mpx.SEC
rt=mpx.route_timing(rt0["route"],ts,wait=False);last=rt["arrival"]+(k-1)*mpx.SEC
fls={x:fl(x,HUB,Tb) for x in D};clear=Tb+max(fls.values())
bk=mpx.route_timing(["Earth","Relay A","Neptune"],clear,wait=True)
ft=np.prod([1-h["p_loss"] for h in rt["hops"]])
rows=[("Client order reaches Operator D (local access)",F(ts-1/3600)),("Cutoff $T_b$",F(Tb)),
 (f"Operator D launches copy 1 of {k} toward Relay A",F(rt["hops"][0]["te"])),
 ("Copy 1 reaches Relay A",F(rt["hops"][0]["ta"])),("Copy 1 leaves Relay A",F(rt["hops"][1]["te"])),
 (f"Copy {k} reaches the hub",F(last)),(f"Deadline $T_b+\\ell^{{\\mathrm{{floor}}}}$ ({floor:.4f}\\,h)",F(dl)),
 ("Batch clears (farthest floor: Neptune)",F(clear)),("Fill report leaves the hub",F(bk["hops"][0]["te"])),
 ("Fill report reaches Neptune",F(bk["arrival"]))]
L.append("%S1\n"+"\n".join(f"{a} & {b} \\\\" for a,b in rows))
L.append(f"%S1 facts: route {rt['route']} delay {rt0['delay_h']:.4f} excess {(rt0['delay_h']-floor)*60:.1f} min; first-try per copy {ft:.4f}; on-time {1-(1-ft)**k:.4f}; send-to-fill {bk['arrival']-ts:.3f} h; margin before deadline {(dl-last)*3600:.1f} s; floors "+", ".join(f"{x} {v:.3f}" for x,v in fls.items()))
# ---------- E3
t0=10.0;v=mpx.best_route(HUB,"Neptune",t0);a=mpx.best_route("Neptune",HUB,v["arrival"]+mpx.SEC)
ab=mpx.route_abandon_bound([h["p_loss"] for h in v["hops"]])[0];Re=2*mpx.T0(v["route"],t0)+24
rows=[("Earth ledger debits 10 credits and locks them","10:00:00.0","on Earth (locked)"),
 ("Transfer launched Earth $\\to$ Relay A",F(v["hops"][0]["te"]),"in flight"),
 ("Relay A forwards",F(v["hops"][1]["te"]),"in flight"),
 ("Neptune ledger credits",F(v["arrival"]),"on Neptune"),
 ("Confirmation launched Neptune $\\to$ Relay A",F(a["hops"][0]["te"]),"on Neptune"),
 ("Earth receives confirmation, lock released",F(a["arrival"]),"on Neptune (final)")]
L.append("%E3\n"+"\n".join(f"{x} & {y} & {z} \\\\" for x,y,z in rows))
L.append(f"%E3 facts: first-try {np.prod([1-h['p_loss'] for h in v['hops']]):.4f} abandon {ab:.4f} Re {Re:.2f}h lossA {v['hops'][0]['p_loss']:.4f} lossN {v['hops'][1]['p_loss']:.4f}")
open("out/tex_tables.txt","w").write("\n".join(L));print("\n".join(L))
