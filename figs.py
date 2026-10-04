import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, scenarios as S
plt.rcParams.update({"font.family":"serif","font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7})
A=np.load("out/tiers_arrays.npz"); g=A["g"]; SET=["Mercury","Venus","Mars","Ceres","Jupiter","Saturn","Uranus","Neptune"]
yr=(g/8766).astype(int); Y=200
H=np.array([[A["DI_"+x][yr==k].max() for k in range(Y)] for x in SET])
fig,ax=plt.subplots(figsize=(6.4,1.9)); im=ax.imshow(H,aspect="auto",cmap="viridis",interpolation="nearest",extent=[0,Y,len(SET)-0.5,-0.5])
ax.set_yticks(range(len(SET)));ax.set_yticklabels(SET);ax.set_xlabel("Julian years after the epoch")
cb=fig.colorbar(im,ax=ax,pad=0.01);cb.set_label("worst one-way to Earth (h)",fontsize=7)
ax.set_title("Best-route delay to the Earth hub, worst hour of each year: always finite, Neptune always hardest")
fig.tight_layout();fig.savefig("out/fig_heatmap.pdf");plt.close()
P=json.load(open("out/phase.json"));LQ=P["LQ"];TR=P["TR"]
fig,axs=plt.subplots(1,2,figsize=(6.4,2.2),sharey=True)
for ax,d in zip(axs,("rising","falling")):
    M=np.array(P[d]); im=ax.imshow(M,origin="lower",cmap="RdYlGn_r",vmin=0,vmax=1,aspect="auto")
    ax.set_xticks(range(len(LQ)));ax.set_xticklabels([f"{int(l*100)}%" for l in LQ]);ax.set_yticks(range(len(TR)));ax.set_yticklabels([f"{int(t*100)}%" for t in TR])
    for i in range(len(TR)):
        for j in range(len(LQ)): ax.text(j,i,f"{M[i,j]*100:.0f}" if M[i,j]>=0.01 else ("<1" if M[i,j]>0 else "0"),ha="center",va="center",fontsize=6.5,color="black")
    ax.add_patch(plt.Rectangle((3-0.5,-0.5),1,len(TR),fill=False,lw=1.6,ec="black"))
    ax.set_title(f"{d.capitalize()} prices: P(default) %");ax.set_xlabel("member liquidity (share of notional)")
axs[0].set_ylabel("trend over 72 h\n(plus 10% jump, $t_3$, isolation)")
fig.tight_layout();fig.savefig("out/fig_phase.pdf");plt.close()
fig,axs=plt.subplots(1,2,figsize=(6.4,2.0),sharey=True)
for ax,(d,lab) in zip(axs,((+1,"Rising run"),(-1,"Falling run"))):
    s,r=S.price_run(0.0,d); ev=sorted(s.events,key=lambda e:e[0])
    for acct,col,name in (("margin:Neptune","#4F6FD0","Neptune (long)"),("margin:Mars","#C9573A","Mars (short)")):
        ts=[0.0];bs=[0.0];b=0.0
        for t,a,dl,_,_ in ev:
            if a==acct: b+=dl; ts.append(t); bs.append(b)
        ax.step(ts,bs,where="post",color=col,lw=1.3,label=name); ax.axhline(r["IM"][acct.split(":")[1]],color=col,lw=0.6,ls=":")
    ax.axhline(0,color="black",lw=0.5); ax.set_title(lab); ax.set_xlabel("hours after the epoch"); ax.set_xlim(0,30); ax.set_xticks([0,6,12,18,24,30]); ax.ticklabel_format(axis="x",style="plain",useOffset=False)
axs[0].set_ylabel("collateral at hub (credits)"); axs[1].legend(fontsize=6.5,loc="upper right",frameon=False)
fig.tight_layout();fig.savefig("out/fig_runs.pdf");plt.close();print("figures done")
