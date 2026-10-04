"""
S3 access tables, E2 probability machinery, E5 shifted epochs. Writes out/s3_e2_e5.md.
HUB is a design choice (see hub ranking in e4_tables.md / hub_loss.md).
"""
import numpy as np, os, json, itertools
import mpx

HUB = "Earth"
S = [s for s in mpx.SETTLEMENTS]
OUT = os.path.join(mpx.HERE, "out")
L = []


def hop_retry_timer(flight_h):
    return 2 * flight_h + 1.0          # R_h = 2 x one-way flight + 60 min


def route_delivery_cdf(rt, T_list):
    """Exact P(packet delivered end-to-end by elapsed time T) for ONE endpoint attempt,
    no queueing, geometry frozen at the first-launch values (flights change by << R_h between
    retries). Hop k's data is delivered on launch j (j=1..4) w.p. p^(j-1)(1-p); launch j leaves
    (j-1)*R_h after launch 1 (a lost receipt also triggers a retry, but the data copy that already
    arrived is the one forwarded, so delivery time depends only on data loss)."""
    hops = rt["hops"]
    opts = []
    for h in hops:
        p, fl = h["p_loss"], h["ta"] - h["te"]
        Rh = hop_retry_timer(fl)
        opts.append([((j) * Rh, (p ** j) * (1 - p)) for j in range(4)])
    base = rt["delay_h"]
    out = []
    for T in T_list:
        tot = 0.0
        for combo in itertools.product(*opts):
            extra = sum(c[0] for c in combo)
            pr = np.prod([c[1] for c in combo])
            if base + extra <= T:
                tot += pr
        out.append(tot)
    return out


def endpoint_cdf(rt, T, t0_h):
    """Up to 4 endpoint attempts spaced R_e = 2*T0 + 24 h; launches independent (brief Sec 4).
    P(delivered by T) = 1 - prod_k (1 - F(T - k*R_e))."""
    Re = 2 * mpx.T0(rt["route"], t0_h) + 24.0
    miss = 1.0
    for k in range(4):
        if T - k * Re >= 0:
            miss *= 1 - route_delivery_cdf(rt, [T - k * Re])[0]
    return 1 - miss, Re


# ------------------------------------------------------------------ S3
def s3_table(h, maintenance=True):
    rows = []
    for x in S:
        if x == HUB:
            rows.append((x, "local", "local access 1 s", 0.0, 0.0, 1.0, 1.0, "-", "-"))
            continue
        for a, b, tag in ((x, HUB, "to hub"), (HUB, x, "from hub")):
            rt = mpx.best_route(a, b, h, wait=True, maintenance=maintenance)
            av = mpx.route_availability(rt["route"], h, maintenance=maintenance)
            ps = [hp["p_loss"] for hp in rt["hops"]]
            ab_exact, _ = mpx.route_abandon_bound(ps)
            rows.append((x, tag, " -> ".join(r.replace("Relay ", "R") for r in rt["route"]),
                         rt["delay_h"], sum(hp["held_h"] for hp in rt["hops"]), av,
                         mpx.route_first_try(ps), " / ".join(f"{p:.3f}" for p in ps), f"{ab_exact:.2e}"))
    return rows


def s3_md(h, maintenance=True):
    out = [f"### Hour {h}{'' if maintenance else ' (maintenance removed)'}", "",
           "| Settlement | Dir | Best route | One-way (h) | Held (h) | Avail 24 h | First-try success | Loss per launch | P(abandon) |",
           "|---|---|---|---|---|---|---|---|---|"]
    for r in s3_table(h, maintenance):
        out.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]:.4f} | {r[4]:.3f} | {r[5]:.3f} | "
                   f"{r[6]:.4f} | {r[7]} | {r[8]} |")
    return out + [""]


L += [f"# S3, E2, E5 (clearing hub = {HUB})", "",
      "Delays are gateway-to-gateway backbone, no loss, empty queues, 1 s serialization per launch, "
      "1 s relay processing. Add 1 s local access at each end for a client talking to its local "
      "operator. Availability = fraction of the 24 h from h in which every link of that route can "
      "launch (1-minute grid, maintenance applied). P(abandon) = 1 - prod(1 - p_i^4), exact.", "",
      "## S3 access tables", ""]
L += s3_md(0.0) + s3_md(300.0)

# S3 best / median / worst by delay into hub at h0
ranked = sorted([(mpx.best_route(x, HUB, 0.0)["delay_h"], x) for x in S if x != HUB])
picks = [ranked[0], ranked[len(ranked) // 2], ranked[-1]]
L += ["### Best / median / worst settlement by one-way delay into the hub at h 0", ""]
for d, x in picks:
    rt = mpx.best_route(x, HUB, 0.0)
    L += [f"**{x}** ({d:.4f} h)", "", "| Hop | ready | t_e | t_a | d (AU) | p loss | R_h (h) |",
          "|---|---|---|---|---|---|---|"]
    for hp in rt["hops"]:
        L.append(f"| {hp['link'][0]} -> {hp['link'][1]} | {mpx.fmt_h(hp['ready'])} | {mpx.fmt_h(hp['te'])} | "
                 f"{mpx.fmt_h(hp['ta'])} | {hp['d']:.4f} | {hp['p_loss']:.4f} | "
                 f"{hop_retry_timer(hp['ta']-hp['te']):.3f} |")
    L.append("")

# ------------------------------------------------------------------ E2
L += ["## E2 probability machinery", "",
      "### Backbone: delivery-by-deadline into the hub at h 0 (exact, single endpoint attempt, then with endpoint retries)", "",
      "| From | Conditional (no-loss) | P by +0 | P by +R_h1 | P by 24 h | P by 48 h | R_e (h) | P by R_e + cond | P never (4 attempts) |",
      "|---|---|---|---|---|---|---|---|---|"]
for x in S:
    if x == HUB:
        continue
    rt = mpx.best_route(x, HUB, 0.0)
    base = rt["delay_h"]
    Rh1 = hop_retry_timer(rt["hops"][0]["ta"] - rt["hops"][0]["te"])
    c = route_delivery_cdf(rt, [base, base + Rh1, 24.0, 48.0, 1e9])
    Re = 2 * mpx.T0(rt["route"], 0.0) + 24.0
    pRe, _ = endpoint_cdf(rt, Re + base, 0.0)
    q = 1 - c[-1]
    L.append(f"| {x} | {base:.3f} h | {c[0]:.4f} | {c[1]:.4f} | {c[2]:.5f} | {c[3]:.6f} | {Re:.2f} | "
             f"{pRe:.7f} | {q**4:.2e} |")
L += ["", "P never (4 attempts) = (single-attempt failure)^4, valid because every launch is an "
      "independent trial (brief Sec 4) and attempts are separated by R_e > 24 h. It ignores closures "
      "and incidents, which are not independent and are handled by the rulebook, not by probability.", ""]

L += ["### Direct service to the hub at h 0: P(at least one of n copies arrives), copies >= 60 s apart", "",
      "| From | One-way (min) | Sun clear (AU) | p per copy | n=1 | n=2 | n=4 | n=12 (whole daily quota) |",
      "|---|---|---|---|---|---|---|---|"]
for x in S:
    if x == HUB:
        continue
    ps = []
    for k in range(12):
        st = mpx.link_state(x, HUB, 2 * mpx.SEC + k * 60 * mpx.SEC, maintenance=False)
        ps.append(float(mpx.loss_direct(st["d"])))
    st0 = mpx.link_state(x, HUB, 2 * mpx.SEC, maintenance=False)
    def succ(n): return 1 - np.prod(ps[:n])
    L.append(f"| {x} | {(float(st0['ta'])-2*mpx.SEC)*60:.2f} | {float(st0['clr']):.3f} | {ps[0]:.4f} | "
             f"{succ(1):.4f} | {succ(2):.4f} | {succ(4):.4f} | {succ(12):.4f} |")
L += ["", "Each copy's loss is recomputed from its own launch geometry (60 s later). The direct service "
      "has no receipt: the sender learns nothing unless a reply (another lossy flight) comes back.", ""]

# capacity / quota bounds
L += ["### Link capacity and queue checks", "",
      "- Max originations: 600 per rolling 24 h system-wide. Each origination makes at most 4 launches per "
      "hop plus one receipt per arriving copy, i.e. <= 8 launches per link it crosses. Worst case on one "
      "directed link: 600 x 8 = 4,800 launch-seconds per 86,400 s = 5.6% utilisation. Capacity never binds "
      "in steady state.",
      "- Queue cap 10,000: a pinned route held behind a closed link accumulates at most the quota rate, "
      "600/day, so the cap is reached only after 16.7 days of holding. Packets live 30 days. Closures last "
      "up to 4,566 h (190 days, Ceres-Relay A), so a session pinned through a closing link WILL overflow "
      "and drop the newest packets. Rule: re-pin sessions before every predicted closure; the "
      "Two-Relay Lemma guarantees the other relay is open.", ""]

# ------------------------------------------------------------------ E5
L += ["## E5 shifted epochs (maintenance not replayed)", ""]
for yrs in (1, 10, 100):
    h = yrs * 365.25 * 24
    L += [f"### +{yrs} Julian year{'s' if yrs > 1 else ''} (t = {h:.0f} h)", ""]
    L += s3_md(h, maintenance=False)[2:]
summ = json.load(open(os.path.join(OUT, "e4_summary.json")))
L += ["### Difficult epoch candidates (from the 200-year scan)", ""]
T_END = 200 * 365.25 * 24
grid = np.arange(0.0, T_END + 1, 1.0)
# worst delay Neptune -> hub: recompute on hourly grid using scan arrays
SC = os.path.join(OUT, "scan")
def ld(s, r):
    z = np.load(os.path.join(SC, f"{s}__{r}.npz".replace(" ", "_")))
    return z["clr"].astype(float), z["tau"].astype(float)
Lk = {k: ld(*k) for k in mpx.DIRECTED}
def route_arr(route):
    ready = grid.copy(); ok = np.ones(len(grid), bool)
    for k in range(len(route) - 1):
        c, t = Lk[(route[k], route[k + 1])]; te = ready + mpx.SEC
        ok &= np.interp(te, grid, c) >= 0.1; ta = te + np.interp(te, grid, t)
        ready = ta + mpx.SEC if route[k + 1] in mpx.RELAYS else ta
    return ready - grid, ok
best = np.full(len(grid), np.inf); n3 = np.zeros(len(grid), bool)
for rt in mpx.candidate_routes("Neptune", HUB):
    d, ok = route_arr(rt); m = ok & (d < best); best = np.where(m, d, best)
    n3 = np.where(m, len(rt) == 4, n3)
i = int(np.argmax(best))
L += [f"- Worst Neptune -> {HUB} best-route delay in 200 yr: {best[i]:.3f} h at t = {grid[i]:.0f} h "
      f"({grid[i]/8766:.2f} yr). Fraction of time Neptune -> {HUB} is forced onto a 3-link route: {n3.mean():.3f}.", ""]
cl = json.load(open(os.path.join(OUT, "e4_closures.json")))
longest = max(((tuple(x["link"]), c) for x in cl["links"] for c in x["closures"]), key=lambda z: z[1][1] - z[1][0])
L += [f"- Longest single-link closure: {longest[0][0]} -> {longest[0][1]}, {longest[1][0]:.0f} h to "
      f"{longest[1][1]:.0f} h ({(longest[1][1]-longest[1][0])/24:.1f} days, starting {longest[1][0]/8766:.2f} yr).", ""]
first_hub = sorted((c[0], tuple(x["link"])) for x in cl["links"] if HUB in x["link"] for c in x["closures"])[:4]
L += [f"- First {HUB} relay-link closures: " + "; ".join(f"{l[0]}->{l[1]} at {t:.1f} h" for t, l in first_hub), ""]

open(os.path.join(OUT, "s3_e2_e5.md"), "w").write("\n".join(L))
print("\n".join(L))
