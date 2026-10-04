"""
E4 analysis on top of e4_scan.py output.
  1. Two-Relay Lemma check: no gateway ever has both relay links closed by the Sun.
  2. Relay A <-> Relay B never closes (analytic: clearance = 2 AU).
  3. All-pairs best-route one-way delay over 200 years (hourly), open-at-hop-time rule.
  4. Hub choice: minimise worst-case delay to/from every settlement.
  5. Poorest service: worst 24 h window of best fixed-route availability per settlement.
Writes out/e4_tables.md and out/e4_summary.json.
"""
import numpy as np, json, os
import mpx

OUT = os.path.join(mpx.HERE, "out")
SCAN = os.path.join(OUT, "scan")
cl = json.load(open(os.path.join(OUT, "e4_closures.json")))
T_END, DT = cl["T_end_h"], cl["dt_h"]
grid = np.arange(0.0, T_END + DT, DT)
N = len(grid)
S = mpx.SETTLEMENTS


def load(s, r):
    z = np.load(os.path.join(SCAN, f"{s}__{r}.npz".replace(" ", "_")))
    return z["clr"].astype(np.float64), z["tau"].astype(np.float64)


L = {(s, r): load(s, r) for (s, r) in mpx.DIRECTED}
closures = {tuple(x["link"]): x["closures"] for x in cl["links"]}


def intervals_intersect(a, b):
    out, i, j = [], 0, 0
    a, b = sorted(a), sorted(b)
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


lines = ["# E4 long-horizon scan results", "",
         f"Span: {cl['years']} Julian years ({T_END:.0f} h), level-1 step {DT} h, "
         "refinement 1 min, boundaries bisected to 1e-7 h.", ""]

# ---- 1 & 2: link closure table + two-relay lemma
lines += ["## Link closures (Sun exclusion only)", "",
          "| Link | V bound (AU/h) | Screen thr (AU) | Closures | Shortest (h) | Longest (h) | Open fraction | Min clearance (AU) |",
          "|---|---|---|---|---|---|---|---|"]
for x in cl["links"]:
    s, r = x["link"]
    c = x["closures"]
    d = [b - a for a, b, _ in c]
    frac = 1 - sum(d) / T_END
    lines.append(f"| {s} -> {r} | {x['V_au_per_h']:.5f} | {x['threshold']:.5f} | {len(c)} | "
                 f"{min(d) if d else float('nan'):.2f} | {max(d) if d else float('nan'):.2f} | "
                 f"{frac:.5f} | {x['min_clr']:.4f} |")
shortest = min(b - a for x in cl["links"] for a, b, _ in x["closures"])
lines += ["", f"Shortest closure on any link in 200 years: **{shortest:.2f} h**, "
          f"{shortest/DT:.0f}x the scan step.", ""]

lemma = {}
for g in S:
    both_out = intervals_intersect([c[:2] for c in closures[(g, "Relay A")]],
                                   [c[:2] for c in closures[(g, "Relay B")]])
    both_in = intervals_intersect([c[:2] for c in closures[("Relay A", g)]],
                                  [c[:2] for c in closures[("Relay B", g)]])
    lemma[g] = (len(both_out), len(both_in))
lines += ["## Two-Relay Lemma check", "",
          "Hours in which a gateway has BOTH relay links closed (outbound / inbound):", ""]
lines += [f"- {g}: {lemma[g][0]} / {lemma[g][1]} overlapping intervals" for g in S]
lines += [""]

# ---- 3: all-pairs best-route delay (open at each hop's own emission time, no waiting)
SEC = mpx.SEC


def interp(arr, t):
    return np.interp(t, grid, arr)


def route_arrays(route):
    t = grid.copy()
    ready = t.copy()
    ok = np.ones(N, dtype=bool)
    for k in range(len(route) - 1):
        clr, tau = L[(route[k], route[k + 1])]
        te = ready + SEC
        ok &= interp(clr, te) >= mpx.R_EXCL
        ta = te + interp(tau, te)
        ready = ta + SEC if route[k + 1] in mpx.RELAYS else ta
    return ready - t, ok


pair = {}
for x in S:
    for y in S:
        if x == y:
            continue
        best = np.full(N, np.inf)
        for rt in mpx.candidate_routes(x, y):
            d, ok = route_arrays(rt)
            best = np.where(ok & (d < best), d, best)
        pair[(x, y)] = best
        assert np.isfinite(best).all(), f"no open route {x}->{y} at some hour"

lines += ["## Best-route one-way delay ranges over 200 years (hours, gateway to gateway)", "",
          "Every ordered pair had an open route at every sampled hour.", "",
          "| From \\ To | " + " | ".join(S) + " |", "|---" * (len(S) + 1) + "|"]
for x in S:
    row = []
    for y in S:
        if x == y:
            row.append("-")
        else:
            b = pair[(x, y)]
            row.append(f"{b.min():.2f}-{b.max():.2f}")
    lines.append(f"| {x} | " + " | ".join(row) + " |")
lines.append("")

# ---- 4: hub choice
hub_rows = []
for h in S:
    worst_in = max(pair[(x, h)].max() for x in S if x != h)
    worst_rt = max(pair[(x, h)].max() + pair[(h, x)].max() for x in S if x != h)
    mean_rt = np.mean([pair[(x, h)].mean() + pair[(h, x)].mean() for x in S if x != h])
    hub_rows.append((h, worst_in, worst_rt, mean_rt))
hub_rows.sort(key=lambda r: r[2])
lines += ["## Clearing-hub choice", "",
          "| Hub | Worst one-way into hub (h) | Worst round trip, any settlement (h) | Mean round trip (h) |",
          "|---|---|---|---|"]
lines += [f"| {h} | {a:.2f} | {b:.2f} | {c:.2f} |" for h, a, b, c in hub_rows]
lines.append("")

# ---- 5: poorest service: worst 24 h window of the best single pinned route, to each hub candidate
HUB = hub_rows[0][0]
W = int(24 / DT)
poor = []
for x in S:
    if x == HUB:
        continue
    best_window = np.zeros(N - W)
    for rt in mpx.candidate_routes(x, HUB):
        _, ok = route_arrays(rt)
        cs = np.concatenate([[0], np.cumsum(ok)])
        avail = (cs[W:N] - cs[:N - W]) / W
        best_window = np.maximum(best_window, avail)
    i = int(np.argmin(best_window))
    poor.append((x, float(best_window.min()), float(grid[i]), float(best_window.mean())))
lines += [f"## Poorest service to hub {HUB}", "",
          "Availability = fraction of a 24 h window in which every link of ONE fixed route is open "
          "(best route chosen per window). Hourly sampling.", "",
          "| Settlement | Worst 24 h availability | Window start (h) | Mean |", "|---|---|---|---|"]
lines += [f"| {x} | {a:.3f} | {t:.0f} ({t/8766:.2f} yr) | {m:.4f} |" for x, a, t, m in poor]
lines.append("")

open(os.path.join(OUT, "e4_tables.md"), "w").write("\n".join(lines))
json.dump(dict(shortest_closure_h=shortest, lemma=lemma, hub_ranking=hub_rows, hub=HUB, poorest=poor,
               pair_ranges={f"{x}->{y}": [float(v.min()), float(v.max()), float(v.mean())]
                            for (x, y), v in pair.items()}),
          open(os.path.join(OUT, "e4_summary.json"), "w"), indent=1)
print("\n".join(lines))
