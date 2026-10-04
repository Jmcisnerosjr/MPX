"""
E4 long-horizon scan: 200 Julian years of all 38 directed backbone links.

Screening argument (rigorous under the brief's model):
  clr(t_e) = distance from the Sun to the segment [p_s(t_e), p_r(t_a(t_e))].
  Every point of that segment is (1-s) p_s + s p_r, so moving the endpoints by
  dp_s, dp_r moves each point by at most max(|dp_s|, |dp_r|). Distance-to-a-set is
  1-Lipschitz, hence
        |d clr / d t_e| <= V := max(v_s, v_r * dt_a/dt_e),
        dt_a/dt_e <= 1 + 2 v_max/(c - v_max)        (FIFO lemma)
  with v = perihelion speed (exact maximum on a Kepler ellipse).
  If a closure exists at t*, the nearest grid sample (|t_k - t*| <= D/2) has
        clr(t_k) <= 0.10 + V D / 2.
  So flagging every sample below 0.10 + V D/2 and refining around it cannot miss a closure.
  Refinement: 1-minute sub-grid on [t_k - D, t_k + D], then bisection of each boundary
  to 1e-7 h. The only feature that can hide inside a 1-minute sub-step is a flip of depth
  < V * (1 min)/2, i.e. < 1.2e-5 AU for the fastest link (Mercury).
"""
import numpy as np, json, os, time
import mpx

YEARS = 200
DT = 1.0                                   # hours, level-1 step
T_END = YEARS * 365.25 * 24
OUT = os.path.join(mpx.HERE, "out", "scan")
os.makedirs(OUT, exist_ok=True)
C_AU_H = 1.0 / mpx.LT_H


def V_bound(s, r):
    vs, vr = mpx.max_speed_au_per_h(s), mpx.max_speed_au_per_h(r)
    vmax = max(vs, vr)
    return max(vs, vr * (1 + 2 * vmax / (C_AU_H - vmax)))


def sun_clr(s, r, te):
    st = mpx.link_state(s, r, te, maintenance=False)
    return st["clr"], st["ta"] - te, st["d"]


def refine_boundary(s, r, t_open, t_closed):
    lo, hi = t_open, t_closed
    while abs(hi - lo) > 1e-7:
        mid = 0.5 * (lo + hi)
        if sun_clr(s, r, mid)[0] >= mpx.R_EXCL:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def scan_link(s, r):
    V = V_bound(s, r)
    thr = mpx.R_EXCL + V * DT / 2
    grid = np.arange(0.0, T_END + DT, DT)
    clr = np.empty_like(grid, dtype=np.float32)
    tau = np.empty_like(grid, dtype=np.float32)
    dist = np.empty_like(grid, dtype=np.float32)
    CH = 250_000
    for i in range(0, len(grid), CH):
        c, t, d = sun_clr(s, r, grid[i:i + CH])
        clr[i:i + CH], tau[i:i + CH], dist[i:i + CH] = c, t, d
    flagged = np.where(clr < thr)[0]
    # group flagged samples into windows, refine each on a 1-minute sub-grid
    closures, grazes = [], 0
    if len(flagged):
        groups = np.split(flagged, np.where(np.diff(flagged) > 1)[0] + 1)
        for g in groups:
            a = max(grid[g[0]] - DT, 0.0)
            b = min(grid[g[-1]] + DT, T_END)
            sub = np.arange(a, b + 1e-12, 1.0 / 60)
            c_sub = sun_clr(s, r, sub)[0]
            closed = c_sub < mpx.R_EXCL
            if not closed.any():
                grazes += 1
                continue
            edges = np.diff(closed.astype(int))
            starts = list(np.where(edges == 1)[0])
            ends = list(np.where(edges == -1)[0])
            if closed[0]:
                starts = [None] + starts
            if closed[-1]:
                ends = ends + [None]
            for st, en in zip(starts, ends):
                t0 = sub[0] if st is None else refine_boundary(s, r, sub[st], sub[st + 1])
                t1 = sub[-1] if en is None else refine_boundary(s, r, sub[en + 1], sub[en])
                closures.append((float(t0), float(t1), float(c_sub[(st or 0):(en or len(sub) - 1) + 1].min())))
    np.savez_compressed(os.path.join(OUT, f"{s}__{r}.npz".replace(" ", "_")),
                        clr=clr, tau=tau, d=dist)
    return dict(link=[s, r], V_au_per_h=V, threshold=thr, n_flag_windows=int(len(flagged) and len(groups)),
                grazes_no_closure=grazes, closures=closures,
                min_clr=float(clr.min()))


if __name__ == "__main__":
    t0 = time.time()
    res = []
    for (s, r) in mpx.DIRECTED:
        out = scan_link(s, r)
        n = len(out["closures"])
        durs = [c[1] - c[0] for c in out["closures"]]
        print(f"{s:>8s} -> {r:<8s} V={out['V_au_per_h']:.5f} closures={n:3d} "
              f"dur[h] min={min(durs) if durs else 0:7.2f} max={max(durs) if durs else 0:7.2f} "
              f"grazes={out['grazes_no_closure']} t={time.time()-t0:.0f}s", flush=True)
        res.append(out)
    json.dump(dict(years=YEARS, dt_h=DT, T_end_h=T_END, links=res),
              open(os.path.join(mpx.HERE, "out", "e4_closures.json"), "w"), indent=1)
