"""
mpx.py  Core model for the MultiPlanetary Exchange System brief (v1.3 data).

Every fixed rule here is quoted from the Participant Brief, Section 2-5:
  frame        Sun-centered J2000 ecliptic, AU; independent fixed Kepler ellipses
  light time   t_a - t_e = 8.317 min/AU * |r_recv(t_a) - r_send(t_e)|, tol 1 ms
  exclusion    straight segment sender@t_e -> receiver@t_a must miss r < 0.10 AU
  launch       1 s serialization; t_e = ready + 1 s; relay adds 1 s before next ready
  maintenance  launch fails if [t_e, t_a] overlaps the window (both directions)
  loss         backbone 1-exp(-0.02 d), direct 1-exp(-0.08 d), d = photon path AU

Public time unit is HOURS since epoch (2026-09-22 00:00 TDB). Distances in AU.
"""
import json, math, itertools, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

LIGHT_MIN_PER_AU = 8.317
LT_H = LIGHT_MIN_PER_AU / 60.0          # hours per AU
SEC = 1.0 / 3600.0                      # one second in hours
R_EXCL = 0.10                           # AU
LT_TOL_H = 1e-9                         # 3.6 microseconds, far inside the 1 ms rule

SETTLEMENTS = ["Mercury", "Venus", "Earth", "Mars", "Ceres",
               "Jupiter", "Saturn", "Uranus", "Neptune"]
RELAYS = ["Relay A", "Relay B"]
NODE_ID = {n: i + 1 for i, n in enumerate(SETTLEMENTS + RELAYS)}

_els = {b["name"]: b for b in json.load(open(os.path.join(DATA, "orbital_elements.json")))["bodies"]}
_net = json.load(open(os.path.join(DATA, "network_model.json")))
for r in _net["relays"]:
    _els[r["name"]] = r
assert abs(_net["light_minutes_per_au"] - LIGHT_MIN_PER_AU) < 1e-12
assert abs(_net["solar_exclusion_radius_au"] - R_EXCL) < 1e-12

MAINT = {}                                # frozenset{a,b} -> list of (start_h, end_h)
for m in _net["maintenance"]:
    MAINT.setdefault(frozenset(m["edge"]), []).append((m["start_hours"], m["end_hours"]))

# 19 two-way candidate links -> 38 directed links
LINKS = [(g, r) for g in SETTLEMENTS for r in RELAYS] + [("Relay A", "Relay B")]
DIRECTED = [(a, b) for a, b in LINKS] + [(b, a) for a, b in LINKS]


# ---------------------------------------------------------------- propagation
def _rot(b):
    d = math.radians
    w, i, O = d(b["arg_peri_deg"]), d(b["i_deg"]), d(b["node_deg"])
    cw, sw, ci, si = math.cos(w), math.sin(w), math.cos(i), math.sin(i)
    cO, sO = math.cos(O), math.sin(O)
    # columns: P (perihelion direction), Q (90 deg ahead in orbit plane)
    P = np.array([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si])
    Q = np.array([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si])
    return P, Q

_ROT = {n: _rot(b) for n, b in _els.items()}


def kepler_E(M, e):
    """Solve E - e sin E = M (radians, arrays ok). Newton, converged to 1e-15."""
    M = np.asarray(M, dtype=float)
    E = M + e * np.sin(M) if e < 0.8 else np.full_like(M, math.pi)
    for _ in range(30):
        dE = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E = E - dE
        if np.all(np.abs(dE) < 1e-15):
            break
    return E


def pos(name, t_h):
    """Heliocentric J2000-ecliptic position (AU) at t_h hours since epoch. Shape (..., 3)."""
    b = _els[name]
    t_d = np.asarray(t_h, dtype=float) / 24.0
    M = np.radians(np.mod(b["mean_anomaly_deg"] + b["mean_motion_deg_day"] * t_d, 360.0))
    a, e = b["a_au"], b["e"]
    E = kepler_E(M, e)
    xp = a * (np.cos(E) - e)
    yp = a * math.sqrt(1 - e * e) * np.sin(E)
    P, Q = _ROT[name]
    return xp[..., None] * P + yp[..., None] * Q


def max_speed_au_per_h(name):
    """Perihelion speed n*a*sqrt((1+e)/(1-e)), an exact upper bound on |v| for the ellipse."""
    b = _els[name]
    n = math.radians(b["mean_motion_deg_day"]) / 24.0
    return n * b["a_au"] * math.sqrt((1 + b["e"]) / (1 - b["e"]))


def elements(name):
    return dict(_els[name])


# ---------------------------------------------------------------- light time
def light_time(sender, receiver, te_h):
    """Solve t_a for emission t_e (hours). Returns (t_a, d_AU, p_send, p_recv). Arrays ok.
    Fixed-point iteration; contraction constant <= v_recv/c ~ 1e-4, so it converges fast."""
    te = np.asarray(te_h, dtype=float)
    ps = pos(sender, te)
    ta = te + LT_H * np.linalg.norm(pos(receiver, te) - ps, axis=-1)
    for _ in range(20):
        pr = pos(receiver, ta)
        new = te + LT_H * np.linalg.norm(pr - ps, axis=-1)
        done = np.max(np.abs(new - ta)) < LT_TOL_H
        ta = new
        if done:
            break
    pr = pos(receiver, ta)
    d = np.linalg.norm(pr - ps, axis=-1)
    return ta, d, ps, pr


def seg_clearance(A, B):
    """Closest distance from the origin (Sun) to segment A->B. Arrays (...,3)."""
    D = B - A
    dd = np.einsum("...i,...i->...", D, D)
    s = np.clip(-np.einsum("...i,...i->...", A, D) / dd, 0.0, 1.0)
    return np.linalg.norm(A + s[..., None] * D, axis=-1)


def link_state(sender, receiver, te_h, maintenance=True):
    """For emission time(s) t_e: dict with t_a, d, clearance, sun_ok, maint_ok, ok."""
    ta, d, ps, pr = light_time(sender, receiver, te_h)
    clr = seg_clearance(ps, pr)
    sun_ok = clr >= R_EXCL
    te = np.asarray(te_h, dtype=float)
    m_ok = np.ones_like(sun_ok, dtype=bool)
    if maintenance:
        for (s0, s1) in MAINT.get(frozenset((sender, receiver)), []):
            # closed flight interval [t_e, t_a] overlaps half-open window [s0, s1)
            m_ok &= ~((te < s1) & (ta >= s0))
    return dict(te=te, ta=ta, d=d, clr=clr, sun_ok=sun_ok, maint_ok=m_ok, ok=sun_ok & m_ok)


def loss_backbone(d):
    return 1.0 - np.exp(-0.02 * np.asarray(d))


def loss_direct(d):
    return 1.0 - np.exp(-0.08 * np.asarray(d))


# ---------------------------------------------------------------- waiting / boundaries
def next_launch(sender, receiver, t_h, horizon_h=24 * 400, step_h=1.0 / 60, maintenance=True):
    """Earliest t_e >= t_h at which the link can launch (geometry + known maintenance).
    Coarse 1-minute march, then bisection to 1e-7 h (0.36 ms). None if nothing in horizon."""
    if link_state(sender, receiver, t_h, maintenance)["ok"]:
        return float(t_h)
    chunk = 24 * 30
    t0 = t_h
    while t0 - t_h < horizon_h:
        grid = t0 + np.arange(0, chunk, step_h)
        ok = link_state(sender, receiver, grid, maintenance)["ok"]
        idx = np.argmax(ok)
        if ok[idx]:
            lo, hi = (grid[idx - 1] if idx > 0 else t0), grid[idx]
            while hi - lo > 1e-7:
                mid = 0.5 * (lo + hi)
                if link_state(sender, receiver, mid, maintenance)["ok"]:
                    hi = mid
                else:
                    lo = mid
            return float(hi)
        t0 = grid[-1] + step_h
    return None


# ---------------------------------------------------------------- routes
def candidate_routes(x, y):
    """All simple backbone routes of at most 3 links from gateway x to gateway y."""
    if x == y:
        return []
    return [(x, "Relay A", y), (x, "Relay B", y),
            (x, "Relay A", "Relay B", y), (x, "Relay B", "Relay A", y)]


def route_timing(route, t_ready_h, wait=True, maintenance=True):
    """Propagate one packet along a pinned backbone route, no loss, empty queues.
    Per hop: t_e = ready + 1 s (serialization); if wait and the link cannot launch then,
    hold until next_launch. At a relay, ready = t_a + 1 s. Returns dict with hops and arrival.
    Backbone launches start at the gateway, so no local-access second is added."""
    ready = float(t_ready_h)
    hops = []
    for k in range(len(route) - 1):
        s, r = route[k], route[k + 1]
        te = ready + SEC
        st = link_state(s, r, te, maintenance)
        held = 0.0
        if not st["ok"]:
            if not wait:
                return None
            nt = next_launch(s, r, te, maintenance=maintenance)
            if nt is None:
                return None
            held, te = nt - te, nt
            st = link_state(s, r, te, maintenance)
        ta = float(st["ta"])
        hops.append(dict(link=(s, r), ready=ready, te=te, ta=ta, d=float(st["d"]),
                         clr=float(st["clr"]), held_h=held, p_loss=float(loss_backbone(st["d"]))))
        ready = ta + SEC if r in RELAYS else ta
    return dict(route=route, t_ready=float(t_ready_h), arrival=hops[-1]["ta"],
                delay_h=hops[-1]["ta"] - t_ready_h, hops=hops)


def T0(route, t_enq_h):
    """Brief's predicted empty-queue one-way delay at enqueue: flight + serialization +
    relay processing, blocked links treated as open (pure light time)."""
    ready, t = float(t_enq_h), float(t_enq_h)
    for k in range(len(route) - 1):
        te = ready + SEC
        ta = float(light_time(route[k], route[k + 1], te)[0])
        ready = ta + SEC if route[k + 1] in RELAYS else ta
    return ta - t


def best_route(x, y, t_h, wait=True, maintenance=True):
    best = None
    for rt in candidate_routes(x, y):
        r = route_timing(rt, t_h, wait=wait, maintenance=maintenance)
        if r is not None and (best is None or r["arrival"] < best["arrival"]):
            best = r
    return best


def route_open_now(route, t_h, maintenance=True):
    """True if every link on the route can launch at emission time t_h (S3 availability rule)."""
    return all(bool(link_state(route[k], route[k + 1], t_h, maintenance)["ok"])
               for k in range(len(route) - 1))


def route_availability(route, h0, span_h=24.0, step_h=1.0 / 60, maintenance=True):
    """Fraction of [h0, h0+span) in which every link on the route can launch (1-minute grid)."""
    grid = h0 + np.arange(0, span_h, step_h)
    ok = np.ones_like(grid, dtype=bool)
    for k in range(len(route) - 1):
        ok &= link_state(route[k], route[k + 1], grid, maintenance)["ok"]
    return float(ok.mean())


# ---------------------------------------------------------------- probabilities (E2)
def hop_abandon_prob(p):
    """A hop is abandoned only if all 4 launches are lost: p^4 (independent launches)."""
    return p ** 4


def hop_receipt_fail_prob(p_fwd, p_rev):
    """Prob. a single launch+receipt round fails (data lost OR receipt lost) -> sender retries."""
    return 1 - (1 - p_fwd) * (1 - p_rev)


def route_first_try(ps):
    return float(np.prod([1 - p for p in ps]))


def route_abandon_bound(ps):
    """P(some hop abandoned) <= sum p_i^4 (union bound); exact = 1 - prod(1 - p_i^4)."""
    exact = 1 - np.prod([1 - p ** 4 for p in ps])
    return float(exact), float(sum(p ** 4 for p in ps))


def direct_success(p, n_copies):
    """At least one of n independent direct copies arrives: 1 - p^n."""
    return 1 - p ** n_copies


def fmt_h(h):
    sgn = "-" if h < 0 else ""
    h = abs(h)
    H = int(h); m = (h - H) * 60; M = int(m); s = (m - M) * 60
    return f"{sgn}{H}h{M:02d}m{s:06.3f}s"
