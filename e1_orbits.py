"""E1 orbital evidence. Writes out/e1_tables.md."""
import numpy as np, os, math
import mpx

REF = {  # Brief Section 3, epoch check positions (AU)
    "Mercury": (-0.213426, -0.410281, -0.013955), "Venus": (0.679067, -0.257950, -0.042725),
    "Earth": (1.003581, -0.023693, -0.000003), "Mars": (0.247476, 1.526869, 0.025929),
    "Ceres": (0.375786, 2.653941, 0.014791), "Jupiter": (-3.438179, 4.038309, 0.060149),
    "Saturn": (9.271221, 1.717485, -0.398962), "Uranus": (8.962475, 17.253919, -0.052129),
    "Neptune": (29.839098, 1.351785, -0.715470), "Relay A": (2.0, 2.0, 0.0), "Relay B": (-2.0, 2.0, 0.0)}
BODIES = list(REF)
L = ["# E1 orbital calculations", "",
     "Method: independent fixed Kepler ellipses from orbital_elements.json / network_model.json. "
     "M(t) = M0 + n t (days) mod 360; Kepler's equation by Newton to |dE| < 1e-15 rad; "
     "orbit-plane (x', y') rotated by omega, i, Omega into the Sun-centred J2000 ecliptic frame (AU). "
     "Time: hours since 2026-09-22 00:00 TDB, 86,400 s days. Light time solved by fixed-point "
     "iteration on t_a = t_e + 8.317 min/AU * |r_recv(t_a) - r_send(t_e)| until successive "
     "iterates differ by < 1e-9 h (3.6 us), 280x tighter than the brief's 1 ms.", "",
     "## Epoch check (t = 0)", "", "| Body | x | y | z | max abs error vs brief (AU) | Pass (<1e-5) |",
     "|---|---|---|---|---|---|"]
worst = 0
for b in BODIES:
    p = mpx.pos(b, 0.0)
    err = float(np.max(np.abs(p - np.array(REF[b]))))
    worst = max(worst, err)
    L.append(f"| {b} | {p[0]:.6f} | {p[1]:.6f} | {p[2]:.6f} | {err:.2e} | {'yes' if err < 1e-5 else 'NO'} |")
L += ["", f"Worst error: {worst:.2e} AU (limit 1e-5). Errors are at the 6-decimal rounding of the brief's table.", ""]

L += ["## Period return and radius bounds", "",
      "| Body | Period (d) | abs(r(P) - r(0)) (AU) | r min sampled | a(1-e) | r max sampled | a(1+e) |",
      "|---|---|---|---|---|---|---|"]
for b in BODIES:
    el = mpx.elements(b)
    P = el["period_days"] * 24
    ret = float(np.linalg.norm(mpx.pos(b, P) - mpx.pos(b, 0.0)))
    t = np.linspace(0, P, 200001)
    r = np.linalg.norm(mpx.pos(b, t), axis=-1)
    a, e = el["a_au"], el["e"]
    L.append(f"| {b} | {el['period_days']:.4f} | {ret:.1e} | {r.min():.6f} | {a*(1-e):.6f} | "
             f"{r.max():.6f} | {a*(1+e):.6f} |")
L.append("")

L += ["## Worked positions at later times (AU)", "", "| Body | h 300 | +1 yr | +10 yr | +100 yr |",
      "|---|---|---|---|---|"]
for b in BODIES:
    cells = []
    for t in (300, 8766, 87660, 876600):
        p = mpx.pos(b, t)
        cells.append(f"({p[0]:.5f}, {p[1]:.5f}, {p[2]:.5f})")
    L.append(f"| {b} | " + " | ".join(cells) + " |")
L.append("")


def receiver_motion(s, r, te):
    """Arrival with moving receiver vs. receiver frozen at its emission-time position."""
    ps = mpx.pos(s, te)
    frozen = te + mpx.LT_H * float(np.linalg.norm(mpx.pos(r, te) - ps))
    iters, ta = [frozen], frozen
    for _ in range(6):
        ta = te + mpx.LT_H * float(np.linalg.norm(mpx.pos(r, ta) - ps))
        iters.append(ta)
    solved = float(mpx.light_time(s, r, te)[0])
    return frozen, solved, iters


L += ["## Receiver motion changes arrival time", ""]
for (s, r) in [("Neptune", "Relay A"), ("Earth", "Relay A"), ("Mercury", "Relay B")]:
    te = mpx.SEC
    frozen, solved, it = receiver_motion(s, r, te)
    L += [f"**{s} -> {r}, t_e = 1 s:** frozen-receiver arrival {mpx.fmt_h(frozen)}, "
          f"moving-receiver arrival {mpx.fmt_h(solved)}, difference **{(solved-frozen)*3600:+.3f} s** "
          f"({abs(solved-frozen)*3.6e6:.0f}x the 1 ms tolerance).", "",
          "Fixed-point iterates (s after t_e): " +
          ", ".join(f"{(x-te)*3600:.6f}" for x in it[:5]), ""]

L += ["## Forward and return analysed separately", "",
      "Round trip Earth gateway -> Relay A -> Earth at t_ready = 0: forward packet, relay 1 s processing "
      "+ 1 s serialization, return launch from the relay's arrival-time position.", "",
      "| Leg | t_e | t_a | flight (s) | d (AU) | Sun clearance (AU) |", "|---|---|---|---|---|---|"]
te1 = mpx.SEC
f = mpx.link_state("Earth", "Relay A", te1)
te2 = float(f["ta"]) + 2 * mpx.SEC
g = mpx.link_state("Relay A", "Earth", te2)
for name, st, te in (("Earth -> Relay A", f, te1), ("Relay A -> Earth", g, te2)):
    L.append(f"| {name} | {mpx.fmt_h(te)} | {mpx.fmt_h(float(st['ta']))} | "
             f"{(float(st['ta'])-te)*3600:.3f} | {float(st['d']):.6f} | {float(st['clr']):.4f} |")
L += ["", f"Forward and return flights differ by {abs((float(f['ta'])-te1)-(float(g['ta'])-te2))*3600:.3f} s "
      "because both endpoints move during the round trip; a symmetric-path assumption would be wrong by this much.", ""]

open(os.path.join(mpx.HERE, "out", "e1_tables.md"), "w").write("\n".join(L))
print("\n".join(L))
