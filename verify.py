"""Independent cross-checks of the toolkit. Every check prints PASS/FAIL."""
import numpy as np, json, os, math, random
from scipy.optimize import brentq
import mpx

ok_all = True
def check(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")

# 1. light time vs. brentq root of f(ta) = ta - te - LT*|r(ta)-s(te)|
rng = np.random.default_rng(1)
worst = 0
for _ in range(200):
    s, r = mpx.DIRECTED[rng.integers(len(mpx.DIRECTED))]
    te = float(rng.uniform(0, 200 * 8766))
    ps = mpx.pos(s, te)
    f = lambda ta: ta - te - mpx.LT_H * float(np.linalg.norm(mpx.pos(r, ta) - ps))
    root = brentq(f, te, te + 10, xtol=1e-12)
    worst = max(worst, abs(root - float(mpx.light_time(s, r, te)[0])))
check("light time matches brentq root (200 random links/times)", worst * 3600 < 1e-3, f"worst {worst*3.6e6:.4f} ms")

# 2. closure boundaries: open just outside, closed just inside (1 s)
cl = json.load(open(os.path.join(mpx.HERE, "out", "e4_closures.json")))
bad = n = 0
for x in cl["links"]:
    s, r = x["link"]
    for (a, b, _) in random.Random(0).sample(x["closures"], min(5, len(x["closures"]))):
        n += 1
        c = lambda t: float(mpx.link_state(s, r, t, maintenance=False)["clr"])
        if a > 1e-3 and not (c(a - mpx.SEC) >= 0.1 and c(a + mpx.SEC) < 0.1): bad += 1
        if not (c(b + mpx.SEC) >= 0.1 and c(b - mpx.SEC) < 0.1): bad += 1
check("closure boundaries bracket the 0.10 AU crossing to 1 s", bad == 0, f"{n} closures sampled, {bad} bad")

# 3. Lipschitz bound respected numerically (finite differences, 1-minute)
worst_ratio = 0
for (s, r) in mpx.DIRECTED:
    t = np.sort(rng.uniform(0, 200 * 8766, 20000))
    h = 1 / 60
    c1 = mpx.link_state(s, r, t, maintenance=False)["clr"]
    c2 = mpx.link_state(s, r, t + h, maintenance=False)["clr"]
    rate = np.max(np.abs(c2 - c1)) / h
    V = max(mpx.max_speed_au_per_h(s), mpx.max_speed_au_per_h(r)) * 1.001
    worst_ratio = max(worst_ratio, rate / V)
check("observed |d clr/dt| <= V bound on every link", worst_ratio <= 1.0, f"max observed/bound = {worst_ratio:.3f}")

# 4. Two-Relay Lemma analytic margin: half-width of the blocked cone vs. 45 deg
r_min = min(mpx.elements(b)["a_au"] * (1 - mpx.elements(b)["e"]) for b in mpx.SETTLEMENTS)
R = mpx.elements("Relay A")["a_au"]
half = math.degrees(math.asin(0.1 * (1 / r_min + 1 / R)))
check("Two-Relay Lemma: blocked half-angle < 45 deg (relays 90 deg apart)", half < 45,
      f"half-angle <= {half:.2f} deg (Mercury perihelion)")

# 5. Relay A-B clearance constant ~2 AU
t = rng.uniform(0, 200 * 8766, 5000)
c = mpx.link_state("Relay A", "Relay B", t, maintenance=False)["clr"]
check("Relay A<->B clearance stays ~2 AU", c.min() > 1.99, f"min {c.min():.5f}, max {c.max():.5f}")

# 6. Monte Carlo of Neptune->Earth single-attempt delivery CDF vs exact
import importlib.util, sys
rt = mpx.best_route("Neptune", "Earth", 0.0)
hops = rt["hops"]
Rh = [2 * (h["ta"] - h["te"]) + 1 for h in hops]
M = 400_000
extra = np.zeros(M); fail = np.zeros(M, bool)
for h, R_ in zip(hops, Rh):
    k = rng.geometric(1 - h["p_loss"], M) - 1      # failed launches before first success
    fail |= k >= 4
    extra += np.minimum(k, 3) * R_
exact_never = 1 - np.prod([1 - h["p_loss"] ** 4 for h in hops])
check("Monte Carlo route abandonment matches exact", abs(fail.mean() - exact_never) < 4 * math.sqrt(exact_never / M),
      f"MC {fail.mean():.5f} vs exact {exact_never:.5f}")
T = rt["delay_h"] + Rh[0]
mc = np.mean(~fail & (rt["delay_h"] + extra <= T))
spec = importlib.util.spec_from_file_location("s3", os.path.join(mpx.HERE, "s3_e2_e5.py"))
src = open(os.path.join(mpx.HERE, "s3_e2_e5.py")).read().split("# ------------------------------------------------------------------ S3")[0]
ns = {}; exec(src, ns)
ex = ns["route_delivery_cdf"](rt, [T])[0]
check("Monte Carlo P(delivered by cond + R_h1) matches exact", abs(mc - ex) < 0.003, f"MC {mc:.4f} vs exact {ex:.4f}")

# 7. maintenance rule: B->Neptune launch whose flight overlaps [2,26) must fail; just after must pass
st_in = mpx.link_state("Relay B", "Neptune", 1.0)
st_after = mpx.link_state("Relay B", "Neptune", 26.0)
check("maintenance overlap rule (B->Neptune at h1 blocked, at h26 open)",
      (not bool(st_in["ok"])) and bool(st_after["ok"]))

print("\nALL PASS" if ok_all else "\nSOME CHECKS FAILED")
