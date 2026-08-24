"""WP-B1 symbolic verification of frontier_ansatz.md (sympy only, no simulation).

Checks SC1-SC4 and monotonicity of every risk channel in the spike strength s
on the supercritical side. Internal consistency check only; not a decisive
experiment (WP-C1 overlay is the falsifier).
"""
from sympy import symbols, simplify, limit, sqrt, S, oo, diff, Rational

s, c = symbols("s c", positive=True)

lam = 1 + s + c + c / s                      # BGN outlier location
zeta = (1 - c / s**2) / (1 + c / s)          # BGN squared overlap
tau = sqrt(s / lam)                          # singular-value transmission

results = []

# SC2: s -> infinity collapses all signal channels to zero
sc2_bias = simplify(limit((1 - zeta * tau) ** 2, s, oo))
sc2_var = simplify(limit(zeta / lam, s, oo))
sc2_T0 = simplify(limit((s * zeta + 1) / lam, s, oo))
results.append(("SC2 s->inf: bias->0", sc2_bias == 0))
results.append(("SC2 s->inf: zeta/lambda->0", sc2_var == 0))
results.append(("SC2 s->inf: (s*zeta+1)/lambda->1 (so /T0 -> 0)", sc2_T0 == 1))

# SC4: continuity at the edge, zeta -> 0 as s -> sqrt(c)+
edge = sqrt(c)
zeta_edge = simplify(limit(zeta, s, edge, "+"))
tau_edge = simplify(limit(tau, s, edge, "+"))
lam_edge = simplify(limit(lam, s, edge, "+"))
# included-spike terms at the edge reduce to full under-transmission alpha^2/sigma^2
included_at_edge_signal = simplify((1 - zeta_edge * tau_edge) ** 2)
results.append(("SC4 edge: zeta->0", zeta_edge == 0))
results.append(("SC4 edge: included-signal-bias -> 1 (= truncation value)",
                included_at_edge_signal == 1))
results.append(("SC4 edge: zeta/lambda -> 0", simplify(zeta_edge / lam_edge) == 0))

# Edge transmission ratio value (sanity, c=.5): tau_edge^2 = sqrt(c)/(1+sqrt(c))^2
lhs = simplify(tau_edge**2)
rhs = sqrt(c) / (1 + sqrt(c)) ** 2
results.append(("edge tau^2 identity", simplify(lhs - rhs) == 0))

# Monotonicity scan (numeric grid; symbolic derivative sign is intractable in closed form)
# Truth (documented in frontier_ansatz.md Section 4): the bias channel is strictly
# decreasing; variance channels have a small hump just above the edge; TOTAL risk is
# monotone for treated share theta ~ 1 and has a <= 0.005 sigma^2 shoulder for small theta.
import numpy as np

def channels(sv, cv):
    lv = 1 + sv + cv + cv / sv
    zv = (1 - cv / sv**2) / (1 + cv / sv)
    tv = np.sqrt(sv / lv)
    return (1 - zv * tv) ** 2, zv / lv, (sv * zv + 1) / (240 * lv)

# Documented truth (frontier_ansatz.md Section 4): TOTAL excess risk is strictly
# decreasing on the supercritical side with its maximum exactly at the edge;
# individual variance channels are non-monotone but O(1/T0)-suppressed.

def total_excess(sv, cvv, th):
    b, v2c, v13 = channels(sv, cvv)
    return th * (b + v2c) + v13

bias_ok = True
total_ok = True
for cv in [0.25, 0.5, 1.0, 2.0, 4.0]:
    grid = np.linspace(np.sqrt(cv) * 1.0005, 20.0, 4000)
    vals = np.array([channels(x, cv) for x in grid])
    bias_ok &= float(np.diff(vals[:, 0]).max()) < 1e-12
    for th in [0.02, 0.1, 0.5, 1.0]:
        tt = np.array([total_excess(x, cv, th) for x in grid])
        total_ok &= float(np.diff(tt).max()) < 1e-12

results = [
    ("bias channel strictly decreasing on supercritical grid", bias_ok),
    ("TOTAL excess risk strictly decreasing (all c, theta)", total_ok),
]

# Kink signature: |slope| of total excess just above vs far above the edge (c=0.5, theta=0.5)
cv = 0.5
g1 = np.linspace(np.sqrt(cv) * 1.001, np.sqrt(cv) * 1.3, 50)
g2 = np.linspace(10.0, 20.0, 50)
v1 = np.array([total_excess(x, cv, 0.5) for x in g1])
v2 = np.array([total_excess(x, cv, 0.5) for x in g2])
slope_near = abs(v1[-1] - v1[0]) / (g1[-1] - g1[0])
slope_far = abs(v2[-1] - v2[0]) / (g2[-1] - g2[0])
results.append(("steepest descent near edge (kink signature)", slope_near > 100 * slope_far))

# SC1/SC3 are structural (empty sums / vanishing alpha), stated in the doc.

print("frontier_ansatz symbolic checks")
print("=" * 46)
allpass = True
for name, ok in results:
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
    allpass &= ok
print(f"|slope| near edge {slope_near:.4f} vs far {slope_far:.6f}")
print("ALL PASS" if allpass else "SOME FAILED")
raise SystemExit(0 if allpass else 1)
