"""C2(iv) repair analysis (post-hoc deviation, documented in gate_g3_memo.md).

The preregistration froze Z_boot/Z_tw as pseudo-cutoff statistics computed
ONLY inside the observed pre window (valid for size calibration in C2(iii),
vacuous for break DETECTION in C2(iv), whose DGPs differ only in the post
window). This script recomputes the correctly-windowed statistic:

    z_post = resid_statistic(donor PRE basis, REAL donor POST window)

on bitwise-reproduced C2(iv) panels (seeds 10000+j, prereg Section 4.4),
calibrated by the same simulated finite-n null. Uses only donor pre+post
outcomes (leakage-safe; treated post never touched).
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from scm_frontier.diagnostics import resid_statistic, simulated_null_z, z_tw_pvalue  # noqa
from scm_frontier.dgps import generate_panel  # noqa


def main() -> dict:
    n, T0, Tp, r, m = 160, 160, 80, 1, 2.0
    c = n / T0
    s = m * np.sqrt(c)
    cells = [("control_d0_aligned", None, "first", 1.0),
             ("break_d1_aligned", 1.0, "first", 1.0),
             ("break_d2_aligned", 2.0, "first", 1.0),
             ("break_d2_orthogonal", 2.0, "none", 0.0)]
    null_z = simulated_null_z(n - 1, T0, Tp, 1.0, (n - 1) / T0, G=300,
                              seed=7_770_001)
    out = {}
    for cid, delta, align, th in cells:
        ps = []
        for j in range(500):
            pan = generate_panel(seed=10_000 + j, n=n, T0=T0, T_post=Tp, r=r,
                                 spike_strengths=(s,), treated_share=(th,),
                                 alignment=align, sigma=1.0,
                                 structural_break=delta)
            Y = pan["Y"]
            z, _ = resid_statistic(Y[1:, :T0], Y[1:, T0:], 1.0, (n - 1) / T0)
            ps.append(z_tw_pvalue(z, null_z))
        ps = np.array(ps)
        out[cid] = {"rej_5": float((ps < 0.05).mean()),
                    "rej_10": float((ps < 0.10).mean())}
        print(cid, out[cid])
    dest = ROOT / "figures" / "c2iv_repair_analysis.json"
    dest.write_text(json.dumps(out, indent=1))
    print("saved", dest)
    return out


if __name__ == "__main__":
    import json
    main()
