"""WP-B2 smoke run: 20 reps x one cell x all methods, timed.

Pass target (plan WP-B2): completes in under 1 minute on an idle laptop;
this script prints per-method wall time so deviations are attributable.
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import scm_frontier as sf

REPS = 20
CFG = dict(n=81, T0=160, T_post=60, r=1, spike_strengths=(6.0,),
           treated_share=(1.0,), alignment="first")
METHODS = ("donor_mean", "scm_simplex", "ridge_sc", "spectral_sc", "mc_nn_cv", "sdid")


def main():
    t_start = time.perf_counter()
    rmse = {m: [] for m in METHODS}
    times = {m: 0.0 for m in METHODS}
    for j in range(REPS):
        pan = sf.generate_panel(**CFG, seed=7000 + j)
        donors, y1_pre = pan["Y"][1:], pan["Y"][0, : CFG["T0"]]
        y_star = pan["Y"][0, CFG["T0"] :]
        for m in METHODS:
            t0 = time.perf_counter()
            pred = getattr(sf, m)(donors, y1_pre)
            times[m] += time.perf_counter() - t0
            rmse[m].append(np.sqrt(np.mean((pred - y_star) ** 2)))
    total = time.perf_counter() - t_start

    print(f"smoke cell: {CFG}, reps={REPS}")
    print(f"{'method':<14}{'wall_s':>9}{'ms/rep':>9}{'RMSE':>9}")
    for m in METHODS:
        print(f"{m:<14}{times[m]:>9.2f}{1000 * times[m] / REPS:>9.0f}"
              f"{np.mean(rmse[m]):>9.3f}")
    print(f"TOTAL {total:.2f}s ({'PASS' if total < 60 else 'OVER'} vs 60s target)")


if __name__ == "__main__":
    main()
