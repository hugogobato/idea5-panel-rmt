"""WP-B3 pilot: one-seed-class cost measurement + LOCAL/COLAB classification.

Measures wall time and RSS per replication on a production-size decisive cell
(WP-C1 mid-grid point) plus a size-scaling probe, fits a log-log cost model,
and prints the LOCAL/COLAB classification table for the planned Phase C
experiments. Run context (system load) is recorded alongside.
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import psutil

import scm_frontier as sf

METHODS = ("donor_mean", "scm_simplex", "ridge_sc", "spectral_sc", "mc_nn_cv", "sdid")
REPS_PILOT = 10


def rep_cost(donors, y1_pre, y_star):
    row = {}
    rss_peak = 0.0
    proc = psutil.Process()
    for m in METHODS:
        t0 = time.perf_counter()
        pred = getattr(sf, m)(donors, y1_pre)
        dt = time.perf_counter() - t0
        err = float(np.sqrt(np.mean((pred - y_star) ** 2)))
        rss_peak = max(rss_peak, proc.memory_info().rss / 2**30)
        row[m] = {"wall_s": dt, "rmse": err}
    return row, rss_peak


def pilot_cell(n, T0, T_post, s, theta, reps):
    per_method = {m: [] for m in METHODS}
    rmses = {m: [] for m in METHODS}
    rss_max = 0.0
    for j in range(reps):
        pan = sf.generate_panel(n=n, T0=T0, T_post=T_post, r=1,
                                spike_strengths=(s,), treated_share=(theta,),
                                alignment="first", sigma=1.0, seed=20_000 + j)
        Y = pan["Y"]
        donors, y1_pre, y_star = Y[1:], Y[0, :T0], Y[0, T0:]
        row, rss = rep_cost(donors, y1_pre, y_star)
        rss_max = max(rss_max, rss)
        for m in METHODS:
            per_method[m].append(row[m]["wall_s"])
            rmses[m].append(row[m]["rmse"])
    return {
        "per_method_wall_s": {m: list(map(float, v)) for m, v in per_method.items()},
        "per_method_rmse_mean": {m: float(np.mean(v)) for m, v in rmses.items()},
        "rss_gib": float(rss_max),
        "wall_total_s": float(sum(sum(v) for v in per_method.values())),
    }


def main():
    load1 = psutil.getloadavg()
    print(f"system load (1/5/15min): {load1[0]:.1f} {load1[1]:.1f} {load1[2]:.1f}")
    print(f"cpu_count: {psutil.cpu_count()} | mem total {psutil.virtual_memory().total / 2**30:.1f} GiB")

    print("\n[pilot] decisive cell n=250 T0=250 T_post=125 s=2.0 theta=0.5 "
          f"(c=1.0), {REPS_PILOT} reps")
    res = pilot_cell(250, 250, 125, 2.0, 0.5, REPS_PILOT)
    print(f"{'method':<14}{'mean_ms':>9}{'max_ms':>9}{'RMSE':>9}")
    for m in METHODS:
        v = res["per_method_wall_s"][m]
        print(f"{m:<14}{1000*np.mean(v):>9.0f}{1000*max(v):>9.0f}"
              f"{res['per_method_rmse_mean'][m]:>9.3f}")
    print(f"peak RSS across reps: {res['rss_gib']:.3f} GiB")

    print("\n[scaling probe] spectral/scm/sdid/mc wall seconds, 2 reps each")
    sizes = [(125, 125, 63), (250, 250, 125), (400, 400, 200)]
    probe = {}
    for (n, T0, Tp) in sizes:
        pan = sf.generate_panel(n=n, T0=T0, T_post=Tp, r=1, spike_strengths=(2.0,),
                                treated_share=(0.5,), alignment="first", seed=31_000)
        Y = pan["Y"]
        donors, y1_pre = Y[1:], Y[0, :T0]
        entry = {}
        for m in ("spectral_sc", "scm_simplex", "sdid", "mc_nn_cv"):
            ts = []
            for _ in range(2):
                t0 = time.perf_counter()
                getattr(sf, m)(donors, y1_pre)
                ts.append(time.perf_counter() - t0)
            entry[m] = float(np.mean(ts))
        probe[f"{n}x{T0}"] = entry
        print(f"n={n:4d}: " + " ".join(f"{m}={entry[m]:.2f}s" for m in entry))

    # log-log scaling exponent vs n*T0 (use largest two sizes)
    keys = list(probe)
    x0 = np.log(np.prod([int(k) for k in keys[-2].split("x")]))
    x1 = np.log(np.prod([int(k) for k in keys[-1].split("x")]))
    slopes = {}
    for m in ("spectral_sc", "scm_simplex", "sdid", "mc_nn_cv"):
        y0, y1v = np.log(max(probe[keys[-2]][m], 1e-4)), np.log(max(probe[keys[-1]][m], 1e-4))
        slopes[m] = float((y1v - y0) / (x1 - x0))
    print("log-log exponents vs n*T0:", {m: round(v, 2) for m, v in slopes.items()})

    # classification table for planned experiments (workers=8 cap, headroom rule)
    mean_rep = {m: float(np.mean(res["per_method_wall_s"][m])) for m in METHODS}
    sum_rep = sum(mean_rep.values())
    cells_c1 = 5 * 15 * 3 * 2
    experiments = [
        ("C1 decisive grid (450 cells x 500 reps)", cells_c1 * 500 * sum_rep),
        ("C1 single-c slice (15 x 500)", 15 * 500 * sum_rep),
        ("C2(i) null battery (6 cells x 500)", 6 * 500 * sum_rep),
        ("C2(ii) baseline-favorable (4 x 500)", 4 * 500 * sum_rep),
        ("C2(iv) structural-break (4 x 500)", 4 * 500 * sum_rep),
        ("C2(iii) AR(1)+bootstrap (2 rho x 500 x boot~50x diag cost)",
         2 * 500 * (sum_rep + 50 * (mean_rep["spectral_sc"] + mean_rep["ridge_sc"]))),
    ]
    print("\nLOCAL/COLAB classification (rule: LOCAL iff projected < 2 h AND "
          "< 4 GiB; workers <= 8)")
    print(f"{'experiment':<58}{'serial_h':>9}{'8wkr_h':>9}  class")
    for name, secs in experiments:
        h_serial = secs / 3600
        h_workers = h_serial / 8
        cls = "LOCAL" if (h_workers < 2 and res["rss_gib"] < 4) else "COLAB"
        print(f"{name:<58}{h_serial:>9.2f}{h_workers:>9.2f}  {cls}")

    import json
    out = Path(__file__).resolve().parents[1] / "figures"
    out.mkdir(exist_ok=True)
    payload = {"load": load1, "pilot": res, "probe": probe, "slopes": slopes,
               "sum_rep_serial_s": sum_rep}
    (out / "pilot_costs.json").write_text(json.dumps(payload, indent=1))
    print(f"\nsaved figures/pilot_costs.json")


if __name__ == "__main__":
    main()
