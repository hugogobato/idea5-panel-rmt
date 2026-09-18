"""WP-D2b Wave-2 driver (docs/preregistrations/preregistration_d2b_addendum.md, frozen pre-run).

W1: Basque panel through the unchanged Wave-1 lens.
W2: donor-pool degradation battery on smoking (positive-sensitivity
control, lite labels, bootstrap disabled per deviation D-W2).

Writes results_d2/wave2_*.json|csv and figures/fig_d2_wave2.png. Gate
outcomes are recorded in the artifacts and judged in the wave-2 memo;
nothing is tuned here.
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from applications import frontier_d2 as F  # noqa: E402
from applications import loaders as L  # noqa: E402
from applications import synth_adh as S  # noqa: E402

SEEDS = {"align": 60_101, "donor": 60_102, "time": 60_103, "perm": 60_104}
Z_SHIFT_SEED = 8_880_001
POSTTEST_SEED = 7_770_001
P1_THRESH_BASQUE = -0.45
DEG_SIZES = [38, 24, 16, 12, 8, 6, 4, 3]
DEG_DRAWS = 20
OUT = ROOT / "results_d2"
FIG = ROOT / "figures"


def _clean(obj):
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


def outcome_only_scm_pre(Yd_pre: np.ndarray, y1_pre: np.ndarray):
    """Frozen Wave-2 spec: outcome-only classic SCM pre-fit (deviation D-W1)."""
    w = S._nnls_simplex_ridge(Yd_pre.T, y1_pre, ridge2=0.0, center=True)
    path = w @ Yd_pre
    rmse = float(np.sqrt(np.mean((y1_pre - path) ** 2)))
    return rmse


def basque_windows(panel: L.Panel):
    years = panel.years
    i61 = int(np.where(years == 1961)[0][0])
    T0 = int(np.where(years == panel.treat_year)[0][0])
    return slice(i61, T0), slice(T0, len(years))


def run_basque():
    print("=== W1: Basque ===", flush=True)
    panel = L.load_basque()
    pre_sl, post_sl = basque_windows(panel)
    don = np.arange(len(panel.units)) != panel.treated_idx
    Yd_pre, y1_pre = panel.Y[don][:, pre_sl], panel.Y[panel.treated_idx, pre_sl]
    cutoff_idx = [8, 10]
    t = F.analyze_unit(Yd_pre, y1_pre, SEEDS, idx=0)
    t["date_alarms"] = F.date_alarm_scan(Yd_pre, cutoff_idx,
                                         Z_SHIFT_SEED, t["sigma2_1"])
    t["alarm_rate"] = float(np.mean([a["alarm"] for a in t["date_alarms"]]))
    t["P3_pass"] = bool(t["alarm_rate"] <= F.ALARM_MAX)
    t["posttest_control"] = F.posttest_control(
        Yd_pre, panel.Y[don][:, post_sl], t["sigma2_1"], POSTTEST_SEED)

    units = [u for u in panel.units if u != L.BASQUE_TREATED]
    rows = []
    for i, name in enumerate(units):
        keep = [u for u in panel.units
                if u not in (name, L.BASQUE_TREATED)]
        m = np.array([panel.units.index(u) for u in keep])
        i_name = panel.units.index(name)
        Yn_pre = panel.Y[m][:, pre_sl]
        y1_i = panel.Y[i_name, pre_sl]
        rmse = outcome_only_scm_pre(Yn_pre, y1_i)
        diag = F.analyze_unit(Yn_pre, y1_i, SEEDS, idx=i + 1)
        rows.append({"unit": name, "rmse": rmse,
                     **{k: diag[k] for k in
                        ["k", "d", "p_align", "r_align", "label",
                         "boot_donor", "stability_pass"]}})
        print(f"  placebo {i + 1}/{len(units)} {name}: d={diag['d']:.3f} "
              f"k={diag['k']} p={diag['p_align']:.3f} "
              f"rmse={rmse:.3f} {diag['label']}", flush=True)

    treated_rmse = outcome_only_scm_pre(Yd_pre, y1_pre)
    q75 = float(np.percentile([r["rmse"] for r in rows], 75))
    d_vec = np.array([r["d"] for r in rows])
    rmse_vec = np.array([r["rmse"] for r in rows])
    labels = [r["label"] for r in rows]
    rho, perm_p = F.spearman_with_perm(d_vec, rmse_vec, SEEDS["perm"])
    nf = F.nf_arms(t, treated_rmse, [r["unit"] for r in rows],
                   d_vec, rmse_vec, labels)
    nf["P1_fit_ordering"] = {"rho": rho, "perm_p": perm_p,
                             "threshold": P1_THRESH_BASQUE,
                             "pass": bool(rho <= P1_THRESH_BASQUE
                                          and perm_p <= 0.05)}
    nf["P2_treated_separation"] = bool(t["d"] >= float(
        np.percentile(d_vec, 90))) if len(d_vec) else False
    nf["P2_q90_placebo_d"] = float(np.percentile(d_vec, 90))
    nf["P3_date_alarms_pass"] = t["P3_pass"]
    nf["stability_pass"] = t["stability_pass"]
    fragile_ok = ("FRAGILE-MISALIGNED" == t["label"]
                  or t["label"] in ("INCONCLUSIVE-SUBEDGE",
                                    "FRAGILE-SUBEDGE"))
    nf["W1_PASS"] = bool(fragile_ok and t["P3_pass"] and t["stability_pass"])
    out = {"treated_label": t["label"], "treated_d": t["d"],
           "treated_boot_donor": t["boot_donor"],
           "treated_p_align": t["p_align"],
           "treated_pre_rmse": treated_rmse, "placebo_rmse_q75": q75,
           "nf": nf, "treated": _clean(t)}
    with open(OUT / "wave2_basque.json", "w") as fh:
        json.dump(_clean(out), fh, indent=2)
    keys = ["unit", "rmse", "k", "d", "p_align", "label"]
    with open(OUT / "wave2_basque_placebos.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows:
            w.writerow([r["unit"], r["rmse"], r["k"], r["d"],
                        r["p_align"], r["label"]])
    return out, rows, t


def run_degradation():
    print("=== W2: degradation battery (smoking) ===", flush=True)
    panel_a, _ = L.load_smoking()
    Ta = panel_a.T0
    don = np.arange(len(panel_a.units)) != panel_a.treated_idx
    Yc = F.center_rows(panel_a.Y[don][:, :Ta])
    y1c = F.center_rows(panel_a.Y[-1, :Ta][None, :])[0]
    n_full = Yc.shape[0]
    out = []
    for rank, size in enumerate(DEG_SIZES):
        draws = []
        for j in range(DEG_DRAWS):
            if size >= n_full:
                sel = np.arange(n_full)
                idx = 300 + rank * DEG_DRAWS + j
            else:
                rng = np.random.default_rng(61301 + 100 * rank + j)
                sel = rng.choice(n_full, size=size, replace=False)
                idx = 300 + rank * DEG_DRAWS + j
            pool = np.random.default_rng(
                60101 + idx * 10**6).standard_normal((F.G_NULL, Ta))
            res = F.analyze_space(Yc[sel], y1c, pool, SEEDS, idx,
                                  with_bootstrap=False)
            draws.append({"k": res["k"], "d": res["d"],
                          "p_align": res["p_align"]})
        ks = np.array([dr["k"] for dr in draws])
        ds = np.array([dr["d"] for dr in draws])
        out.append({"size": size, "median_d": float(np.median(ds)),
                    "share_silent": float((ks == 0).mean()),
                    "share_misaligned": float(np.mean([
                        dr["k"] >= 1 and dr["p_align"] >= F.ALIGN_P
                        for dr in draws]))})
        print(f"  n_d={size}: median d={out[-1]['median_d']:.3g} "
              f"silent={out[-1]['share_silent']:.2f}", flush=True)
    sizes_ge4 = [o for o in out if o["size"] >= 4]
    sensitive = any(o["share_silent"] >= 0.5 for o in sizes_ge4)
    full_loud = out[0]["share_silent"] < 0.5
    summary = {"curve": out, "full_panel_loud": bool(full_loud),
               "transition_by_size_ge4": bool(sensitive)}
    with open(OUT / "wave2_degradation.json", "w") as fh:
        json.dump(_clean(summary), fh, indent=2)
    return summary


def write_figure(basque_out, basque_rows, deg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    e = np.asarray(basque_out["treated"]["evals_top12"])
    axes[0].semilogy(range(1, len(e) + 1), e, "ko-", ms=4)
    axes[0].axhline(basque_out["treated"]["edge"], color="crimson",
                    ls="--", label=f"MP edge ({basque_out['treated']['edge']:.3g})")
    axes[0].set_title(
        f"Basque donor spectrum: k={basque_out['treated']['k']}, "
        f"d={basque_out['treated']['d']:.2f}, "
        f"{basque_out['treated_label']}", fontsize=9)
    axes[0].legend(fontsize=8)

    ds = np.array([max(r["d"], 1e-3) for r in basque_rows])
    rs = np.array([r["rmse"] for r in basque_rows])
    frag = np.array([r["label"].startswith("FRAGILE")
                     or r["label"].startswith("INCONCLUSIVE")
                     for r in basque_rows])
    axes[1].scatter(ds[frag], rs[frag], s=26, c="tab:orange",
                    label="FRAGILE/INCONCLUSIVE")
    axes[1].scatter(ds[~frag], rs[~frag], s=26, c="tab:green",
                    label="RECOVERABLE")
    axes[1].axvline(max(basque_out["treated_d"], 1e-3), color="k", lw=2,
                    label=f"treated d={basque_out['treated_d']:.2f}")
    axes[1].set_xscale("symlog", linthresh=0.01)
    axes[1].set_xlabel("d (BBP-edge units)")
    axes[1].set_ylabel("outcome-only SCM pre-fit RMSE")
    axes[1].set_title("Basque placebos", fontsize=10)
    axes[1].legend(fontsize=8)

    sizes = [o["size"] for o in deg["curve"]]
    med = [max(o["median_d"], 1e-3) for o in deg["curve"]]
    sil = [o["share_silent"] for o in deg["curve"]]
    axes[2].plot(sizes, med, "ko-")
    axes[2].set_xscale("log")
    axes[2].set_xlabel("donor-pool size n_d")
    axes[2].set_ylabel("median d-hat")
    ax2 = axes[2].twinx()
    ax2.plot(sizes, sil, "b^--", alpha=0.7)
    ax2.set_ylabel("share silent (k=0)", color="b")
    axes[2].set_title("Degradation battery (smoking)", fontsize=10)
    fig.suptitle("WP-D2b Wave-2 (preregistration_d2b_addendum)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_d2_wave2.png", dpi=150)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    basque_out, basque_rows, _ = run_basque()
    deg = run_degradation()
    write_figure(basque_out, basque_rows, deg)
    print("=" * 76)
    nf = basque_out["nf"]
    print(f"BASQUE: label={basque_out['treated_label']} "
          f"d={basque_out['treated_d']:.3f} "
          f"CI=[{basque_out['treated_boot_donor']['q05']:.2f},"
          f"{basque_out['treated_boot_donor']['q95']:.2f}] "
          f"p_align={basque_out['treated_p_align']:.4f} "
          f"preRMSE={basque_out['treated_pre_rmse']:.3f}")
    print(f"  N1={nf['N1_fragile_treated']} N2={nf['N2_certified_safe']} "
          f"N3={nf['N3_certification_asymmetry']} "
          f"P1={nf['P1_fit_ordering']['pass']} "
          f"(rho={nf['P1_fit_ordering']['rho']:.2f}) "
          f"P2={nf['P2_treated_separation']} P3={nf['P3_date_alarms_pass']} "
          f"stab={nf['stability_pass']} -> W1_PASS={nf['W1_PASS']}")
    print(f"DEGRADATION: full loud={deg['full_panel_loud']} "
          f"transition by size>=4: {deg['transition_by_size_ge4']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
