"""WP-D2 decisive driver: distance-to-frontier on the two canonical panels.

Implements preregistration_d2_addendum.md exactly (frozen 2026-08-26
BEFORE this run). Outputs land in results_d2/ and figures/; this script is
the ONLY writer of those artifacts. Exit code is 0 iff the script ran to
completion; gate outcomes are recorded in results_d2/summary.json and
judged in gate_g4_memo.md, never tuned here.
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
from scm_frontier.estimators import mc_nn_cv, spectral_sc  # noqa: E402

SEEDS = {"align": 60_101, "donor": 60_102, "time": 60_103, "perm": 60_104}
PLACEBO_V_SEED = 60_001
Z_SHIFT_SEED = 8_880_001
POSTTEST_SEED = 7_770_001
CUTOFF_YEARS = {
    "smoking_prop99": [1978, 1980, 1982],
    "german_reunification": [1972, 1976, 1980, 1984],
}
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


def year_indices(panel, years):
    return [int(np.where(panel.years == y)[0][0]) for y in years]


def treated_analysis(name, Yd_pre, y1_pre, Yd_post, cutoff_idx):
    print(f"[{name}] treated-unit spectral pipeline ...", flush=True)
    res = F.analyze_unit(Yd_pre, y1_pre, SEEDS, idx=0)
    res["date_alarms"] = F.date_alarm_scan(
        Yd_pre, cutoff_idx, Z_SHIFT_SEED, res["sigma2_1"])
    res["alarm_rate"] = float(np.mean([a["alarm"] for a in res["date_alarms"]]))
    res["P3_pass"] = bool(res["alarm_rate"] <= F.ALARM_MAX)
    res["posttest_control"] = F.posttest_control(
        Yd_pre, Yd_post, res["sigma2_1"], POSTTEST_SEED)
    return res


def smoking_placebos(panel, cdf):
    T0 = panel.T0
    X = L.smoking_adh_predictors(panel, cdf)
    order = list(panel.units)
    donors_all = [u for u in order if u != L.SMOKING_TREATED]
    rows = []
    for i, name in enumerate(donors_all):
        keep = [u for u in order if u not in (name, L.SMOKING_TREATED)]
        idx_map = [order.index(u) for u in keep + [name]]
        Yn = panel.Y[idx_map]
        fit = S.synth_adh_smoking(X[idx_map], Yn, T0, seed=PLACEBO_V_SEED)
        path = fit["w"] @ Yn[:-1]
        rmse = float(np.sqrt(np.mean((Yn[-1, :T0] - path[:T0]) ** 2)))
        diag = F.analyze_unit(Yn[:-1, :T0], Yn[-1, :T0], SEEDS, idx=i + 1)
        rows.append({"unit": name, "rmse": rmse,
                     **{k: diag[k] for k in
                        ["k", "d", "p_align", "r_align", "label",
                         "boot_donor", "stability_pass"]}})
        print(f"  placebo {i + 1}/38 {name}: d={diag['d']:.3f} "
              f"k={diag['k']} p={diag['p_align']:.3f} "
              f"rmse={rmse:.2f} {diag['label']}", flush=True)
    return rows


def germany_placebos(panel, raw):
    T0 = panel.T0
    order = list(panel.units)
    Xt_all = L.germany_predictors(raw, "training")
    Xm_all = L.germany_predictors(raw, "main")
    donors_all = [u for u in order if u != L.GERMANY_TREATED]
    rows = []
    for i, name in enumerate(donors_all):
        keep = [u for u in order if u not in (name, L.GERMANY_TREATED)]
        idx_map = [order.index(u) for u in keep + [name]]
        Yn = panel.Y[idx_map]
        Z0 = Yn[:-1][:, T0 - 9:T0 + 1]
        Z1 = Yn[-1, T0 - 9:T0 + 1]
        fit = S.synth_adh_germany(Xt_all[idx_map], len(idx_map) - 1,
                                  Z0, Z1, Xm_all[idx_map],
                                  seed=PLACEBO_V_SEED)
        path = fit["w"] @ Yn[:-1]
        lvl = float(Yn[-1, :T0].mean())
        rmse = float(100.0 * np.sqrt(
            np.mean((Yn[-1, :T0] - path[:T0]) ** 2)) / lvl)
        diag = F.analyze_unit(Yn[:-1, :T0], Yn[-1, :T0], SEEDS, idx=i + 1)
        rows.append({"unit": name, "rmse": rmse,
                     **{k: diag[k] for k in
                        ["k", "d", "p_align", "r_align", "label",
                         "boot_donor", "stability_pass"]}})
        print(f"  placebo {i + 1}/16 {name}: d={diag['d']:.3f} "
              f"k={diag['k']} p={diag['p_align']:.3f} "
              f"rmse={rmse:.3f}%lvl {diag['label']}", flush=True)
    return rows


def incumbents_smoking(panel, cdf, spec):
    T0 = panel.T0
    X = L.smoking_adh_predictors(panel, cdf)
    fit = S.synth_adh_smoking(X, panel.Y, T0, seed=PLACEBO_V_SEED)
    path = fit["w"] @ panel.Y[:-1]
    scm_pre_rmse = float(np.sqrt(np.mean(
        (panel.Y[-1, :T0] - path[:T0]) ** 2)))
    pred_mc = mc_nn_cv(panel.Y[:-1], panel.Y[-1, :T0])
    att_mc = float((panel.Y[-1, T0:] - pred_mc).mean())
    pred_sp = spectral_sc(panel.Y[:-1], panel.Y[-1, :T0],
                          sigma=spec["sigma2_1"] ** 0.5, c=spec["c"])
    att_sp = float((panel.Y[-1, T0:] - pred_sp).mean())
    return {"scm_classic_pre_rmse": scm_pre_rmse, "mc_nn_att": att_mc,
            "spectral_gated_att": att_sp}


def incumbents_germany(panel, raw, spec):
    T0 = panel.T0
    Xt = L.germany_predictors(raw, "training")
    Xm = L.germany_predictors(raw, "main")
    mask = np.arange(len(panel.units)) != panel.treated_idx
    Z0 = panel.Y[mask][:, T0 - 9:T0 + 1]
    Z1 = panel.Y[panel.treated_idx, T0 - 9:T0 + 1]
    fit = S.synth_adh_germany(Xt, panel.treated_idx, Z0, Z1, Xm,
                              seed=PLACEBO_V_SEED)
    path = fit["w"] @ panel.Y[mask]
    lvl = float(panel.Y[panel.treated_idx, :T0].mean())
    scm_pre_rmse_pct = float(100.0 * np.sqrt(np.mean(
        (panel.Y[panel.treated_idx, :T0] - path[:T0]) ** 2)) / lvl)
    pred_mc = mc_nn_cv(panel.Y[mask], panel.Y[panel.treated_idx, :T0])
    att_mc = float((panel.Y[panel.treated_idx, T0:] - pred_mc).mean())
    pred_sp = spectral_sc(panel.Y[mask], panel.Y[panel.treated_idx, :T0],
                          sigma=spec["sigma2_1"] ** 0.5, c=spec["c"])
    att_sp = float((panel.Y[panel.treated_idx, T0:] - pred_sp).mean())
    return {"scm_classic_pre_rmse_pct_lvl": scm_pre_rmse_pct,
            "mc_nn_att": att_mc, "spectral_gated_att": att_sp}


def assemble_panel(treated, placebos, treated_rmse, poor_looking, name):
    treated["poor_looking"] = bool(poor_looking)
    d_vec = np.array([r["d"] for r in placebos])
    rmse_vec = np.array([r["rmse"] for r in placebos])
    labels = [r["label"] for r in placebos]
    rho, perm_p = F.spearman_with_perm(d_vec, rmse_vec, SEEDS["perm"])
    thr = F.RHO_THRESH[name]
    p1 = {"rho": rho, "perm_p": perm_p, "threshold": thr,
          "pass": bool(rho <= thr and perm_p <= 0.05)}
    nf = F.nf_arms(treated, treated_rmse, [r["unit"] for r in placebos],
                   d_vec, rmse_vec, labels)
    nf["P1_fit_ordering"] = p1
    nf["P3_date_alarms_pass"] = treated["P3_pass"]
    nf["stability_pass"] = treated["stability_pass"]
    nf["placebo_label_counts"] = {
        lab: int(sum(1 for l in labels if l == lab))
        for lab in sorted(set(labels))}
    return nf


def write_figure(t_a, t_b, pl_a, pl_b, nf_a, nf_b, rmse_a, rmse_b):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))
    for ax, res, ttl in [(axes[0, 0], t_a, "Smoking donor spectrum"),
                         (axes[0, 1], t_b, "Germany donor spectrum")]:
        e = np.asarray(res["evals_top12"])
        ax.semilogy(range(1, len(e) + 1), e, "ko-", ms=4)
        ax.axhline(res["edge"], color="crimson", ls="--",
                   label=f"MP edge ({res['edge']:.3g})")
        ax.set_xlabel("eigenvalue rank")
        ax.set_title(f"{ttl}: k={res['k']}, d={res['d']:.2f}, "
                     f"{res['label']}", fontsize=10)
        ax.legend(fontsize=8)
    for ax, pl, nf, rt, ttl in [
            (axes[1, 0], pl_a, nf_a, rmse_a,
             "Smoking placebos: d vs pre-fit RMSE"),
            (axes[1, 1], pl_b, nf_b, rmse_b,
             "Germany placebos: d vs pre-fit RMSE (% of level)")]:
        ds = np.array([max(r["d"], 1e-3) for r in pl])
        rs = np.array([r["rmse"] for r in pl])
        frag = [r["label"].startswith("FRAGILE") or r["label"].startswith(
            "INCONCLUSIVE") for r in pl]
        ax.scatter(ds[frag], rs[frag], s=24, c="tab:orange",
                   label="FRAGILE/INCONCLUSIVE")
        ax.scatter(ds[~np.array(frag)], rs[~np.array(frag)], s=24,
                   c="tab:green", label="RECOVERABLE")
        ax.axvline(max(rt, 1e-3), color="k", lw=2,
                   label=f"treated d={rt:.2f}")
        ax.axvline(nf["P2_q90_placebo_d"], color="gray", ls=":",
                   label=f"Q90(placebo d)={nf['P2_q90_placebo_d']:.2f}")
        ax.set_xscale("symlog", linthresh=0.01)
        ax.set_xlabel("d (BBP-edge units)")
        ax.set_ylabel("classic-SCM pre-fit RMSE")
        ax.set_title(ttl, fontsize=10)
        ax.legend(fontsize=8)
    fig.suptitle("WP-D2 distance-to-frontier "
                 "(preregistration_d2_addendum, frozen 2026-08-26)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_d2_distance_to_frontier.png", dpi=150)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)

    print("=== Panel A: California smoking ===", flush=True)
    panel_a, cdf_a = L.load_smoking()
    Ta = panel_a.T0
    don_a = np.arange(len(panel_a.units)) != panel_a.treated_idx
    t_a = treated_analysis("smoking_prop99", panel_a.Y[don_a][:, :Ta],
                           panel_a.Y[-1, :Ta], panel_a.Y[don_a][:, Ta:],
                           year_indices(panel_a, CUTOFF_YEARS["smoking_prop99"]))
    inc_a = incumbents_smoking(panel_a, cdf_a, t_a)
    print("  placebos ...", flush=True)
    pl_a = smoking_placebos(panel_a, cdf_a)

    print("=== Panel B: German reunification ===", flush=True)
    panel_b = L.load_germany()
    raw_b = L.germany_raw()
    Tb = panel_b.T0
    don_b = np.arange(len(panel_b.units)) != panel_b.treated_idx
    t_b = treated_analysis("german_reunification", panel_b.Y[don_b][:, :Tb],
                           panel_b.Y[panel_b.treated_idx, :Tb],
                           panel_b.Y[don_b][:, Tb:],
                           year_indices(panel_b,
                                        CUTOFF_YEARS["german_reunification"]))
    inc_b = incumbents_germany(panel_b, raw_b, t_b)
    print("  placebos ...", flush=True)
    pl_b = germany_placebos(panel_b, raw_b)

    rmse_ta = inc_a["scm_classic_pre_rmse"]
    rmse_tb = inc_b["scm_classic_pre_rmse_pct_lvl"]
    q75a = float(np.percentile([r["rmse"] for r in pl_a], 75))
    q75b = float(np.percentile([r["rmse"] for r in pl_b], 75))
    nf_a = assemble_panel(t_a, pl_a, rmse_ta, rmse_ta > q75a,
                          "smoking_prop99")
    nf_b = assemble_panel(t_b, pl_b, rmse_tb, rmse_tb > q75b,
                          "german_reunification")

    ident_ok = all(t["label"] != "FRAGILE-INVISIBLE" for t in (t_a, t_b))
    applied_any = any(
        nf["NF"] and nf["P2_treated_separation"]
        and nf["P3_date_alarms_pass"] and nf["stability_pass"]
        for nf in (nf_a, nf_b))
    verdict_full = applied_any and any(
        nf["NF"] and nf["P1_fit_ordering"]["pass"]
        for nf in (nf_a, nf_b) if
        (nf["NF"] and nf["P2_treated_separation"]
         and nf["P3_date_alarms_pass"] and nf["stability_pass"]))

    out_summary = {
        "identification_control_pass": ident_ok,
        "applied_value_pass_any_panel": bool(applied_any),
        "verdict_full_pass": bool(verdict_full),
        "panels": {
            "smoking_prop99": {
                "treated_label": t_a["label"],
                "treated_d": t_a["d"],
                "treated_boot_donor": t_a["boot_donor"],
                "nf": nf_a, "incumbents": inc_a,
                "treated_pre_rmse": rmse_ta,
                "placebo_rmse_q75": q75a},
            "german_reunification": {
                "treated_label": t_b["label"],
                "treated_d": t_b["d"],
                "treated_boot_donor": t_b["boot_donor"],
                "nf": nf_b, "incumbents": inc_b,
                "treated_pre_rmse": rmse_tb,
                "placebo_rmse_q75": q75b},
        },
    }

    for tag, res in [("smoking", t_a), ("germany", t_b)]:
        with open(OUT / f"distance_{tag}.json", "w") as fh:
            json.dump(_clean(res), fh, indent=2)
    for tag, rows in [("smoking", pl_a), ("germany", pl_b)]:
        keys = ["unit", "rmse", "k", "d", "p_align", "r_align", "label",
                "stability_pass", "boot_donor_q05", "boot_donor_q95"]
        with open(OUT / f"placebos_{tag}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(keys)
            for r in rows:
                w.writerow([
                    r["unit"], r["rmse"], r["k"], r["d"], r["p_align"],
                    r["r_align"], r["label"], r["stability_pass"],
                    r["boot_donor"]["q05"], r["boot_donor"]["q95"]])
    with open(OUT / "summary.json", "w") as fh:
        json.dump(_clean(out_summary), fh, indent=2)

    write_figure(t_a, t_b, pl_a, pl_b, nf_a, nf_b, rmse_ta, rmse_tb)

    print("=" * 76)
    print(f"A: label={t_a['label']} d={t_a['d']:.3f} "
          f"CI=[{t_a['boot_donor']['q05']:.2f},"
          f"{t_a['boot_donor']['q95']:.2f}] "
          f"p_align={t_a['p_align']:.4f} NF={nf_a['NF']} "
          f"P1={nf_a['P1_fit_ordering']['pass']} P2="
          f"{nf_a['P2_treated_separation']} P3="
          f"{nf_a['P3_date_alarms_pass']}")
    print(f"B: label={t_b['label']} d={t_b['d']:.3f} "
          f"CI=[{t_b['boot_donor']['q05']:.2f},"
          f"{t_b['boot_donor']['q95']:.2f}] "
          f"p_align={t_b['p_align']:.4f} NF={nf_b['NF']} "
          f"P1={nf_b['P1_fit_ordering']['pass']} P2="
          f"{nf_b['P2_treated_separation']} P3="
          f"{nf_b['P3_date_alarms_pass']}")
    print(f"identification_control_pass={ident_ok} "
          f"applied_value_pass_any_panel={out_summary['applied_value_pass_any_panel']} "
          f"verdict_full_pass={out_summary['verdict_full_pass']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
