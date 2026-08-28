"""Phase E decisive driver: Bonander Florida SYG (prereg_bonander_addendum.md).

Frozen 2026-08-28 BEFORE decisive run (except the single exploratory probe that
forced the BBP lam/sigma2 fix). This is the ONLY writer of results_e/bonander/.
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

SEEDS = {"align": 70_001, "donor": 70_002, "time": 70_003, "perm": 70_004}
Z_SHIFT_SEED = 8_880_001
POSTTEST_SEED = 7_770_001
# pseudo-cutoffs inside pre 1:81  -> 27≈2001-03, 54≈2003-06
CUTOFF_IDX = [27, 54]
OUT = ROOT / "results_e" / "bonander"
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


def simple_rmse(y1: np.ndarray, Yd: np.ndarray) -> float:
    # OLS w = argmin ||y1 - w'Yd||_2  (no simplex, frozen proxy for pre-fit)
    # Yd n x T0, y1 T0
    w, *_ = np.linalg.lstsq(Yd.T, y1, rcond=None)
    pred = w @ Yd
    return float(np.sqrt(np.mean((y1 - pred) ** 2)))


def treated_analysis(Yd_pre, y1_pre, Yd_post, cutoff_idx):
    print("[bonander] treated spectral pipeline ...", flush=True)
    res = F.analyze_unit(Yd_pre, y1_pre, SEEDS, idx=0)
    res["date_alarms"] = F.date_alarm_scan(Yd_pre, cutoff_idx, Z_SHIFT_SEED, res["sigma2_1"])
    res["alarm_rate"] = float(np.mean([a["alarm"] for a in res["date_alarms"]]))
    res["P3_pass"] = bool(res["alarm_rate"] <= F.ALARM_MAX)
    res["posttest_control"] = F.posttest_control(Yd_pre, Yd_post, res["sigma2_1"], POSTTEST_SEED)
    return res


def bonander_placebos(panel):
    T0 = int(panel.meta["T0"])
    order = list(panel.units)
    treated = L.BONANDER_TREATED
    rows = []
    for i, name in enumerate(L.BONANDER_DONORS):
        # donor pool for placebo i = all 16 except name
        keep = [u for u in order if u != name]
        idx = [order.index(u) for u in keep + [name]]
        Yn = panel.Y[idx]  # (16) x 192, last is placebo
        Yd_pre = Yn[:-1, :T0]
        y1_pre = Yn[-1, :T0]
        rmse = simple_rmse(y1_pre, Yd_pre)
        diag = F.analyze_unit(Yd_pre, y1_pre, SEEDS, idx=i + 1)
        rows.append({"unit": name, "rmse": rmse,
                     **{k: diag[k] for k in ["k", "d", "p_align", "r_align", "label", "boot_donor", "stability_pass", "evals_top12"]}})
        print(f"  placebo {i+1}/15 {name}: d={diag['d']:.2f} k={diag['k']} p={diag['p_align']:.3f} rmse={rmse:.3f} {diag['label']}", flush=True)
    return rows


def assemble_panel(treated, placebos, treated_rmse, poor_looking):
    treated["poor_looking"] = bool(poor_looking)
    d_vec = np.array([r["d"] for r in placebos])
    rmse_vec = np.array([r["rmse"] for r in placebos])
    labels = [r["label"] for r in placebos]
    rho, perm_p = F.spearman_with_perm(d_vec, rmse_vec, SEEDS["perm"])
    thr = -0.30
    p1 = {"rho": rho, "perm_p": perm_p, "threshold": thr, "pass": bool(rho <= thr and perm_p <= 0.05)}
    nf = F.nf_arms(treated, treated_rmse, [r["unit"] for r in placebos], d_vec, rmse_vec, labels)
    nf["P1_fit_ordering"] = p1
    nf["P3_date_alarms_pass"] = treated["P3_pass"]
    nf["stability_pass"] = treated["stability_pass"]
    nf["placebo_label_counts"] = {lab: int(sum(1 for l in labels if l == lab)) for lab in sorted(set(labels))}
    return nf


def write_figure(t_res, placebos, nf, treated_rmse):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    e = np.asarray(t_res["evals_top12"])
    ax = axes[0]
    ax.semilogy(range(1, len(e) + 1), e, "ko-", ms=4)
    ax.axhline(t_res["edge"], color="crimson", ls="--", label=f"MP edge {t_res['edge']:.3g}")
    ax.set_xlabel("eigenvalue rank")
    ax.set_title(f"Bonander HomicideRates donor spectrum: k={t_res['k']}, d={t_res['d']:.2f}, {t_res['label']}", fontsize=10)
    ax.legend(fontsize=8)

    ax = axes[1]
    ds = np.array([max(r["d"], 1e-3) for r in placebos])
    rs = np.array([r["rmse"] for r in placebos])
    frag = np.array([r["label"].startswith("FRAGILE") or r["label"].startswith("INCONCLUSIVE") for r in placebos])
    ax.scatter(ds[frag], rs[frag], s=36, c="tab:orange", label="FRAGILE/INCONCLUSIVE")
    ax.scatter(ds[~frag], rs[~frag], s=36, c="tab:green", label="RECOVERABLE")
    ax.axvline(max(treated_rmse, 1e-3) if False else max(t_res["d"], 1e-3), color="k", lw=2, label=f"treated d={t_res['d']:.2f}")
    # actual treated d line
    ax.axvline(max(t_res["d"], 1e-3), color="k", lw=2, label=f"treated d={t_res['d']:.2f}")
    ax.axvline(nf["P2_q90_placebo_d"] if np.isfinite(nf["P2_q90_placebo_d"]) else 0, color="gray", ls=":", label=f"Q90 d={nf['P2_q90_placebo_d']:.2f}")
    ax.set_xscale("symlog", linthresh=0.01)
    ax.set_xlabel("d (BBP-edge units)")
    ax.set_ylabel("OLS pre-fit RMSE (HomicideRates p100k)")
    ax.set_title("Bonander placebos: d vs pre-fit RMSE", fontsize=10)
    ax.legend(fontsize=8)
    fig.suptitle("Bonander Florida SYG — distance-to-frontier (prereg_bonander_addendum 2026-08-28)")
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "fig_bonander_distance_to_frontier.png", dpi=150)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("=== Bonander Florida SYG (primary HomicideRates) ===", flush=True)
    panel = L.load_bonander(outcome="HomicideRates")
    T0 = int(panel.meta["T0"])
    treat_time = int(panel.meta["treat_time"])
    don = np.arange(len(panel.units)) != panel.treated_idx
    # split at T0=81 (time 1:81 pre, 82:192 post)
    Yd_pre = panel.Y[don][:, :T0]
    y1_pre = panel.Y[panel.treated_idx, :T0]
    Yd_post = panel.Y[don][:, T0:]  # actually time 82:192 is post, but Y columns 0-indexed: :81 pre, 81: post
    # correction: panel.Y columns are time 1..192 indexed 0..191, so T0=81 means pre 0:81, post 81:192
    # Already done: Yd_pre = :T0 (0:81 = time1:81), Yd_post = T0: (81:192 = time82:192)
    t_res = treated_analysis(Yd_pre, y1_pre, Yd_post, CUTOFF_IDX)
    print(f"  treated: k={t_res['k']} d={t_res['d']:.2f} p={t_res['p_align']:.4f} label={t_res['label']} q05={t_res['boot_donor']['q05']:.2f}", flush=True)
    # OLS pre-fit for treated for NF poor-looking
    rmse_t = simple_rmse(y1_pre, Yd_pre)
    print(f"  treated OLS rmse={rmse_t:.4f}", flush=True)
    print("  placebos ...", flush=True)
    pl = bonander_placebos(panel)
    q75 = float(np.percentile([r["rmse"] for r in pl], 75))
    nf = assemble_panel(t_res, pl, rmse_t, rmse_t > q75)
    # secondaries descriptive
    secondaries = {}
    for outc in ["p100khomicide_firearm", "p100khomicide_exclfirearm", "Firearm.Suicide.Rates", "Homicide_count"]:
        pan2 = L.load_bonander(outcome=outc)
        Yd2 = pan2.Y[don][:, :T0]
        y1_2 = pan2.Y[pan2.treated_idx, :T0]
        r2 = F.analyze_unit(Yd2, y1_2, SEEDS, idx=100 + hash(outc) % 1000)
        secondaries[outc] = {k: r2[k] for k in ["k", "d", "p_align", "r_align", "label", "sigma2_1"]}
        print(f"  secondary {outc}: k={r2['k']} d={r2['d']:.2f} p={r2['p_align']:.3f} {r2['label']}", flush=True)

    ident_ok = t_res["label"] != "FRAGILE-INVISIBLE"
    applied_pass = bool(nf["NF"] and nf["P2_treated_separation"] and nf["P3_date_alarms_pass"] and nf["stability_pass"])
    verdict_full = bool(applied_pass and nf["P1_fit_ordering"]["pass"])

    summary = {
        "panel": "bonander_florida_syg",
        "treat_time": treat_time, "T0": T0, "n_d": int(np.sum(don)),
        "treated": {"label": t_res["label"], "d": t_res["d"], "boot_donor": t_res["boot_donor"], "p_align": t_res["p_align"], "rmse": rmse_t, "q75_placebo_rmse": q75},
        "nf": nf, "identification_control_pass": ident_ok,
        "applied_value_pass": applied_pass, "verdict_full_pass": verdict_full,
        "secondaries": secondaries,
    }

    with open(OUT / "distance_bonander.json", "w") as fh:
        json.dump(_clean(t_res), fh, indent=2)
    with open(OUT / "placebos_bonander.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["unit", "rmse", "k", "d", "p_align", "r_align", "label", "stability_pass", "boot_donor_q05", "boot_donor_q95"])
        for r in pl:
            w.writerow([r["unit"], r["rmse"], r["k"], r["d"], r["p_align"], r["r_align"], r["label"], r["stability_pass"], r["boot_donor"]["q05"], r["boot_donor"]["q95"]])
    with open(OUT / "summary_bonander.json", "w") as fh:
        json.dump(_clean(summary), fh, indent=2)
    # also write raw panel copy for provenance? already in data/raw
    write_figure(t_res, pl, nf, rmse_t)

    print("=" * 76)
    print(f"Florida label={t_res['label']} d={t_res['d']:.2f} CI=[{t_res['boot_donor']['q05']:.2f},{t_res['boot_donor']['q95']:.2f}] p={t_res['p_align']:.4f} NF={nf['NF']} P1={nf['P1_fit_ordering']['pass']} P2={nf['P2_treated_separation']} P3={nf['P3_date_alarms_pass']} stability={t_res['stability_pass']}")
    print(f"ident={ident_ok} applied_pass={applied_pass} full={verdict_full}")
    print(f"placebo labels {nf['placebo_label_counts']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
