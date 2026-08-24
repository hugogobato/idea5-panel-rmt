"""Phase C analysis + figures (runs locally AFTER merge_shards.py passes).

Produces in figures/:
  fig_threshold_kink.png    gated spectral RMSE vs m per c column + F overlay
  fig_de_overlay.png        DE ansatz F vs simulation (full arm), all methods
  fig_size_calibration.png  empirical size/power of Z_tw / Z_boot / t-test
  fig_baseline_favorable.png C2(ii) dense-weak battery
  fig_rank_accuracy.png     rank-selector head-to-head on C1 cells
Plus analysis_outputs.json with every preregistered number for gate_g3_memo.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


def frontier_rmse(m: float, c: float, T0: int, theta: float) -> float:
    """Conjectured DE risk F (frontier_ansatz.md Section 3), r=1, K = {s>sqrt(c)}."""
    s = m * np.sqrt(c)
    if s <= np.sqrt(c):
        return float(np.sqrt(1.0 + theta))
    lam = 1.0 + s + c + c / s
    zeta = (1.0 - c / s**2) / (1.0 + c / s)
    tau = np.sqrt(s / lam)
    excess = theta * ((1.0 - zeta * tau) ** 2 + zeta / lam) + (s * zeta + 1.0) / (T0 * lam)
    return float(np.sqrt(1.0 + excess))


def kink_location(ms: np.ndarray, rmse: np.ndarray) -> tuple[float, np.ndarray]:
    """Preregistered estimator: argmax discrete curvature within [0.6, 1.6]."""
    mask = (ms >= 0.6 - 1e-9) & (ms <= 1.6 + 1e-9)
    idx = np.where(mask)[0]
    curv = np.full(len(ms), np.nan)
    for i in idx[1:-1]:
        curv[i] = rmse[i + 1] - 2 * rmse[i] + rmse[i - 1]
    best = idx[np.nanargmax(curv[idx[1:-1]]) + 1]
    return float(ms[best]), curv


def main() -> None:
    out: dict[str, object] = {}
    p_c1 = ROOT / "results_c1" / "risk_curves.parquet"
    if not p_c1.exists():
        print(f"missing {p_c1}; run scripts/merge_shards.py after shards arrive")
        return
    df = pd.read_parquet(p_c1)

    # ---------------- WP-C1 kink criterion ---------------------------------
    spec = df[(df.experiment == "c1") & (df.method == "spectral_gated")
              & (df.arm == "full") & (df.r == 1)]
    kink_table = []
    for c, sub in spec.groupby("c"):
        g = sub.groupby("m")["rmse"].agg(["mean", "sem"]).reset_index()
        ms = g["m"].to_numpy(float)
        mk, _ = kink_location(ms, g["mean"].to_numpy(float))
        ok = abs(mk - 1.0) <= 0.15
        kink_table.append(dict(c=float(c), m_kink=mk, within_15pct=bool(ok),
                               n_reps=int(len(sub) / len(g))))
        print(kink_table[-1])
    share = np.mean([k["within_15pct"] for k in kink_table])
    out["kink"] = dict(table=kink_table,
                       share_within_15pct=float(share),
                       criterion_pass=bool(share >= 0.8))
    out["onset"] = json.loads((ROOT / "results_raw" / "onset" /
                               "onset_verdict.json").read_text()) \
        if (ROOT / "results_raw" / "onset" / "onset_verdict.json").exists() \
        else None

    # ---------------- figures: threshold kink + DE overlay ------------------
    fig, axes = plt.subplots(1, len(kink_table), figsize=(4.2 * len(kink_table), 3.8),
                             squeeze=False)
    for ax, row in zip(axes[0], kink_table):
        c = row["c"]
        T0 = int(spec[spec.c == c]["T0"].iloc[0])
        sub = spec[spec.c == c].groupby("m")["rmse"]
        g = sub.agg(["mean", "sem"]).reset_index()
        ax.errorbar(g.m, g["mean"], yerr=1.96 * g["sem"], fmt="o-", ms=3.5,
                    lw=1, color="#1f77b4", label="spectral SC (gated)")
        mm = np.linspace(g.m.min(), g.m.max(), 300)
        ax.plot(mm, [frontier_rmse(v, c, T0, 1.0) for v in mm], "--",
                color="crimson", label="ansatz $F$ (theta=1)")
        for meth, col in (("scm_simplex", "#7f7f7f"), ("mc_nn_cv", "#2ca02c"),
                          ("ridge_sc", "#9467bd"), ("donor_mean", "#c7c7c7")):
            gm = df[(df.experiment == "c1") & (df.method == meth)
                    & (df.arm == "full") & (df.r == 1)
                    & (df.c == c)].groupby("m")["rmse"].mean()
            ax.plot(gm.index, gm.values, ":", color=col, lw=1.1, alpha=0.85,
                    label=meth)
        ax.axvline(1.0, color="k", lw=0.7, alpha=0.5)
        ax.axvline(row["m_kink"], color="orange", lw=0.8, ls="-.",
                   label=f"kink est {row['m_kink']:.2f}")
        ax.set_title(f"c = {c} ({'PASS' if row['within_15pct'] else 'FAIL'})")
        ax.set_xlabel("spike multiplier m"); ax.set_ylabel("RMSE / sigma")
        ax.legend(fontsize=6)
    fig.suptitle("WP-C1: risk vs spike strength, full-alignment arm, r = 1")
    fig.tight_layout(); fig.savefig(FIG / "fig_de_overlay.png", dpi=150); plt.close(fig)

    # ---------------- bite criterion ---------------------------------------
    incs = df[(df.experiment == "c1") & (df.method.isin(
        ["scm_simplex", "ridge_sc", "mc_nn_cv"])) & (df.c >= 1)]
    worst = incs.groupby(["method", "arm", "r"])["rmse"].mean().max()
    flag_sub = df[(df.experiment == "c1") & (df.method == "_diag")
                  & (((df.m <= 0.9) & (df.theta > 0)) | (df.arm == "orthogonal"))]
    flag_ok = df[(df.experiment == "c1") & (df.method == "_diag")
                 & (df.m >= 1.5) & (df.arm != "orthogonal")]
    power = float((flag_sub["k_selected"].astype(int) == 0).mean())
    false_alarm = float((flag_ok["k_selected"].astype(int) == 0).mean())
    out["bite"] = dict(worst_incumbent_rmse=float(worst),
                       diagnostic_power_subfrontier=power,
                       false_alarm_rate_recoverable=false_alarm,
                       bite_pass=bool(worst >= 2.0 and power >= 0.8
                                      and false_alarm <= 0.2))

    # ---------------- rank accuracy head-to-head ---------------------------
    dg = df[df.method == "_diag"].copy()
    dg["k_true"] = dg["r"]
    sup = dg[(dg.m >= 1.2) & (dg.arm != "orthogonal")]
    acc = pd.DataFrame({
        "gated": (sup.k_selected.astype(int) == sup.k_true).groupby(
            sup.c).mean(),
        "cv_rank": (sup.cv_rank.astype(int) == sup.k_true).groupby(
            sup.c).mean(),
        "ungated": (sup.k_selected.astype(int) == sup.k_true).groupby(
            sup.c).mean(),
    })
    acc.plot.bar(figsize=(7, 3.8))
    plt.ylabel("rank accuracy (m >= 1.2, aligned)")
    plt.tight_layout(); plt.savefig(FIG / "fig_rank_accuracy.png", dpi=150)
    plt.close(fig)
    mid = acc.loc[[c for c in (0.5, 1.0, 2.0) if c in acc.index]] \
        if len(acc) else acc
    out["rank_head_to_head"] = dict(
        accuracy_by_c=acc.reset_index().to_dict("records"),
        c3_pass=bool(len(mid) and ((mid.gated - mid.cv_rank) >= 0.10).all()))

    # ---------------- C2 batteries ----------------------------------------
    for fam, fname in (("c2iii", "fig_size_calibration.png"),
                       ("c2iv", None), ("c2ii", "fig_baseline_favorable.png")):
        p = ROOT / "results_c2" / f"{fam}.parquet"
        if not p.exists():
            continue
        d2 = pd.read_parquet(p)
        if fam == "c2iii":
            dd = d2[d2.method == "_diag"]
            size_tab = dd.assign(rej_zb=dd.z_boot_p < 0.05,
                                 rej_ztw=dd.z_tw_p < 0.05,
                                 rej_t=dd.trend_p < 0.05) \
                .groupby("cell_id")[["rej_zb", "rej_ztw", "rej_t"]].mean()
            ax = size_tab.plot.bar(figsize=(8, 3.8))
            ax.axhline(0.05, color="k", lw=0.7); ax.axhspan(0.03, 0.08,
                                                            color="orange", alpha=0.15)
            ax.set_ylabel("reject rate @ 5% nominal")
            plt.tight_layout(); plt.savefig(FIG / fname, dpi=150); plt.close(fig)
            out["size_calibration"] = size_tab.reset_index().to_dict("records")
        if fam == "c2ii":
            piv = d2[d2.method != "_diag"].pivot_table(
                index="cell_id", columns="method", values="rmse")
            piv.plot.bar(figsize=(8, 3.8))
            plt.ylabel("RMSE / sigma"); plt.title("dense weak factors")
            plt.tight_layout(); plt.savefig(FIG / fname, dpi=150); plt.close(fig)

    (FIG / "analysis_outputs.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "kink"},
                     indent=1)[:1200])


if __name__ == "__main__":
    main()
