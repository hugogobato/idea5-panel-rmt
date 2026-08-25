"""C5 analysis (runs after merge_shards.py validates the C5 families).

Produces figures/fig_c5*.png and figures/memo_c5_inputs.json with every
preregistration_c5_addendum verdict number.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "scripts"))
from scm_frontier import kink_breakpoint  # noqa: E402
from make_figures import frontier_rmse  # noqa: E402

FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


def main() -> None:
    out: dict[str, object] = {}

    # ---------------- C5a ---------------------------------------------------
    pa = ROOT / "results_c1" / "c5a_kink_confirm.parquet"
    if pa.exists():
        dfa = pd.read_parquet(pa)
        spec = dfa[dfa.method == "spectral_gated"]
        rows = []
        fig, axes = plt.subplots(1, 5, figsize=(21, 3.8), squeeze=False)
        for ax, (c, sub) in zip(axes[0], spec.groupby("c")):
            g = sub.groupby("m")["rmse"].agg(["mean", "sem"]).reset_index()
            ms, mu = g["m"].to_numpy(float), g["mean"].to_numpy(float)
            b = kink_breakpoint(ms, mu)
            T0 = int(sub["T0"].iloc[0])
            rows.append(dict(c=float(c), breakpoint=float(b),
                             within_15pct=bool(abs(b - 1.0) <= 0.15)))
            ax.errorbar(ms, mu, yerr=1.96 * g["sem"], fmt="o-", ms=3, lw=1,
                        color="#1f77b4")
            mm = np.linspace(ms.min(), ms.max(), 200)
            ax.plot(mm, [frontier_rmse(v, float(c), T0, 1.0) for v in mm],
                    "--", color="crimson")
            ax.axvline(1.0, color="k", lw=0.7, alpha=0.6)
            ax.axvline(b, color="orange", ls="-.", lw=1)
            ax.set_title(f"c={c}: b={b:.2f}")
        fig.tight_layout()
        fig.savefig(FIG / "fig_c5a_kink_confirm.png", dpi=140)
        plt.close(fig)
        share = float(np.mean([r["within_15pct"] for r in rows]))
        out["c5a"] = dict(table=rows, share=share,
                          criterion_pass=bool(share >= 0.8))

    # ---------------- C5b ---------------------------------------------------
    pb = ROOT / "results_c1" / "c5b_bite_extension.parquet"
    if pb.exists():
        dfb = pd.read_parquet(pb)
        incs = ["scm_simplex", "ridge_sc", "mc_nn_cv"]
        sub = dfb[(dfb.m <= 0.8)]
        plateaus = sub[sub.method.isin(incs)].groupby(
            ["method", "c"])["rmse"].mean().unstack().round(3)
        spec_sub = dfb[(dfb.method == "spectral_gated") & (dfb.m <= 0.8)]
        sup = dfb[(dfb.method == "_diag") & (dfb.m >= 1.5)]
        diag_sub = dfb[(dfb.method == "_diag") & (dfb.m <= 0.8)]
        flag_power = float((diag_sub.k_selected.astype(int) == 0).mean())
        false_alarm = float((sup.k_selected.astype(int) == 0).mean())
        crit1 = bool((plateaus >= 1.95).all().all())
        crit2 = bool(flag_power >= 0.80)
        crit3 = bool(false_alarm <= 0.20)
        piv = dfb[dfb.method != "_diag"].pivot_table(
            index=["c", "m"], columns="method", values="rmse")
        fig, ax = plt.subplots(figsize=(8.5, 4))
        for meth, col in (("spectral_gated", "#1f77b4"),
                          ("scm_simplex", "#7f7f7f"),
                          ("mc_nn_cv", "#2ca02c"), ("ridge_sc", "#9467bd")):
            gm = dfb[dfb.method == meth].groupby("m")["rmse"].mean()
            ax.plot(gm.index, gm.values, "o-", ms=3, lw=1.1, color=col,
                    label=meth)
        ax.axhline(2.0, color="crimson", ls="--", lw=1, label="2 sigma line")
        ax.axhline(np.sqrt(4.0), color="gray", ls=":", lw=1,
                   label="F plateau sqrt(1+theta)")
        ax.set_xlabel("spike multiplier m"); ax.set_ylabel("RMSE / sigma")
        ax.set_title("C5b bite extension: theta = 3 arm (pooled c)")
        ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(FIG / "fig_c5b_bite_extension.png",
                                        dpi=140); plt.close(fig)
        out["c5b"] = dict(plateau_means=plateaus.reset_index().to_dict("records"),
                          flag_power_subedge=flag_power,
                          false_alarm_supedge=false_alarm,
                          criterion=dict(crit1_all_incumbents_ge195=crit1,
                                         crit2_flag_power=crit2,
                                         crit3_false_alarm=crit3),
                          bite_pass=bool(crit1 and crit2 and crit3))

    # ---------------- C5c ---------------------------------------------------
    pc = ROOT / "results_raw" / "c5c" / "c5c_diagv2.csv"
    if pc.exists():
        dc = pd.read_csv(pc)
        dc["state"] = dc["state"].fillna("null")  # None -> empty field
        dc["rej_shift"] = dc.z_shift_p < 0.05
        dc["rej_perm"] = dc.z_perm_p < 0.05
        dc["rej_ztw"] = dc.z_tw_p < 0.05
        dc["rej_t"] = dc.trend_p < 0.05
        size = dc[dc.state == "null"].groupby("law")[
            ["rej_shift", "rej_perm", "rej_ztw", "rej_t",
             "gate_lrv_k", "gate_mp_k"]].mean().round(3)
        powr = dc[dc.state.str.startswith("break")].groupby(
            ["state", "law"])[["rej_shift", "rej_perm"]].mean().round(3)
        spiked = dc[dc.state == "spiked"].groupby("law")[["rej_shift"]].mean().round(3)
        laws_ok = ["ar1_r03", "ar1_r07", "het"]
        ok_size = bool(((size.loc[laws_ok, "rej_shift"] >= 0.03) &
                        (size.loc[laws_ok, "rej_shift"] <= 0.08)).all())
        ok_gate = bool((size["gate_lrv_k"] <= 0.06).all())
        det_d2 = powr.xs("break_d2", level="state")["rej_shift"]
        ok_pow = bool((det_d2 >= 0.80).all())
        fig, axes = plt.subplots(1, 2, figsize=(13, 4))
        size[["rej_shift", "rej_perm", "rej_ztw", "rej_t"]].plot.bar(
            ax=axes[0]); axes[0].axhline(0.05, color="k", lw=0.8)
        axes[0].axhspan(0.03, 0.08, color="orange", alpha=0.15)
        axes[0].set_title("size @5% by law"); axes[0].set_ylim(0, 1.05)
        powr.unstack("state")["rej_shift"].plot.bar(ax=axes[1])
        axes[1].axhline(0.80, color="crimson", ls="--", lw=1)
        axes[1].set_title("z_shift detection by break delta")
        fig.tight_layout(); fig.savefig(FIG / "fig_c5c_size_power.png",
                                        dpi=140); plt.close(fig)
        out["c5c"] = dict(size=size.reset_index().to_dict("records"),
                          power=powr.reset_index().to_dict("records"),
                          spiked_no_break=spiked.reset_index().to_dict("records"),
                          criterion=dict(size_window=ok_size, gate_lrv=ok_gate,
                                         detection_delta2=ok_pow),
                          pass_=bool(ok_size and ok_gate and ok_pow))

    # ---------------- C5d ---------------------------------------------------
    pd_path = ROOT / "results_raw" / "c5d" / "c5d_break_formal.csv"
    if pd_path.exists():
        dd = pd.read_csv(pd_path)
        dose = dd.assign(rej=dd.p_post < 0.05).groupby("cell_id").agg(
            power=("rej", "mean"), rmse=("rmse_spectral", "mean")).round(3)
        ok_pow = bool(dose.loc["delta2.0", "power"] >= 0.80) if "delta2.0" in dose.index else False
        ok_size = bool(dose.loc["delta0", "power"] <= 0.15) if "delta0" in dose.index else False
        fig, ax = plt.subplots(figsize=(6.5, 4))
        ax.bar(dose.index, dose.power, color="#1f77b4")
        ax.axhline(0.05, color="k", ls=":", lw=1)
        ax.axhspan(0.03, 0.08, color="orange", alpha=0.15)
        ax.set_title("C5d post-window break test, rejection @5%")
        fig.tight_layout(); fig.savefig(FIG / "fig_c5d_break.png", dpi=140)
        plt.close(fig)
        out["c5d"] = dict(dose=dose.reset_index().to_dict("records"),
                          criterion=dict(power_delta2=ok_pow, size_delta0=ok_size),
                          pass_=bool(ok_pow and ok_size))

    (FIG / "memo_c5_inputs.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:2400])


if __name__ == "__main__":
    main()
