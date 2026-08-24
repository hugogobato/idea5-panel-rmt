"""Phase C cell runner (shared by all Colab notebooks and local validation).

Row contract (results_schema.yaml): one row per (rep, method) plus one
"_diag" row per rep carrying the diagnostic battery. All estimators see the
same generated panel within a replication; seeds follow preregistration
Section 3 (seed = 10000 + rep index).

Diagnostic levels (preregistration Sections 4-5):
  "none"  : estimators only.
  "light" : + gated/ungated rank, CV-rank comparator, trend t-test, scree.
  "full"  : light + Z_boot (block bootstrap) and Z_tw (simulated null,
            cached once per cell config).
"""

from __future__ import annotations

import time

import numpy as np

from .diagnostics import (
    cv_rank_selector,
    classical_trend_ttest,
    simulated_null_z,
    spike_estimates,
    z_boot_pvalue,
    z_tw_pvalue,
)
from .dgps import generate_panel
from .estimators import (
    donor_mean,
    mc_nn_cv,
    ridge_sc,
    scm_simplex,
    sdid,
    spectral_sc_full,
)

METHODS = (
    "donor_mean",
    "scm_simplex",
    "ridge_sc",
    "spectral_gated",
    "spectral_ungated",
    "mc_nn_cv",
    "sdid",
)

META_COLS = (
    "experiment", "shard", "cell_id", "c", "n", "T0", "T_post", "r", "m",
    "arm", "theta", "delta", "noise", "rho", "het_ratio", "rep_seed",
)
ROW_COLS = META_COLS + (
    "method", "rmse", "att_bias", "k_selected", "ridge_lambda", "mc_lambda",
    "solver_ok", "ci_cover", "z_boot_p", "z_tw_p", "trend_p", "cv_rank",
    "lambda1", "edge", "wall_ms",
)


def _meta(cell: dict, seed: int) -> dict:
    return {
        "experiment": cell["experiment"],
        "shard": cell.get("shard", ""),
        "cell_id": cell["cell_id"],
        "c": cell["c"],
        "n": cell["n"],
        "T0": cell["T0"],
        "T_post": cell["T_post"],
        "r": cell["r"],
        "m": cell["m"],
        "arm": cell["arm"],
        "theta": cell["theta"],
        "delta": cell.get("delta") if cell.get("delta") is not None else "",
        "noise": cell.get("noise", "gaussian"),
        "rho": cell.get("rho", 0.5),
        "het_ratio": cell.get("het_ratio", 4.0),
        "rep_seed": seed,
    }


def panel_kwargs(cell: dict) -> dict:
    arm = cell["arm"]
    if arm == "orthogonal":
        alignment = "none"
        share = tuple([0.0] * cell["r"])
    elif arm == "spread":
        alignment = "all"
        share = tuple([cell["theta"] / cell["r"]] * cell["r"])
    else:  # "partial" / "full" -> concentrated on factor 0
        alignment = "first"
        share = tuple([cell["theta"]] * cell["r"])
    if "spike_strengths" in cell:
        strengths = tuple(cell["spike_strengths"])
    else:
        s = cell["m"] * np.sqrt(cell["c"])
        strengths = tuple([float(s)] * cell["r"])
    return dict(
        n=cell["n"], T0=cell["T0"], T_post=cell["T_post"], r=cell["r"],
        spike_strengths=strengths,
        treated_share=share,
        alignment=alignment, sigma=1.0,
        noise=cell.get("noise", "gaussian"),
        rho=cell.get("rho", 0.5),
        het_ratio=cell.get("het_ratio", 4.0),
        structural_break=cell.get("delta"),
    )


def _empty_diag() -> dict:
    return {"z_boot_p": "", "z_tw_p": "", "trend_p": "", "cv_rank": "",
            "lambda1": "", "edge": ""}


def run_rep(cell: dict, seed: int, methods: tuple[str, ...],
            diag_level: str = "none", null_z: np.ndarray | None = None) -> list[dict]:
    """One replication -> schema rows.

    `null_z` is the per-cell cached simulated null for Z_tw; it is required
    when diag_level == "full" (compute via simulated_null_z once per cell).
    """
    pan = generate_panel(seed=seed, **panel_kwargs(cell))
    Y = pan["Y"]
    T0, T_post = cell["T0"], cell["T_post"]
    donors, y1_pre, y_star = Y[1:], Y[0, :T0], Y[0, T0:]
    c = float(cell["c"])
    sigma = 1.0
    meta = _meta(cell, seed)
    att_true = float(y_star.mean())
    rows: list[dict] = []

    need_spectral = any(m.startswith("spectral") for m in methods) or diag_level != "none"
    spec = spectral_sc_full(donors, y1_pre, sigma=sigma, c=c) if need_spectral else None

    for mname in methods:
        t0 = time.perf_counter()
        rl: object = ""
        ml: object = ""
        ok: object = ""
        k_sel: object = ""
        cover: object = ""
        if mname == "spectral_gated":
            pred, lo, hi = spec["pred_gated"], spec["ci_lo"], spec["ci_hi"]
            k_sel = spec["k_gated"]
            cover = float(np.mean((y_star >= lo) & (y_star <= hi)))
        elif mname == "spectral_ungated":
            pred = spec["pred_ungated"]
            k_sel = spec["k_ungated"]
        else:
            fn = {"donor_mean": donor_mean, "scm_simplex": scm_simplex,
                  "ridge_sc": ridge_sc, "mc_nn_cv": mc_nn_cv, "sdid": sdid}[mname]
            info: dict = {}
            pred = fn(donors, y1_pre, info=info)
            rl = info.get("ridge_lambda", "")
            ml = info.get("mc_lambda", "")
            ok = info.get("scm_success", info.get("sdid_solver_ok", ""))
        wall = (time.perf_counter() - t0) * 1000.0
        rows.append({**meta, **_empty_diag(), "method": mname,
                     "rmse": round(float(np.sqrt(np.mean((pred - y_star) ** 2))), 6),
                     "att_bias": round(float(np.mean(pred) - att_true), 6),
                     "k_selected": k_sel, "ridge_lambda": rl, "mc_lambda": ml,
                     "solver_ok": ok, "ci_cover": cover, "wall_ms": round(wall, 3)})

    if diag_level != "none":
        t0 = time.perf_counter()
        dg = _empty_diag()
        dg["cv_rank"] = cv_rank_selector(donors[:, :T0], y1_pre)
        dg["trend_p"] = round(classical_trend_ttest(y1_pre), 6)
        sc = spike_estimates(donors[:, :T0], sigma, c)
        dg["lambda1"] = round(float(sc["eigs"][0]), 6)
        dg["edge"] = round(float((1 + np.sqrt(c)) ** 2), 6)
        dg["k_selected"] = spec["k_gated"] if spec is not None else ""
        if diag_level == "full":
            pb, zb, _ = z_boot_pvalue(donors[:, :T0], sigma, T_post)
            dg["z_boot_p"] = round(pb, 6)
            dg["z_tw_p"] = round(z_tw_pvalue(zb, null_z), 6)
        rows.append({**meta, **dg, "method": "_diag", "rmse": "",
                     "att_bias": "", "ridge_lambda": "", "mc_lambda": "",
                     "solver_ok": "", "ci_cover": "",
                     "wall_ms": round((time.perf_counter() - t0) * 1000.0, 3)})
    return rows


def cell_null_z(cell: dict, G: int = 300, seed: int = 7_770_001) -> np.ndarray:
    """Simulated iid null of Z for this cell's geometry (computed once)."""
    n_d = cell["n"] - 1
    Tp_eff = int(min(cell["T_post"], cell["T0"] // 2))
    Tb = cell["T0"] - Tp_eff
    return simulated_null_z(n_d, Tb, Tp_eff, 1.0, n_d / Tb, G=G, seed=seed)


def run_cell(cell: dict, methods: tuple[str, ...], reps: int, seed_base: int,
             on_chunk, diag_level: str = "none") -> None:
    """Run a full cell, flushing rows to `on_chunk(rows)` every 25 reps."""
    null_z = cell_null_z(cell) if diag_level == "full" else None
    done = 0
    while done < reps:
        take = min(25, reps - done)
        chunk: list[dict] = []
        for j in range(done, done + take):
            chunk.extend(run_rep(cell, seed_base + j, methods, diag_level, null_z))
        on_chunk(chunk)
        done += take
