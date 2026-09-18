"""WP-D2 distance-to-frontier machinery (docs/preregistrations/preregistration_d2_addendum.md).

Every function implements a frozen rule from the addendum; thresholds are
constants here and must not be tuned after seeing real-panel outputs.
Leakage: only donor-pre data and the treated-pre row ever enter Sections
3-5 statistics; donor-post enters solely through pre_trends_post_test.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from scm_frontier.diagnostics import (
    cv_rank_selector,
    gate_lrv,
    gated_rank,
    invert_bbp,
    pre_trends_post_test,
    select_rank_gap_ratio,
    z_shift_pvalue,
)

GATE_TOL = 1.05
K_MAX = 4
ALIGN_P = 0.05
G_NULL = 2000
B_DONOR = 500
B_TIME = 200
TIME_BLOCK = 4
TRIM = 2
ALARM_MAX = 0.20
RHO_THRESH = {"smoking_prop99": -0.30, "german_reunification": -0.42}
N_PERM = 999
COMPARABLE_FRAC = 0.20


def center_rows(Y: np.ndarray) -> np.ndarray:
    return Y - Y.mean(axis=1, keepdims=True)


@lru_cache(maxsize=64)
def mp_median(c: float) -> float:
    """Median of the unit-mean Marchenko-Pastur law with parameter c."""
    sq = math.sqrt(c)
    a, b = (1.0 - sq) ** 2, (1.0 + sq) ** 2
    span = b - a

    def g(u: float) -> float:
        x = a + span * u * u
        return math.sqrt(max((b - x) * (x - a), 0.0)) / max(x, 1e-300) \
            * 2.0 * u * span

    total = quad(g, 0.0, 1.0, limit=200)[0]
    target = 0.5 * total

    def h(u: float) -> float:
        return quad(g, 0.0, u, limit=200)[0] - target

    u_star = brentq(h, 0.0, 1.0, xtol=1e-12)
    return float(a + span * u_star * u_star)


def spectrum_pipeline(evals: np.ndarray, n_d: int, T0: int) -> dict:
    """Addendum Section 3: bulk-refined sigma, gate, BBP inversion, d."""
    evals = np.maximum(np.asarray(evals, dtype=float), 0.0)
    c = n_d / T0
    floor = 1e-8 * max(float(evals[0]), 1e-12)
    sigma2_0 = max(float(np.median(evals) / mp_median(c)), floor)
    k0 = gated_rank(evals, math.sqrt(sigma2_0), c, K_MAX)
    if n_d - k0 >= max(10, math.ceil(n_d / 4)):
        sigma2_1 = max(float(np.mean(evals[k0:])), floor)
    else:
        sigma2_1 = sigma2_0
    k = int(gated_rank(evals, math.sqrt(sigma2_1), c, K_MAX))
    # BBP law is on sigma^2=1 scale: invert lam/sigma^2
    lam_scaled = evals[: max(k, 1)] / max(sigma2_1, floor)
    s_hat = invert_bbp(lam_scaled, c)[:k] if k > 0 else np.array([])
    m_hat = s_hat / math.sqrt(c) if k > 0 else np.array([])
    valid = m_hat[np.isfinite(m_hat)] if k > 0 else np.array([])
    d = float(valid.max()) if valid.size else 0.0
    return {
        "c": c, "sigma2_0": sigma2_0, "sigma2_1": sigma2_1, "k": k,
        "s_hat": [float(x) for x in s_hat],
        "m_hat": [float(x) for x in m_hat],
        "edge": sigma2_1 * (1.0 + math.sqrt(c)) ** 2, "d": d,
    }


def _energy(Vt: np.ndarray, y1: np.ndarray, k: int) -> float:
    if k == 0:
        return 0.0
    proj = Vt[:k] @ y1
    return float(proj @ proj / (y1 @ y1))


def alignment_test(
    Vt: np.ndarray, y1: np.ndarray, k: int, pool: np.ndarray
) -> tuple[float, float, float]:
    """Addendum Section 4: scale-free pure-noise alignment null."""
    e_obs = _energy(Vt, y1, k)
    e_null = ((pool @ Vt[:k].T) ** 2).sum(axis=1) / (pool ** 2).sum(axis=1)
    p = float((int((e_null >= e_obs).sum()) + 1.0) / (len(e_null) + 1.0))
    r_align = float(e_obs / max(np.median(e_null), 1e-12))
    return e_obs, p, r_align


def classify(k: int, p_align: float, q05: float, q95: float) -> str:
    """Addendum Section 5 labels."""
    if k == 0:
        return "FRAGILE-INVISIBLE"
    if p_align >= ALIGN_P:
        return "FRAGILE-MISALIGNED"
    if q05 >= 1.0:
        return "RECOVERABLE"
    if q95 < 1.0:
        return "FRAGILE-SUBEDGE"
    return "INCONCLUSIVE-SUBEDGE"


def _block_indices(rng: np.random.Generator, T0: int, block: int) -> np.ndarray:
    starts = rng.integers(0, T0, size=int(math.ceil(T0 / block)))
    return np.concatenate(
        [np.arange(s, s + block) % T0 for s in starts])[:T0]


def _select_rank(selector: str, evals: np.ndarray, Ycen: np.ndarray) -> int:
    if selector == "gate_lrv":
        return int(gate_lrv(Ycen, K_MAX))
    c = Ycen.shape[0] / Ycen.shape[1]
    sigma2 = float(np.median(evals) / mp_median(c))
    return int(gated_rank(evals, math.sqrt(sigma2), c, K_MAX))


def analyze_space(
    Yd_raw: np.ndarray,
    y1_raw: np.ndarray,
    pool: np.ndarray,
    seeds: dict,
    idx: int,
    selector: str = "primary",
    with_bootstrap: bool = True,
) -> dict:
    """Full pipeline of Sections 3-6 on one (space, window) configuration."""
    Yc = center_rows(Yd_raw)
    y1c = center_rows(y1_raw[None, :])[0]
    n_d, T0 = Yc.shape
    evals = np.sort(np.linalg.eigvalsh((Yc @ Yc.T) / T0))[::-1]

    def run(ev: np.ndarray, mat: np.ndarray) -> dict:
        sp = spectrum_pipeline(ev, mat.shape[0], mat.shape[1])
        if selector != "primary":
            k_alt = _select_rank(selector, ev, mat)
            sp["k"] = k_alt
            # scale by bulk variance for BBP inversion
            floor_alt = 1e-8 * max(float(ev[0]), 1e-12) if len(ev) else 1e-12
            lam_s = ev[:max(k_alt, 1)] / max(sp["sigma2_1"], floor_alt)
            sh = invert_bbp(lam_s, sp["c"])[:k_alt] \
                if k_alt > 0 else np.array([])
            mh = sh / math.sqrt(sp["c"]) if k_alt > 0 else np.array([])
            vh = mh[np.isfinite(mh)] if k_alt > 0 else np.array([])
            sp["d"] = float(vh.max()) if vh.size else 0.0
        return sp

    spec = run(evals, Yc)
    _, _, Vt = np.linalg.svd(Yc, full_matrices=False)
    e_obs, p_align, r_align = alignment_test(Vt, y1c, spec["k"], pool)
    out = {**spec, "e_obs": e_obs, "p_align": p_align, "r_align": r_align,
           "evals_top12": [float(x) for x in evals[:12]]}
    if not with_bootstrap:
        return out

    rng_d = np.random.default_rng(seeds["donor"] + idx * 10**6)
    d_boot, k_boot = [], []
    for _ in range(B_DONOR):
        rid = rng_d.integers(0, n_d, n_d)
        Yb = Yc[rid]
        ev = np.sort(np.linalg.eigvalsh((Yb @ Yb.T) / T0))[::-1]
        sp = run(ev, Yb)
        _, _, Vtb = np.linalg.svd(Yb, full_matrices=False)
        _, pb, _ = alignment_test(Vtb, y1c, sp["k"], pool)
        d_boot.append(sp["d"])
        k_boot.append(sp["k"])
    q05, q50, q95 = np.percentile(d_boot, [5, 50, 95])

    rng_t = np.random.default_rng(seeds["time"] + idx * 10**6)
    d_tb = []
    for _ in range(B_TIME):
        cols = _block_indices(rng_t, T0, TIME_BLOCK)
        Yb, y1b = Yc[:, cols], y1c[cols]
        ev = np.sort(np.linalg.eigvalsh((Yb @ Yb.T) / len(cols)))[::-1]
        d_tb.append(run(ev, Yb)["d"])
    t05, _, t95 = np.percentile(d_tb, [5, 50, 95])

    out.update({
        "boot_donor": {"q05": float(q05), "q50": float(q50),
                       "q95": float(q95),
                       "share_k0": float(np.mean(np.array(k_boot) == 0))},
        "boot_time": {"q05": float(t05), "q50": float(np.median(d_tb)),
                      "q95": float(t95)},
        "label": classify(spec["k"], p_align, float(q05), float(q95)),
    })
    return out


def analyze_unit(
    Yd_pre: np.ndarray,
    y1_pre: np.ndarray,
    seeds: dict,
    idx: int,
) -> dict:
    """Primary space + diff-space decomposition + Section 9 sensitivity."""
    T0 = y1_pre.shape[0]
    rng_a = np.random.default_rng(seeds["align"] + idx * 10**6)
    pool = rng_a.standard_normal((G_NULL, T0))

    prim = analyze_space(Yd_pre, y1_pre, pool, seeds, idx)

    D = np.diff(np.vstack([Yd_pre, y1_pre[None, :]]), axis=1)
    pool_d = np.random.default_rng(
        seeds["align"] + idx * 10**6 + 1).standard_normal((G_NULL, T0 - 1))
    dif = analyze_space(D[:-1], D[-1], pool_d, seeds, idx,
                        with_bootstrap=False)

    sens = {}
    for win in ["full", "trimmed"]:
        sl = slice(None) if win == "full" else slice(TRIM, T0 - TRIM)
        Yw, y1w = Yd_pre[:, sl], y1_pre[sl]
        Tw = y1w.shape[0]
        pool_w = np.random.default_rng(
            seeds["align"] + idx * 10**6 + 2).standard_normal((G_NULL, Tw))
        for sel in ["primary", "gate_lrv"]:
            sens[f"{sel}_{win}"] = analyze_space(
                Yw, y1w, pool_w, seeds, idx, selector=sel)

    base_label = sens["primary_full"]["label"]
    agree = sum(
        sens[key]["label"] == base_label
        for key in ["primary_full", "gate_lrv_full",
                    "primary_trimmed", "gate_lrv_trimmed"])
    prim.update({
        "diff_space": {kk: dif[kk] for kk in
                       ["k", "d", "p_align", "sigma2_1"]},
        "sensitivity": {key: {kk: sens[key][kk] for kk in
                              ["k", "d", "label"]}
                        for key in sens},
        "stability_agree": int(agree),
        "stability_pass": bool(agree >= 3),
        "cv_rank_k": int(cv_rank_selector(center_rows(Yd_pre),
                                          center_rows(y1_pre[None, :])[0],
                                          K_MAX)),
        "ungated_k": int(select_rank_gap_ratio(
            np.sort(np.linalg.eigvalsh(
                (center_rows(Yd_pre) @ center_rows(Yd_pre).T) / T0))[::-1],
            K_MAX)),
    })
    return prim


def spearman_with_perm(
    x: np.ndarray, y: np.ndarray, seed: int
) -> tuple[float, float]:
    """Spearman rho plus one-sided permutation p (H1: rho < 0)."""
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rho = float(np.corrcoef(rx, ry)[0, 1])
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(N_PERM):
        rz = np.corrcoef(rx, ry[rng.permutation(len(ry))])[0, 1]
        if rz <= rho:
            ge += 1
    return rho, float((ge + 1.0) / (N_PERM + 1.0))


@dataclass
class Seeds:
    align: int
    donor: int
    time: int
    perm: int

    def as_dict(self) -> dict:
        return {"align": self.align, "donor": self.donor, "time": self.time}


def date_alarm_scan(
    Yd_pre: np.ndarray, cutoff_idx: list[int], seed: int, sigma2: float
) -> list[dict]:
    """P3: shipped z_shift instrument at pseudo-cutoffs inside the pre window."""
    rows = []
    for ci in cutoff_idx:
        p, _, _ = z_shift_pvalue(center_rows(Yd_pre[:, :ci]),
                                 math.sqrt(sigma2),
                                 T_post=6, B=200, seed=seed)
        rows.append({"cutoff_index": int(ci), "p": float(p),
                     "alarm": bool(p < ALIGN_P)})
    return rows


def posttest_control(
    Yd_pre: np.ndarray, Yd_post: np.ndarray, sigma2: float, seed: int
) -> dict:
    """Descriptive donor-post factor-law control (addendum Section 8)."""
    mu = Yd_pre.mean(axis=1, keepdims=True)
    res = pre_trends_post_test(center_rows(Yd_pre), center_rows(Yd_post - mu),
                               math.sqrt(sigma2), G=300, seed=seed)
    return {"z": res["z"], "p": res["p"], "k": res["k"]}


def nf_arms(
    treated: dict, treated_rmse: float, placebo_units: list[str],
    placebo_d: np.ndarray, placebo_rmse: np.ndarray,
    placebo_labels: list[str],
) -> dict:
    """Predeclared novel-finding arms (addendum Section 10)."""
    n1 = treated["label"].startswith("FRAGILE")
    n2 = treated["label"] == "RECOVERABLE" and treated.get("poor_looking", False)
    comp = np.abs(placebo_rmse - treated_rmse) <= COMPARABLE_FRAC * abs(treated_rmse)
    n3 = bool(comp.any()) and (
        treated["label"] == "RECOVERABLE"
        and any(l.startswith("FRAGILE") for l, c in zip(placebo_labels, comp) if c)
        or treated["label"].startswith("FRAGILE")
        and any(l == "RECOVERABLE" for l, c in zip(placebo_labels, comp) if c))
    q90 = float(np.percentile(placebo_d, 90)) if len(placebo_d) else float("nan")
    p2 = bool(treated["d"] >= q90) if len(placebo_d) else False
    return {
        "N1_fragile_treated": bool(n1),
        "N2_certified_safe": bool(n2),
        "N3_certification_asymmetry": n3,
        "NF": bool(n1 or n2 or n3),
        "P2_treated_separation": p2,
        "P2_q90_placebo_d": q90,
        "n_comparable_placebos": int(comp.sum()),
    }
