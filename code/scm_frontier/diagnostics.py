"""Diagnostic suite for the Idea 5 project (WP-B2).

Spectral diagnostics use pre-period data only (leakage rule). Calibration
conventions follow model_card.md: unit-space scatter, sigma^2 units,
c = n / T0.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from .estimators import select_rank_gap_ratio, unit_scatter


def scree(donors_pre: np.ndarray) -> np.ndarray:
    """Descending eigenvalues of the unit-space scatter."""
    return np.sort(np.linalg.eigvalsh(unit_scatter(donors_pre)))[::-1]


def tw_mu_nu(rows: int, cols: int) -> tuple[float, float]:
    """Johnstone (2001) centering/scaling for the largest eigenvalue of a
    rows x cols iid-Gaussian matrix on the (1/cols)-scaled convention."""
    mu = (np.sqrt(rows - 1.0) + np.sqrt(cols)) ** 2 / cols
    nu = (
        (np.sqrt(rows - 1.0) + np.sqrt(cols))
        / cols
        * (1.0 / np.sqrt(rows - 1.0) + 1.0 / np.sqrt(cols)) ** (1.0 / 3.0)
    )
    return float(mu), float(nu)


def tw_statistic(donors_pre: np.ndarray, sigma: float) -> float:
    """Standardized largest eigenvalue (Tracy-Widom scaling, Johnstone 2001).

    For an n_d x T0 iid-Gaussian matrix, (lambda1 - mu)/(nu) converges to TW1,
    with mu/nu computed on the (1/T0)-scaled convention used here.
    """
    n_d, T0 = donors_pre.shape
    evals = scree(donors_pre)
    mu, nu = tw_mu_nu(n_d, T0)
    return float((evals[0] / sigma**2 - mu) / nu)


def gated_rank(evals_desc: np.ndarray, sigma: float, c: float, k_max: int = 4) -> int:
    """Silence gate (top eigenvalue <= 1.05 x MP edge => 0) plus gap ratio."""
    edge = sigma**2 * (1.0 + np.sqrt(c)) ** 2
    if evals_desc[0] <= 1.05 * edge:
        return 0
    return select_rank_gap_ratio(evals_desc, k_max)


def cv_rank_selector(
    donors_pre: np.ndarray,
    y1_pre: np.ndarray,
    k_max: int = 4,
    folds: int = 4,
    cv_seed: int = 1234,
) -> int:
    """CV-rank incumbent comparator (preregistration Section 5.2).

    Chooses k in 0..k_max minimizing contiguous-block CV MSE of the k-PC
    regression of y1_pre on donor scores (same fold stream as ridge_cv_seed).
    """
    T0 = len(y1_pre)
    U, d, Vt = np.linalg.svd(donors_pre, full_matrices=False)
    km = int(min(k_max, len(d)))
    Z_all = Vt * d[:, None]  # k x T0 score series
    rng = np.random.default_rng(cv_seed)
    perm = rng.permutation(T0)
    blocks = np.array_split(perm, folds)
    best_k, best_err = 0, np.inf
    for k in range(0, km + 1):
        err = 0.0
        Z = Z_all[:k]  # k x T0 score series
        for b in blocks:
            tr = np.setdiff1d(perm, b)
            if k == 0:
                mu_tr = y1_pre[tr].mean()
                pred = np.full(len(b), mu_tr)
            else:
                Xt = Z[:, tr]
                beta, *_ = np.linalg.lstsq(Xt.T, y1_pre[tr], rcond=None)
                pred = Z[:, b].T @ beta
            err += float(np.sum((pred - y1_pre[b]) ** 2))
        if err < best_err - 1e-12:
            best_err, best_k = err, k
    return int(best_k)


def classical_trend_ttest(y1_pre: np.ndarray) -> float:
    """Classical incumbent: two-sided OLS t-test of the linear-trend
    coefficient of the treated pre-period row (nominal level reference)."""
    T0 = len(y1_pre)
    t = np.arange(T0, dtype=float)
    X = np.column_stack([np.ones(T0), t])
    beta, *_ = np.linalg.lstsq(X, y1_pre, rcond=None)
    resid = y1_pre - X @ beta
    dof = T0 - 2
    s2 = float(resid @ resid) / dof
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
    return float(2.0 * stats.t.sf(abs(beta[1] / se), dof))


def resid_statistic(
    basis_window: np.ndarray,
    post_window: np.ndarray,
    sigma: float,
    c_basis: float,
    k_max: int = 4,
) -> tuple[float, int]:
    """Post-residual TW statistic Z (preregistration Section 8.8).

    Projects the post window off the gated spike basis of the basis window;
    returns the standardized top eigenvalue of the residual unit-space
    scatter, (lam1/sigma^2 - mu)/nu with mu/nu at aspect n_d/Tp, plus the
    gated rank used. Large Z = factor-law instability evidence.
    """
    n_d, Tb = basis_window.shape
    Tp = post_window.shape[1]
    U, d, _ = np.linalg.svd(basis_window, full_matrices=False)
    evals = d**2 / Tb
    k = gated_rank(evals, sigma, c_basis, k_max)
    if k > 0:
        R = post_window - U[:, :k] @ (U[:, :k].T @ post_window)
    else:
        R = post_window
    lam = float(np.linalg.eigvalsh((R @ R.T) / Tp)[-1]) / sigma**2
    mu, nu = tw_mu_nu(n_d, Tp)
    return float((lam - mu) / nu), int(k)


def simulated_null_z(
    n_d: int,
    Tb: int,
    Tp: int,
    sigma: float,
    c_basis: float,
    k_max: int = 4,
    G: int = 300,
    seed: int = 7_770_001,
) -> np.ndarray:
    """Sorted simulated finite-n null of the Z statistic (iid Gaussian),
    reproducible given (shapes, seed); cache ONE draw set per cell config."""
    rng = np.random.default_rng(seed)
    zs = np.empty(G)
    for g in range(G):
        B = rng.normal(0.0, sigma, size=(n_d, Tb))
        P = rng.normal(0.0, sigma, size=(n_d, Tp))
        zs[g], _ = resid_statistic(B, P, sigma, c_basis, k_max)
    return np.sort(zs)


def z_tw_pvalue(z_obs: float, null_z_sorted: np.ndarray) -> float:
    """Parametric-simulation p-value: P(null Z >= observed), +1 corrected."""
    G = len(null_z_sorted)
    exceed = int(np.searchsorted(null_z_sorted, z_obs, side="left"))
    return float((G - exceed + 1.0) / (G + 1.0))


def z_boot_pvalue(
    donors_pre: np.ndarray,
    sigma: float,
    T_post: int,
    k_max: int = 4,
    B: int = 200,
    block: int = 10,
    seed: int = 8_880_001,
) -> tuple[float, float, int]:
    """Circular block-bootstrap pre-trends test (preregistration Section 8.8).

    Observed statistic splits the pre window into a basis part (first T0 -
    Tp_eff columns) and a pseudo-post part (last Tp_eff), Tp_eff =
    min(T_post, T0 // 2). Null: B circular block resamples (length `block`)
    of the pre time index, same split. Returns (p_value, z_obs, k_used).
    """
    Y = donors_pre
    n_d, T0 = Y.shape
    Tp = int(min(T_post, T0 // 2))
    Tb = T0 - Tp
    z_obs, k_used = resid_statistic(Y[:, :Tb], Y[:, Tb:], sigma, n_d / Tb, k_max)
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(B):
        starts = rng.integers(0, T0, size=int(np.ceil(T0 / block)))
        idx = np.concatenate([(np.arange(s, s + block) % T0) for s in starts])[:T0]
        Ystar = Y[:, idx]
        z_star, _ = resid_statistic(Ystar[:, :Tb], Ystar[:, Tb:], sigma, n_d / Tb, k_max)
        if z_star >= z_obs:
            ge += 1
    return float((ge + 1.0) / (B + 1.0)), float(z_obs), int(k_used)


def invert_bbp(lam: np.ndarray, c: float):
    """Map outlier location(s) back to spike strengths via the BBP/BGN law.

    Solves lambda = 1 + s + c + c/s for s; values at or below the bulk edge
    map to nan. Vectorized over lam.
    """
    lam = np.asarray(lam, dtype=float)
    edge = (1.0 + np.sqrt(c)) ** 2
    b = lam - 1.0 - c
    disc = b**2 - 4.0 * c
    s = np.where(disc > 0.0, (b + np.sqrt(np.maximum(disc, 0.0))) / 2.0, np.nan)
    return np.where(lam > edge, s, np.nan)


def alignment_energy(donors_pre: np.ndarray, y1_pre: np.ndarray, k: int) -> float:
    """Fraction of the treated row's pre-period energy captured by the top-k
    sample subspace, in excess of the pure-noise expectation sigma^2-normalized.

    Returns proj^2 energy share in [0, 1]: sum_{j<=k} <y1, v_j>^2 / ||y1||^2.
    """
    _, d, Vt = np.linalg.svd(donors_pre, full_matrices=False)
    kk = min(k, len(d))
    proj = Vt[:kk] @ y1_pre
    return float(proj @ proj / (y1_pre @ y1_pre))


def spike_estimates(donors_pre: np.ndarray, sigma: float, c: float, k_max: int = 4) -> dict:
    """Estimate supercritical spike strengths and the m = s/sqrt(c) ratios."""
    evals = scree(donors_pre)[: max(1, k_max)]
    s_hat = invert_bbp(evals, c)
    m_hat = s_hat / np.sqrt(c)
    return {"eigs": evals, "s_hat": s_hat, "m_hat": m_hat}
