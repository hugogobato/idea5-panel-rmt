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


def gap_statistic(
    basis_window: np.ndarray,
    post_window: np.ndarray,
    sigma: float,
    c_basis: float,
    k_max: int = 4,
) -> tuple[float, int]:
    """SELF-NORMALIZED post-residual gap (C5 deviation D5v2).

    g = lambda1 / median(eigenvalues) of the residual unit-space scatter.
    Dividing by the panel's own bulk center cancels the panel-level random
    scale, making the conditional (within-panel rotation) null coincide
    with the marginal law -- the unstudentized z failed exactly this way
    (conditional sd 0.76 vs marginal 1.31). Robust to sigma misspecification
    by construction; large g = factor-law instability evidence.
    """
    n_d, Tb = basis_window.shape
    Tp = post_window.shape[1]
    U, d, _ = np.linalg.svd(basis_window, full_matrices=False)
    evals_b = d**2 / Tb
    k = gated_rank(evals_b, sigma, c_basis, k_max)
    if k > 0:
        R = post_window - U[:, :k] @ (U[:, :k].T @ post_window)
    else:
        R = post_window
    ev = np.linalg.eigvalsh((R @ R.T) / Tp)
    med = float(np.median(ev))
    return float(ev[-1] / max(med, 1e-12)), int(k)


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

    KNOWN DEFECT (gate_g3_memo Section 3.2): with-replacement resampling
    duplicates blocks and inflates the null eigenvalue; superseded by
    z_shift_pvalue / z_perm_pvalue in the C5 addendum. Kept for the record.
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


# --------------------- C5 addendum instruments (frozen designs) ------------


def robust_row_scales(basis_window: np.ndarray) -> np.ndarray:
    """C5c het handling: per-row robust scale from first differences of the
    BASIS window only (leakage-safe): median(|dY|)/(0.6745*sqrt(2))."""
    d = np.diff(basis_window, axis=1)
    mad = np.median(np.abs(d), axis=1)
    return np.maximum(mad / (0.6745 * np.sqrt(2.0)), 1e-12)


def _split_windows(Y: np.ndarray, T_post: int):
    Tp = int(min(T_post, Y.shape[1] // 2))
    return Y[:, : Y.shape[1] - Tp], Y[:, Y.shape[1] - Tp :]


def z_shift_pvalue(
    donors_pre: np.ndarray,
    sigma: float,
    T_post: int,
    k_max: int = 4,
    B: int = 200,
    seed: int = 8_880_001,
    standardize: bool = True,
) -> tuple[float, float, int]:
    """Circular-shift null (C5 addendum C5c primary instrument).

    Rotates the time index by a uniform random offset; under H0 every split
    position is exchangeable and within-series dependence is preserved
    exactly. When `standardize`, EVERY draw (observed and rotated alike)
    is standardized by ITS OWN basis-window robust scales -- global
    one-shot standardization would break exchangeability (deviation D5).
    Returns (p_value, z_obs, k_used). Validity is claimed ONLY under H0.
    """
    Y = donors_pre
    n_d = Y.shape[0]
    Tb_w, Tp_w = _split_windows(Y, T_post)
    c_basis = n_d / Tb_w.shape[1]

    def stat(basis: np.ndarray, post: np.ndarray):
        # Per-draw standardization (deviation D5): scales from THIS draw's
        # basis window only; global standardization breaks exchangeability.
        if standardize:
            basis, post, sig = _std_windows(basis, post)
        else:
            sig = sigma
        return gap_statistic(basis, post, sig, c_basis, k_max)

    z_obs, k_used = stat(Tb_w, Tp_w)
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(B):
        o = int(rng.integers(Y.shape[1]))
        Bs, Ps = _split_windows(np.roll(Y, o, axis=1), T_post)
        z_star, _ = stat(Bs, Ps)
        if z_star >= z_obs:
            ge += 1
    return float((ge + 1.0) / (B + 1.0)), float(z_obs), int(k_used)


def z_perm_pvalue(
    donors_pre: np.ndarray,
    sigma: float,
    T_post: int,
    k_max: int = 4,
    B: int = 200,
    block: int = 20,
    seed: int = 8_880_001,
    standardize: bool = True,
) -> tuple[float, float, int]:
    """Disjoint-block PERMUTATION null (C5 secondary; no duplication).
    Standardization is per-draw (see z_shift_pvalue, deviation D5)."""
    Y = donors_pre
    n_d, T0 = Y.shape
    Tb_w, Tp_w = _split_windows(Y, T_post)
    c_basis = n_d / Tb_w.shape[1]

    def stat(basis: np.ndarray, post: np.ndarray):
        if standardize:
            basis, post, sig = _std_windows(basis, post)
        else:
            sig = sigma
        return gap_statistic(basis, post, sig, c_basis, k_max)

    z_obs, k_used = stat(Tb_w, Tp_w)
    rng = np.random.default_rng(seed)
    nb = int(np.ceil(T0 / block))
    starts = np.arange(nb) * block
    ge = 0
    for _ in range(B):
        order = rng.permutation(nb)
        idx = np.concatenate(
            [np.arange(starts[o], starts[o] + block) % T0 for o in order]
        )[:T0]
        Bs, Ps = _split_windows(Y[:, idx], T_post)
        z_star, _ = stat(Bs, Ps)
        if z_star >= z_obs:
            ge += 1
    return float((ge + 1.0) / (B + 1.0)), float(z_obs), int(k_used)


def _lrv_rows(Y: np.ndarray) -> np.ndarray:
    """Per-row Newey-West long-run variance (Bartlett, lag 4(T0/100)^(2/9)).

    NOTE (C5 deviation D5): biased low at short lags; retained only as a
    diagnostic. The gate uses level-variance inflation instead (see
    gate_lrv), which is exact for AR(1): phi * sigma^2 = Var(y).
    """
    n, T0 = Y.shape
    L = int(4.0 * (T0 / 100.0) ** (2.0 / 9.0))
    out = np.empty(n)
    for i, row in enumerate(Y):
        x = row - row.mean()
        ac0 = float(x @ x) / T0
        s = ac0
        for l in range(1, L + 1):
            w = 1.0 - l / (L + 1.0)
            s += 2.0 * w * float(x[l:] @ x[:-l]) / T0
        out[i] = max(s, 1e-12)
    return out


def gate_lrv(
    donors_pre: np.ndarray, k_max: int = 4, standardize: bool = True,
    tol: float = 1.05,
) -> int:
    """Silence gate v2 (C5c, amended design D5v2).

    Robust difference-based row-standardization removes cross-sectional
    heteroskedasticity; the MP-edge sigma^2 is then the median per-row
    long-run variance estimated as LRV_i = Vy_i * (2 Vy_i - Vd_i) / Vd_i,
    where Vy_i is the demeaned level variance and Vd_i the first-difference
    variance of row i. This estimator is EXACT for AR(1)
    (LRV = sigma_y^2 (1+rho)/(1-rho)) and reduces to Vy under iid, avoiding
    the Newey-West small-sample downward bias (deviation D5).
    Returns k = 0 when lambda1 <= tol*(1+sqrt(c))^2*sigma2_eff, else the
    largest-gap rank among the top k_max.
    """
    Y = donors_pre
    scales = robust_row_scales(Y)
    # Adaptive standardization (deviation D5v3): engage only when clear
    # cross-sectional scale heterogeneity is detected, so that under
    # (near-)homoskedastic laws the gate reduces to the known-sigma MP rule.
    # When engaged, use the Gaussian-efficient difference SD (low-noise
    # scaling); fall back to the robust median scale for heavy tails
    # (median/MAD-vs-SD discrepancy > 2x).
    tol_eff = tol
    if standardize and Y.shape[1] > 8 and np.quantile(scales, 0.9) > 2.0 * np.quantile(scales, 0.1):
        # Triggered path absorbs its own scale-estimation noise via a wider
        # tolerance (deviation D5v3).
        tol_eff = max(tol, 1.15)
        dsd = np.sqrt(np.mean(np.diff(Y, axis=1) ** 2, axis=1))
        mad_ratio = np.median(np.abs(np.diff(Y, axis=1))) / (
            0.6745 * np.maximum(dsd / np.sqrt(2.0), 1e-12))
        eff = np.where(mad_ratio > 2.0,
                       np.maximum(scales, 1e-12),
                       np.maximum(dsd / np.sqrt(2.0), 1e-12))
        Y = Y / eff[:, None]
    evals = scree(Y)
    c = Y.shape[0] / Y.shape[1]
    vy = np.var(Y, axis=1)
    vd = np.var(np.diff(Y, axis=1), axis=1)
    # Var(diff) = 2 sigma_y^2 (1 - rho) => rho_hat = 1 - vd/(2 vy);
    # LRV = sigma_y^2 (1+rho)/(1-rho) = vy (1+rho_hat)/(1-rho_hat).
    # Exact for AR(1); reduces to vy under iid.
    rho_hat = np.clip(1.0 - vd / (2.0 * np.maximum(vy, 1e-12)), -0.9, 0.98)
    phi = (1.0 + rho_hat) / (1.0 - rho_hat)
    lrv = vy * phi
    sigma2_eff = float(np.median(lrv))
    if evals[0] <= tol_eff * (1.0 + np.sqrt(c)) ** 2 * sigma2_eff:
        return 0
    return select_rank_gap_ratio(evals, k_max)


def _std_windows(basis: np.ndarray, post: np.ndarray):
    """Per-draw standardization: scales from THIS draw's basis window only,
    applied to both windows (preserves exchangeability across rotations)."""
    sc = robust_row_scales(basis)
    sigma_eff = float(np.median(sc))
    return basis / sc[:, None], post / sc[:, None], sigma_eff


def pre_trends_post_test(
    donors_pre: np.ndarray,
    donors_post: np.ndarray,
    sigma: float,
    G: int = 300,
    seed: int = 7_770_001,
    null_z: np.ndarray | None = None,
) -> dict:
    """C5d formalized break-detection statistic (preregistration_c5_addendum).

    resid_statistic on the FULL donor-pre basis vs the REAL donor-post
    window, calibrated by the simulated iid finite-n null. Gaussian-noise
    scope; leakage-safe (treated post never touched).
    """
    n_d, T0 = donors_pre.shape
    Tp = donors_post.shape[1]
    z, k = resid_statistic(donors_pre, donors_post, sigma, n_d / T0)
    if null_z is None:
        null_z = simulated_null_z(n_d, T0, Tp, sigma, n_d / T0, G=G, seed=seed)
    return {"z": float(z), "p": z_tw_pvalue(z, null_z), "k": int(k)}


def kink_breakpoint(ms: np.ndarray, rmse_mean: np.ndarray) -> float:
    """Amended primary kink estimator (C5a): plateau+linear two-segment fit.

    Returns argmin_b over arange(0.6, 1.6, 0.05) of SSE(constant | m<=b)
    + SSE(OLS line | m>b).
    """
    best_sse, best_b = np.inf, float("nan")
    for b in np.arange(0.60, 1.6001, 0.05):
        L = ms <= b
        R = ms > b
        if L.sum() < 2 or R.sum() < 3:
            continue
        sse = float(((rmse_mean[L] - rmse_mean[L].mean()) ** 2).sum())
        X = np.column_stack([np.ones(int(R.sum())), ms[R]])
        beta, *_ = np.linalg.lstsq(X, rmse_mean[R], rcond=None)
        sse += float(((rmse_mean[R] - X @ beta) ** 2).sum())
        if sse < best_sse:
            best_sse, best_b = sse, float(b)
    return best_b


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
