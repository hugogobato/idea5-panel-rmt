"""Estimator library for the Idea 5 project (WP-B2).

Leakage rule (G0): estimators receive ONLY the donor matrix over the full
timeline and the treated pre-period row. The treated post-period outcomes are
never passed in, so leakage is structurally impossible; the test suite
additionally verifies behavioral invariance to treated-post perturbations.

Estimand: realized outcome y*_t = L_0,t + E_0,t for t in the post window.
The oracle therefore predicts L_0,post exactly and attains RMSE = sigma.

All estimators share the signature f(donors, y1_pre, *, info=None, **params)
and return a length-T_post vector of predictions. `info`, when given, is an
optional dict filled with method-specific diagnostics.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

DEFAULT_RIDGE_LAMBDAS = tuple(np.logspace(-1.0, 4.0, 11))
DEFAULT_MC_LAMBDAS = tuple(np.logspace(0.0, 3.0, 11))


def donor_mean(donors: np.ndarray, y1_pre: np.ndarray, info: dict | None = None) -> np.ndarray:
    """Uniform-weight donor average, one number per post period."""
    return donors[:, y1_pre.shape[0] :].mean(axis=0)


def oracle_predict(L_post_row: np.ndarray) -> np.ndarray:
    """Oracle uses the true latent trajectory of the treated unit."""
    return L_post_row


def scm_simplex(
    donors: np.ndarray,
    y1_pre: np.ndarray,
    maxiter: int = 500,
    ftol: float = 1e-12,
    info: dict | None = None,
) -> np.ndarray:
    """Abadie SC: least squares under w >= 0, sum(w) = 1 via SLSQP."""
    T0 = y1_pre.shape[0]
    X = donors[:, :T0].T
    n_d = X.shape[1]

    def obj(w):
        r = X @ w - y1_pre
        return float(r @ r)

    def jac(w):
        return 2.0 * X.T @ (X @ w - y1_pre)

    w0 = np.full(n_d, 1.0 / n_d)
    res = minimize(
        obj, w0, jac=jac, method="SLSQP",
        bounds=[(0.0, 1.0)] * n_d,
        constraints=[
            {"type": "eq", "fun": lambda w: w.sum() - 1.0, "jac": lambda w: np.ones(n_d)}
        ],
        options={"maxiter": maxiter, "ftol": ftol},
    )
    if info is not None:
        feasible = abs(res.x.sum() - 1.0) < 1e-6 and (res.x > -1e-8).all()
        info["scm_success"] = bool(res.success or (res.status == 8 and feasible))
        info["scm_status"] = int(res.status)
    return res.x @ donors[:, T0:]


def _cv_ridge_lambda(X, y, lambdas, folds, seed):
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    perm = rng.permutation(n)
    blocks = np.array_split(perm, folds)
    eye = np.eye(X.shape[1])
    best_lam, best_err = lambdas[0], np.inf
    for lam in lambdas:
        err = 0.0
        for b in blocks:
            tr = np.setdiff1d(perm, b)
            mux = X[tr].mean(axis=0)
            muy = y[tr].mean()
            Xt = X[tr] - mux
            yt = y[tr] - muy
            coef = np.linalg.solve(Xt.T @ Xt + lam * eye, Xt.T @ yt)
            pred = (X[b] - mux) @ coef + muy
            err += float(np.sum((pred - y[b]) ** 2))
        if err < best_err:
            best_err, best_lam = err, lam
    return best_lam


def ridge_sc(
    donors: np.ndarray,
    y1_pre: np.ndarray,
    lambdas=DEFAULT_RIDGE_LAMBDAS,
    folds: int = 4,
    cv_seed: int = 1234,
    info: dict | None = None,
) -> np.ndarray:
    """Ridge SC: treated row on donor rows with intercept; penalty by CV over
    pre periods only."""
    T0 = y1_pre.shape[0]
    X = donors[:, :T0].T
    lam = _cv_ridge_lambda(X, y1_pre, list(lambdas), folds, cv_seed)
    mux, muy = X.mean(axis=0), y1_pre.mean()
    Xt, yt = X - mux, y1_pre - muy
    coef = np.linalg.solve(Xt.T @ Xt + lam * np.eye(X.shape[1]), Xt.T @ yt)
    if info is not None:
        info["ridge_lambda"] = float(lam)
    return (donors[:, T0:].T - mux) @ coef + muy


def unit_scatter(Y_pre: np.ndarray) -> np.ndarray:
    """(1/T0) * Y Y' in unit space; model-card calibration convention."""
    T0 = Y_pre.shape[1]
    return (Y_pre @ Y_pre.T) / T0


def select_rank_gap_ratio(eigs_desc: np.ndarray, k_max: int) -> int:
    """Rank by largest successive eigenvalue-gap ratio, k in 1..k_max."""
    km = min(k_max, len(eigs_desc) - 1)
    ratios = eigs_desc[:km] / eigs_desc[1 : km + 1]
    return int(np.argmax(ratios)) + 1


def rank_selector(
    donors_pre: np.ndarray, sigma: float, c: float, k_max: int = 4
) -> int:
    """Gap-ratio selector with a TW-style silence gate at the MP edge.

    Returns 0 while the top normalized eigenvalue sits inside the noise band
    (lambda_hat <= 1.05 * sigma^2 (1+sqrt(c))^2, witness-P5 rule); otherwise
    returns the largest-gap position up to k_max. Proper TW calibration is a
    Phase C deliverable (WP-C3).
    """
    evals = np.sort(np.linalg.eigvalsh(unit_scatter(donors_pre)))[::-1]
    edge = sigma**2 * (1.0 + np.sqrt(c)) ** 2
    if evals[0] <= 1.05 * edge:
        return 0
    return select_rank_gap_ratio(evals, k_max)


def spectral_sc(
    donors: np.ndarray,
    y1_pre: np.ndarray,
    k_max: int = 4,
    sigma: float | None = None,
    c: float | None = None,
    info: dict | None = None,
) -> np.ndarray:
    """Hard-threshold spectral SC (frontier_ansatz.md Section 1 form).

    beta_j = <y1, v_j>/d_j on top-k PC scores; post predictions transport the
    donor cross-sections through the sample left basis. Rank by gap ratio;
    when sigma and c are supplied, a silence gate may force k = 0, in which
    case the forecast falls back to the constant mean(y1_pre).
    """
    T0 = y1_pre.shape[0]
    Xd = donors[:, :T0]
    U, d, Vt = np.linalg.svd(Xd, full_matrices=False)
    if sigma is not None and c is not None:
        eigs_unit = d**2 / T0
        edge = sigma**2 * (1.0 + np.sqrt(c)) ** 2
        k = 0 if eigs_unit[0] <= 1.05 * edge else select_rank_gap_ratio(eigs_unit, k_max)
    else:
        k = select_rank_gap_ratio(d**2, k_max)
    if info is not None:
        info["spectral_k"] = int(k)
    if k == 0:
        return np.full(donors.shape[1] - T0, y1_pre.mean())
    beta = (Vt[:k] @ y1_pre) / d[:k]
    sh = U[:, :k].T @ donors[:, T0:]
    return beta @ sh


def spectral_sc_full(
    donors: np.ndarray,
    y1_pre: np.ndarray,
    sigma: float | None = None,
    c: float | None = None,
    k_max: int = 4,
    z: float = 1.959963985,
) -> dict:
    """Gated + ungated spectral SC from ONE SVD, with the preregistered 95%
    prediction interval for the gated fit (preregistration Section 2.4).

    Returns dict with pred_gated, pred_ungated, k_gated, k_ungated, ci_lo,
    ci_hi and info. When sigma/c are absent the gate is off and both
    variants coincide with the plain gap-ratio rule.
    """
    T0 = y1_pre.shape[0]
    Xd, Xp = donors[:, :T0], donors[:, T0:]
    U, d, Vt = np.linalg.svd(Xd, full_matrices=False)
    evals_unit = d**2 / T0
    k_u = select_rank_gap_ratio(evals_unit, k_max)
    k_g = k_u
    if sigma is not None and c is not None:
        edge = sigma**2 * (1.0 + np.sqrt(c)) ** 2
        k_g = 0 if evals_unit[0] <= 1.05 * edge else k_u
    Tp = Xp.shape[1]

    def _predict(k):
        if k == 0:
            ybar = float(y1_pre.mean())
            s2 = float(np.sum((y1_pre - ybar) ** 2) / T0)
            half = z * np.sqrt(s2 * (1.0 + 1.0 / T0))
            return np.full(Tp, ybar), np.full(Tp, ybar - half), np.full(Tp, ybar + half), s2
        beta = (Vt[:k] @ y1_pre) / d[:k]
        sh = U[:, :k].T @ Xp                      # k x Tp post scores
        pred = beta @ sh
        fitted = Vt[:k].T @ (beta * d[:k])        # in-sample pre fit (time space)
        rss = float(np.sum((y1_pre - fitted) ** 2))
        s2 = rss / max(T0 - k, 1)
        h = np.sum(sh**2 / d[:k, None] ** 2, axis=0)
        half = z * np.sqrt(s2 * (1.0 + h))
        return pred, pred - half, pred + half, s2

    pg, lo, hi, _ = _predict(k_g)
    pu, _, _, _ = _predict(k_u)
    return {
        "pred_gated": pg,
        "pred_ungated": pu,
        "k_gated": int(k_g),
        "k_ungated": int(k_u),
        "ci_lo": lo,
        "ci_hi": hi,
        "info": {"lambda1": float(evals_unit[0]), "edge": float(sigma**2 * (1.0 + np.sqrt(c)) ** 2) if sigma is not None and c is not None else None},
    }


def _soft_impute(M_obs, mask, lam, iters, tol=1e-4, X0=None):
    X = np.where(mask, M_obs, 0.0) if X0 is None else X0.copy()
    denom = np.linalg.norm(np.where(mask, M_obs, 0.0)) + 1e-12
    for _ in range(iters):
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        s_st = np.maximum(s - lam, 0.0)
        X_new = (U * s_st) @ Vt
        X_new[mask] = M_obs[mask]
        delta = np.linalg.norm(X_new - X) / denom
        X = X_new
        if delta < tol:
            break
    return X


def mc_nn_cv(
    donors: np.ndarray,
    y1_pre: np.ndarray,
    lambdas=DEFAULT_MC_LAMBDAS,
    cv_holdout: float = 0.1,
    cv_seed: int = 4321,
    iters: int = 60,
    info: dict | None = None,
) -> np.ndarray:
    """Nuclear-norm MC (soft-impute) with lambda by held-out CV.

    Observed panel stacks the treated pre-row over the donors; only the
    treated post block is missing. CV masks a random fraction of observed
    entries, selects lambda by held-out MSE (warm-started down the grid),
    refits on all observed entries, reads out the completed treated-post block
    (MC-NNM-style block completion, Athey et al. 2021 setting without unit or
    time fixed-effects terms).
    """
    T0 = y1_pre.shape[0]
    T = donors.shape[1]
    M = np.vstack([np.concatenate([y1_pre, np.full(T - T0, np.nan)]), donors])
    obs = ~np.isnan(M)
    rng = np.random.default_rng(cv_seed)
    cand = np.argwhere(obs)
    n_hold = max(1, int(cv_holdout * len(cand)))
    hold_idx = cand[rng.choice(len(cand), size=n_hold, replace=False)]
    cv_mask = obs.copy()
    cv_mask[hold_idx[:, 0], hold_idx[:, 1]] = False
    M0 = np.where(obs, M, 0.0)
    scale = float(np.std(M[obs]))
    best_lam, best_err = None, np.inf
    Xw = None
    for lam in sorted(lambdas, reverse=True):
        Xi = _soft_impute(M0, cv_mask, lam * scale, iters, X0=Xw)
        Xw = Xi
        err = float(
            np.sum((Xi[hold_idx[:, 0], hold_idx[:, 1]] - M[hold_idx[:, 0], hold_idx[:, 1]]) ** 2)
        )
        if err < best_err:
            best_err, best_lam = err, lam
    Xhat = _soft_impute(M0, obs, best_lam * scale, iters)
    if info is not None:
        info["mc_lambda"] = float(best_lam)
    return Xhat[0, T0:]


def _simplex_intercept_ridge(Z, target, row_w, reg):
    """min over p in simplex(|coef|) and intercept a of
    sum_i row_w_i * (Z[i] @ coef + a - target[i])^2 + reg * ||coef||^2."""
    def obj(p):
        coef, a = p[:-1], p[-1]
        r = Z @ coef + a - target
        return float(row_w @ (r * r)) + reg * float(coef @ coef)

    def jac(p):
        coef, a = p[:-1], p[-1]
        r = Z @ coef + a - target
        gw = 2.0 * (Z.T @ (row_w * r)) + 2.0 * reg * coef
        ga = np.array([2.0 * float(row_w @ r)])
        return np.concatenate([gw, ga])

    m = Z.shape[1]
    p0 = np.r_[np.full(m, 1.0 / m), float(target.mean())]
    res = minimize(
        obj, p0, jac=jac, method="SLSQP",
        bounds=[(0.0, 1.0)] * m + [(None, None)],
        constraints=[
            {"type": "eq", "fun": lambda p: p[:-1].sum() - 1.0,
             "jac": lambda p: np.r_[np.ones(m), 0.0]}
        ],
        options={"maxiter": 600, "ftol": 1e-12},
    )
    ok = bool(res.success) or (
        res.status == 8 and abs(res.x[:-1].sum() - 1.0) < 1e-6
    )
    return res.x[:-1], res.x[-1], ok


def sdid(
    donors: np.ndarray,
    y1_pre: np.ndarray,
    reg_scale: float = 0.25,
    sweeps: int = 1,
    info: dict | None = None,
) -> np.ndarray:
    """    Reduced SDID port (Arkhangelsky et al. 2021, Algorithm 3 shape).

    Step 1 fits time weights beta (simplex + intercept) aligning pre-period
    donor levels with each donor's post average; step 2 fits unit weights
    (simplex + intercept) on beta^2-weighted pre periods; predictions are
    intercept + weighted donor outcomes. Weights enter through ROW WEIGHTS on
    residuals (never through design scaling, which would silently rescale the
    ridge penalty relative to the data term). Ridge regularization uses a
    variance heuristic (reg_scale); two alternating sweeps; no inference.
    """
    T0 = y1_pre.shape[0]
    Ypre = donors[:, :T0]
    Ypost = donors[:, T0:]
    n_d = Ypre.shape[0]
    var_y = float(np.var(np.diff(Ypre, axis=1))) + 1e-8

    beta = np.full(T0, 1.0 / T0)
    ok_all = True
    w = np.full(n_d, 1.0 / n_d)
    a = float(y1_pre.mean() - w @ Ypre.mean(axis=1))
    for _ in range(sweeps):
        beta_t, b_t, ok1 = _simplex_intercept_ridge(
            Ypre, Ypost.mean(axis=1), np.ones(n_d), reg_scale * var_y / n_d
        )
        beta = np.clip(beta_t, 1e-9, None)
        beta /= beta.sum()
        ok_all &= ok1
        w_t, a_t, ok2 = _simplex_intercept_ridge(
            Ypre.T, y1_pre, beta**2, reg_scale * var_y * float(beta @ beta)
        )
        w, a = np.clip(w_t, 0.0, None), float(a_t)
        ok_all &= ok2
    if info is not None:
        info["sdid_solver_ok"] = bool(ok_all)
    return w @ Ypost + a
