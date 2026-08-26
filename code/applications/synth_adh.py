"""Classic Synth-style SCM machinery (WP-D1 trusted-benchmark pipeline).

Replicates the published Abadie-Diamond-Hainmueller estimators closely
enough to reproduce their reported point estimates to reporting precision:

  * predictor matrices X0/X1 scaled row-wise by the combined donor+treated
    sample SD (Synth's internal convention);
  * unit weights w solve min (X1 - X0'w)' V (X1 - X0'w) s.t. w >= 0,
    sum(w) = 1 (SLSQP with analytic gradients);
  * predictor weights V solve min_v MSPE(Z1 - Z0' w(v)) over the
    preregistered SSR window, searched by seeded multistart Nelder-Mead on
    a softmax parameterization plus structured starts (uniform, one-hot),
    mirroring Synth's multi-method search without its ipop dependency.

All stochasticity flows through seeds registered in seeds.yaml
(phase_d.wp_d1.v_multistart_seed). Also provides the Arkhangelsky et al.
(2021) SC and SDID point estimators exactly as defined in arXiv:1812.09970v4
(Eqs. 2.1-2.4, Algorithm 1), which serve as the second reproduction
instrument on the outcome-only Prop 99 panel.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize, nnls


def scale_rows(X0: np.ndarray, X1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Synth scaling: divide each predictor row by its combined SD."""
    x1 = np.asarray(X1, float)
    if x1.ndim == 1:
        x1 = x1[:, None]
    big = np.hstack([np.asarray(X0, float), x1])
    sd = big.std(axis=1, ddof=1)
    if not (sd > 0).all():
        raise ValueError("zero-variance predictor row")
    return (np.asarray(X0, float) / sd[:, None],
            (x1 / sd[:, None]).ravel())


def _simplex_qp(V: np.ndarray, X0s: np.ndarray, X1s: np.ndarray,
                lam_aug: float | None = None):
    """min (X0'w-X1)' V (X0'w-X1) over w in simplex; returns (w, ok).

    Convention: X0s is (predictors x units), X1s is (predictors,). Solved
    exactly by Lawson-Hanson NNLS on the augmented system that encodes the
    equality constraint with a large penalty row (SLSQP proved unreliable
    on this problem class: it terminates at the initial point reporting
    success; documented in preprocessing_frozen.md Section 5).
    """
    v = np.asarray(V, float)
    if v.ndim == 2:
        v = np.diag(v)
    root = np.sqrt(v)
    M = root[:, None] * np.asarray(X0s, float)   # (m x n_units) design
    t = root * np.asarray(X1s, float)
    n = X0s.shape[1]
    scale = float(np.abs(M).max()) + 1e-12
    lam = scale * 1e6 if lam_aug is None else lam_aug
    A = np.vstack([M, lam * np.ones((1, n))])
    b = np.concatenate([t, [lam]])
    w, _ = nnls(A, b)
    s = w.sum()
    ok = s > 1e-12
    if ok:
        w = w / s
    return w, bool(ok)


def _mspe_of(v: np.ndarray, X0s: np.ndarray, X1s: np.ndarray,
             Z0: np.ndarray, Z1: np.ndarray) -> float:
    w, _ = _simplex_qp(np.diag(v), X0s, X1s)
    resid = Z1 - Z0.T @ w
    return float(resid @ resid / Z0.shape[1])


def fit_v(X0: np.ndarray, X1: np.ndarray, Z0: np.ndarray, Z1: np.ndarray,
          seed: int = 60_001, n_random: int = 16, maxiter: int = 250,
          fatol: float = 1e-9, polish_rounds: int = 1) -> dict:
    """Seeded multistart search for V minimizing outcome-path MSPE."""
    X0s, X1s = scale_rows(X0, X1)
    m = X0.shape[0]
    rng = np.random.default_rng(seed)

    def expand(u):
        e = np.exp(u - u.max())
        return e / e.sum()

    starts = [np.full(m, 1.0 / m)]
    starts += [np.eye(m)[j] for j in range(m)]
    starts += list(rng.dirichlet(np.ones(m), size=n_random))

    best = (np.inf, None, None)
    for v0 in starts:
        res = minimize(
            lambda u: _mspe_of(expand(u), X0s, X1s, Z0, Z1),
            np.log(np.maximum(v0, 1e-8)), method="Nelder-Mead",
            options={"maxiter": maxiter, "xatol": 1e-6, "fatol": fatol},
        )
        v_hat = expand(res.x)
        mspe = _mspe_of(v_hat, X0s, X1s, Z0, Z1)
        if mspe < best[0]:
            w_hat, _ = _simplex_qp(np.diag(v_hat), X0s, X1s)
            best = (mspe, v_hat, w_hat)
    mspe, v_hat, _ = best
    # Deterministic polish: restart NM at the incumbent optimum.
    u_hat = np.log(np.maximum(v_hat, 1e-8))
    for _ in range(polish_rounds):
        res = minimize(
            lambda u: _mspe_of(expand(u), X0s, X1s, Z0, Z1), u_hat,
            method="Nelder-Mead",
            options={"maxiter": maxiter, "xatol": 1e-8,
                     "fatol": fatol * 0.1},
        )
        v_pol = expand(res.x)
        mspe_pol = _mspe_of(v_pol, X0s, X1s, Z0, Z1)
        if mspe_pol < mspe:
            v_hat, mspe, u_hat = v_pol, mspe_pol, res.x
    w_hat, _ = _simplex_qp(np.diag(v_hat), X0s, X1s)
    return {"v": v_hat, "w": w_hat, "mspe": mspe,
            "scaled": (X0s, X1s)}


def fit_w_given_v(V: np.ndarray, X0: np.ndarray, X1: np.ndarray) -> dict:
    """Stage-2 fit: re-scale predictors, solve w at the frozen V."""
    X0s, X1s = scale_rows(X0, X1)
    Vm = np.diag(np.asarray(V, float)) if np.asarray(V).ndim == 1 \
        else np.asarray(V, float)
    w, ok = _simplex_qp(Vm, X0s, X1s)
    return {"w": w, "solver_ok": ok, "scaled": (X0s, X1s)}


def split_X(X_all: np.ndarray, treated_idx: int):
    """(m x units) predictor block -> (donor block m x n_d, treated vector)."""
    mask = np.arange(X_all.shape[0]) != treated_idx
    return X_all[mask].T, X_all[treated_idx]


def synth_adh_germany(X_train: np.ndarray, treated_idx: int,
                      Z0_train: np.ndarray, Z1_train: np.ndarray,
                      X_main: np.ndarray, seed: int = 60_001) -> dict:
    """ADH 2015 two-stage procedure (rep_updated.r):

    stage 1 picks V by fitting the training-window outcome path; stage 2
    re-solves w on the main predictor windows with that V frozen.
    """
    X0_tr, x1_tr = split_X(X_train, treated_idx)
    st1 = fit_v(X0_tr, x1_tr, Z0_train, Z1_train, seed=seed)
    X0_ma, x1_ma = split_X(X_main, treated_idx)
    st2 = fit_w_given_v(st1["v"], X0_ma, x1_ma)
    return {"v": st1["v"], "w": st2["w"], "stage1_mspe": st1["mspe"],
            "stage2_solver_ok": st2["solver_ok"]}


def synth_adh_smoking(X: np.ndarray, Y: np.ndarray, T0: int,
                      seed: int = 60_001) -> dict:
    """ADH 2010 single-stage procedure: V minimizes the pre-treatment
    outcome MSPE over the FULL pre window (paper Section 2.3)."""
    return fit_v(X[:-1].T, X[-1], Y[:-1, :T0], Y[-1, :T0], seed=seed,
                 n_random=16, maxiter=250, polish_rounds=1)


# --------------------- Arkhangelsky et al. (2021) instruments --------------


def _sdid_sigma_hat(Ydonors_pre: np.ndarray) -> float:
    d = np.diff(Ydonors_pre, axis=1)
    return float(d.std(ddof=1))


def _fw_step(A: np.ndarray, x: np.ndarray, b: np.ndarray,
             eta: float) -> np.ndarray:
    """One exact-line-search Frank-Wolfe step (synthdid R fw.step)."""
    Ax = A @ x
    half_grad = A.T @ (Ax - b) + eta * x
    i = int(np.argmin(half_grad))
    d_x = -x.copy()
    d_x[i] += 1.0
    if np.allclose(d_x, 0.0):
        return x
    d_err = A[:, i] - Ax
    step = -(half_grad @ d_x) / (float(d_err @ d_err) + eta * float(d_x @ d_x))
    step = min(1.0, max(0.0, step))
    return x + step * d_x


def _sparsify(v: np.ndarray) -> np.ndarray:
    """synthdid R sparsify_function: zero weights <= max/4, renormalize."""
    v = v.copy()
    v[v <= v.max() / 4.0] = 0.0
    s = v.sum()
    return v / s if s > 0 else v


def _sc_weight_fw(A: np.ndarray, b: np.ndarray, zeta_sq_eta: float,
                  intercept: bool, noise_level: float,
                  max_iter: int = 10_000,
                  max_iter_pre_sparsify: int = 100,
                  min_decrease_mult: float = 1e-5) -> np.ndarray:
    """Port of synthdid R sc.weight.fw incl. the two-round sparsify.

    A is (obs x k), b the target; minimizes ||Ax-b||^2/n_obs +
    zeta^2 ||x||^2 over the simplex, zeta^2*n_obs = zeta_sq_eta. Starts at
    uniform, runs Frank-Wolfe to max_iter_pre_sparsify, applies the
    sparsify tie-break, then continues to max_iter.
    """
    P = A.copy()
    y = b.copy()
    if intercept:
        P = P - P.mean(axis=0, keepdims=True)
        y = y - y.mean()

    def run(x0, max_it):
        x = x0.copy()
        k = len(x)
        vals = []
        t = 0
        while t < max_it and (t < 2 or vals[t - 2] - vals[t - 1]
                              > (min_decrease_mult * noise_level) ** 2):
            x = _fw_step(P, x, y, zeta_sq_eta)
            err = P @ x - y
            vals.append(float(err @ err) / P.shape[0]
                        + (zeta_sq_eta / P.shape[0]) * float(x @ x))
            t += 1
        return x

    k = P.shape[1]
    x = run(np.full(k, 1.0 / k), max_iter_pre_sparsify)
    x = run(_sparsify(x), max_iter)
    s = x.sum()
    return x / s if s > 0 else np.full(k, 1.0 / k)


def _nnls_simplex_ridge(P: np.ndarray, y: np.ndarray, ridge2: float = 0.0,
                        center: bool = True):
    """min_w sum_t (y_t - [a] - P_t w)^2 [+ ridge2*||w||^2] over w >= 0.

    With center=True the intercept a is profiled out by centering; the
    equality sum(w)=1 is enforced by renormalization. P is (n obs x k).
    Retained for cross-checks; the shipped synthdid instruments use the
    faithful Frank-Wolfe port above.
    """
    Pc = P - P.mean(axis=0, keepdims=True) if center else P
    tc = y - y.mean() if center else y
    k = P.shape[1]
    A = [Pc]
    b = [tc]
    if ridge2 > 0:
        A.append(np.sqrt(ridge2) * np.eye(k))
        b.append(np.zeros(k))
    w, _ = nnls(np.vstack(A), np.concatenate(b))
    s = w.sum()
    return w / s if s > 0 else np.full(k, 1.0 / k)


def synthdid_sc(Y: np.ndarray, T0: int, treated_idx: int) -> dict:
    """SC estimator per the synthdid reference implementation
    (sc_estimate: omega-intercept FALSE, zeta.omega = 1e-6 sigma-hat,
    lambda identically zero; estimate = post-period mean of the omega-gap
    in levels). Weights via the sparsified Frank-Wolfe solver."""
    tr = Y[treated_idx]
    co = np.delete(Y, treated_idx, axis=0)
    Ypre_co, Ypost_co = co[:, :T0], co[:, T0:]
    sig = _sdid_sigma_hat(Ypre_co)
    zeta = 1e-6 * sig
    omega = _sc_weight_fw(Ypre_co.T, tr[:T0], T0 * zeta**2,
                          intercept=False, noise_level=sig)
    tau = float((tr[T0:] - omega @ Ypost_co).mean())
    return {"tau": tau, "w": omega, "sigma_hat": sig}


def synthdid_did(Y: np.ndarray, T0: int, treated_idx: int) -> float:
    """DID comparator: uniform omega and lambda (did_estimate)."""
    tr = Y[treated_idx]
    co = np.delete(Y, treated_idx, axis=0)
    delta_tr = tr[T0:].mean() - tr[:T0].mean()
    delta_co = co[:, T0:].mean(axis=1) - co[:, :T0].mean(axis=1)
    return float(delta_tr - delta_co.mean())


def synthdid_att(Y: np.ndarray, T0: int, treated_idx: int) -> dict:
    """Full SDID point estimate mirroring synthdid_estimate exactly:

    lambda from sparsified FW on {donors x (pre | post-mean)} with
    intercept and zeta.lambda = 1e-6 sigma-hat; omega from sparsified FW
    on {time x (donors | treated)} with intercept and zeta.omega =
    (N_tr*T_post)^(1/4) sigma-hat; estimate =
    c(-omega, 1)' Y c(-lambda, 1/T_post).
    """
    n, T = Y.shape
    tr = Y[treated_idx]
    co = np.delete(Y, treated_idx, axis=0)
    n_co = co.shape[0]
    T1 = T - T0
    Ypre_co, Ypost_co = co[:, :T0], co[:, T0:]
    sig = _sdid_sigma_hat(Ypre_co)

    zeta_lambda = 1e-6 * sig
    lam = _sc_weight_fw(Ypre_co, Ypost_co.mean(axis=1),
                        n_co * zeta_lambda**2, intercept=True,
                        noise_level=sig)

    zeta_omega = ((1 * T1) ** 0.25) * sig
    omg = _sc_weight_fw(Ypre_co.T, tr[:T0], T0 * zeta_omega**2,
                        intercept=True, noise_level=sig)

    g_pre = tr[:T0] - Ypre_co.T @ omg
    g_post = tr[T0:] - Ypost_co.T @ omg
    tau = float(g_post.mean() - lam @ g_pre)
    return {"tau": tau, "omega": omg, "lambda": lam,
            "zeta_omega": zeta_omega, "sigma_hat": sig}


