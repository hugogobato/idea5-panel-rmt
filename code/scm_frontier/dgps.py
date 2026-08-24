"""Factor-panel data-generating processes for the Idea 5 project.

Model (model_card.md): Y_it = L_it + E_it, L_it = sum_k A_ik f_kt with
standardized factor scores, spike strengths s_j = ||a^(j)||^2 / sigma^2,
aspect ratio c = n / T0. Unit 0 is treated; units 1..n-1 are donors.
All spectral diagnostics must use pre-period columns only.
"""

from __future__ import annotations

import numpy as np

NOISE_LAWS = ("gaussian", "ar1", "heteroskedastic")
ALIGNMENTS = ("none", "first", "all")


def _draw_noise(rng, law, shape, sigma, rho, het_ratio):
    if law == "gaussian":
        return rng.normal(0.0, sigma, size=shape)
    if law == "ar1":
        e = np.empty(shape)
        innov_sd = sigma * np.sqrt(1.0 - rho**2)
        e[:, 0] = rng.normal(0.0, sigma, size=shape[0])
        for t in range(1, shape[1]):
            e[:, t] = rho * e[:, t - 1] + rng.normal(0.0, innov_sd, size=shape[0])
        return e
    if law == "heteroskedastic":
        base = rng.normal(0.0, sigma, size=shape)
        hi = rng.random(shape[0]) < 0.5
        scale = np.where(hi, np.sqrt(het_ratio), 1.0 / np.sqrt(het_ratio))
        return base * scale[:, None]
    raise ValueError(f"unknown noise law {law!r}; expected one of {NOISE_LAWS}")


def generate_panel(
    n: int = 121,
    T0: int = 240,
    T_post: int = 100,
    r: int = 1,
    spike_strengths=(6.0,),
    treated_share=(0.5,),
    alignment: str = "first",
    sigma: float = 1.0,
    noise: str = "gaussian",
    rho: float = 0.5,
    het_ratio: float = 4.0,
    structural_break: float | None = None,
    seed: int = 0,
) -> dict:
    """Generate one synthetic panel.

    Parameters
    ----------
    n : total units including the treated unit at row 0.
    T0, T_post : pre- and post-period counts; c = n / T0.
    r : number of factors.
    spike_strengths : donor-carried s_j = ||a^(j)||^2 / sigma^2 per factor.
    treated_share : theta_j = alpha_j^2 / sigma^2 per factor; combined with
        `alignment`: "none" zeroes all alpha, "first" puts leverage only on
        factor 0, "all" spreads it equally across factors.
    structural_break : if not None, post-period factor scores are shifted by
        this amount per unit variance (violates assumption A4).
    seed : numpy default_rng seed; identical seeds give bitwise-identical panels.

    Returns
    -------
    dict with keys Y (n x (T0+T_post)), L, E, F ((T0+T_post) x r), A (n x r),
    and the config echo.
    """
    if noise not in NOISE_LAWS:
        raise ValueError(f"unknown noise law {noise!r}")
    if alignment not in ALIGNMENTS:
        raise ValueError(f"unknown alignment {alignment!r}; expected {ALIGNMENTS}")
    if len(spike_strengths) != r or len(treated_share) != r:
        raise ValueError("spike_strengths and treated_share must have length r")
    rng = np.random.default_rng(seed)
    n_d = n - 1
    T = T0 + T_post

    A = np.zeros((n, r))
    for j, sj in enumerate(spike_strengths):
        A[1:, j] = rng.normal(0.0, np.sqrt(sj * sigma**2 / n_d), size=n_d)
    if alignment == "none":
        pass
    elif alignment == "first":
        A[0, :] = 0.0
        A[0, 0] = sigma * np.sqrt(treated_share[0])
    else:
        A[0, :] = [sigma * np.sqrt(th / r) for th in treated_share]

    F = rng.normal(0.0, 1.0, size=(T, r))
    if structural_break is not None:
        F[T0:, :] += structural_break

    L = A @ F.T
    E = _draw_noise(rng, noise, (n, T), sigma, rho, het_ratio)
    Y = L + E
    config = dict(
        n=n, T0=T0, T_post=T_post, c=n / T0, r=r,
        spike_strengths=tuple(spike_strengths), treated_share=tuple(treated_share),
        alignment=alignment, sigma=sigma, noise=noise, rho=rho,
        het_ratio=het_ratio, structural_break=structural_break, seed=seed,
    )
    return {"Y": Y, "L": L, "E": E, "F": F, "A": A, "config": config}
