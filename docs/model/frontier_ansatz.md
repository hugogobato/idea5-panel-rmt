# Frontier Ansatz (WP-B1): Conjectured Recoverability Frontier for Spectral SC

**Status:** CONJECTURE. Nothing in this document is proved. It states the deterministic-equivalent (DE) ansatz for the counterfactual risk of hard-threshold spectral SC, verifies its special-case algebra symbolically (`code/check_frontier_ansatz.py`), and maps every ingredient to its source result. Rigorous derivation is Phase E target T1 (ridge variant) and T2 (lower bound). Per the plan, a numerical overlay of this formula on simulations is the WP-C1 falsifier, not a Phase B task; the symbolic checks here are internal consistency only.
**Date:** 2026-08-24. **Conventions:** inherit `model_card.md` (spiked covariance calibration, sigma^2 units, c = n/T0).

---

## 1. Estimator class covered (hard-threshold spectral SC)

Estimator (the WP-B2 implementation matches this form exactly):

1. Compute the SVD of the donor pre-period matrix `Y_D` (rows = donors, columns = pre periods); let `d_j >= u_j v_j^T` be sample singular triplets, `j = 1..min(n_d, T0)`.
2. Select the top `k` directions (rank rule: successive eigenvalue-gap ratio, `k <= k_max`; the ansatz treats `k` as given).
3. Regress the treated pre-period row `r = y_1,pre` on the score matrix `S = Y_D V_k` (OLS; equivalently `beta_j = <r, v_j> / d_j`).
4. For each post period `t`, form post scores from the donor cross-section: `sh_jt = u_j^T y_{.,t}^{post}` (pre-period basis transported to post windows; factor persistence, Assumptions A3-A4).
5. Predict `yhat*_t = sum_{j<=k} beta_j sh_jt`.

All risk statements below are for the *realized-outcome* estimand of `model_card.md` Section 4: target `y*_t = L_1t + E_1t`, so every method carries the irreducible floor `rho >= 1` in normalized MSE^2 units (`MSE/sigma^2`). This refines the plan's phrase "risk -> interpolation floor": in our estimand the limiting floor is the noise floor `sigma`, reached because the signal-component risk vanishes.

## 2. Single-spike risk decomposition (r = 1, k = 1)

Write the treated loading as `alpha` (so `L_1t = alpha f_t`, `E f_t^2 = 1`) and the donor-carried spike strength as `s = ||a_D||^2 / sigma^2`. Three population objects drive everything (sources in Section 5):

```
lambda(s, c) = 1 + s + c + c/s          [BBP/BGN outlier location, s > sqrt(c)]
zeta(s, c)   = (1 - c/s^2) / (1 + c/s)  [BGN squared eigenvector overlap]
tau(s, c)    = sqrt(s / lambda)         [singular-value transmission ratio]
```

Expanding the prediction error into (signal bias, coefficient-noise-times-post-signal, post-noise-times-coefficient-signal, noise-times-noise) and taking expectations (details of the expansion are in Section 6), the normalized risk converges to:

```
rho_excess(s, alpha; c) = (alpha^2/sigma^2) * [ (1 - zeta*tau)^2 + zeta/lambda ]
                        + (s*zeta + 1) / (T0 * lambda)
```

with `rho = 1 + rho_excess` and `o(1)` terms suppressed. Interpretation of the four channels:

1. `(alpha^2/sigma^2)(1 - zeta tau)^2`: systematic under-transmission. The pipeline transmits the treated signal component with gain `zeta*tau < 1` (both the sample left/right directions de-align and the singular value inflates). Vanishes iff `zeta -> 1`.
2. `(alpha^2/sigma^2)(zeta/lambda)`: coefficient-signal times post-noise. Vanishes as the spike strengthens.
3. `s*zeta/(T0 lambda)`: coefficient-noise times post-signal (present even at `alpha = 0`, size `O(sigma^2/T0)`).
4. `1/(T0 lambda)`: noise-times-noise, always negligible.

## 3. Multi-spike ansatz (formula F)

For spikes `j = 1..r` with strengths `s_j`, treated loadings `alpha_j`, retained set `K = {j <= k}`:

```
F({s_j}, {alpha_j}, K; c, T0) :=
    1                                  [realized-noise floor]
  + sum_{j in K, s_j > sqrt(c)}  (alpha_j^2/sigma^2) * [ (1 - zeta_j tau_j)^2 + zeta_j/lambda_j ]
  + sum_{j in K, s_j > sqrt(c)}  (s_j*zeta_j + 1) / (T0 * lambda_j)
  + sum_{j not in K}             (alpha_j^2/sigma^2)        [truncation bias]
```

Subcritical spikes (`s_j <= sqrt(c)`) never appear inside `K` terms (no outlier to retain; if forced into the basis they act as noise directions, contributing only channel-3/4-size terms). Cross-spike leakage terms (`O(alpha_j alpha_l * <f_j, v_l>` overlaps) are neglected for well-separated spikes; near-degenerate spike pairs are a tagged limitation (Section 7).

Scaling caveat that matters for Phase C design: with `c` fixed and treated loading held at a *typical-donor* size (`alpha_j^2 ~ kappa * s_j * sigma^2 / n`), the signal-channel terms vanish as `n -> infinity` and every unit becomes recoverable. A non-vanishing frontier requires the treated unit to hold a non-vanishing *share* of the spike: parameterize `theta_j := alpha_j^2/sigma^2` directly (recommended for the WP-C1 grid; e.g. `theta in {0.1, 0.5, 1}`), or work at finite `n` where the `1/T0` channels still bite. Flag for `preregistration.md`.

## 4. Special-case reductions (verified symbolically)

Script: `code/check_frontier_ansatz.py` (sympy; run log in Section 8). All four PASS.

- SC1 (`r = 0` or `K` empty): every sum is empty, `F = 1` exactly. Matches Witness 2 NULL arm (all methods at `1.002-1.030 sigma`).
- SC2 (`s -> inf`, single spike): `zeta -> 1`, `tau -> 1`, `lambda/s -> 1`; hence `(1 - zeta tau)^2 -> 0`, `zeta/lambda -> 0`, `(s zeta + 1)/(T0 lambda) -> 1/T0 -> 0`; `F -> 1`. Signal risk hits the interpolation-free floor; only the realized-noise floor remains.
- SC3 (`alpha = 0` on all spikes): signal channels vanish identically; residual `(s zeta + 1)/(T0 lambda) -> 0` as `T0 -> infinity` at fixed `c`. Matches Witness 2 MISALIGNED arm (paired RMSE deltas `~ 10^-3 sigma`). Channel-2 logic holds exactly in the ansatz.
- SC4 (continuity at the edge): `lim_{s downarrow sqrt(c)} zeta = 0` implies the included-spike terms reduce to `(alpha^2/sigma^2)` (full under-transmission bias, zero variance benefit), which equals the truncation term for an excluded spike. The frontier is therefore CONTINUOUS at `m = 1` with a KINK, not a jump: retaining a just-supercritical direction costs nothing relative to exclusion, and helps only through the derivative. This is the precise sense in which the ansatz predicts a phase transition in *marginal value*, and why Witness 1's coefficient-detectability onset (which turns on through the TW fluctuation band below the edge) sits lower than the risk kink location.

Monotonicity (checked numerically in the script over `s in (sqrt(c), 20]`, five `c` values, four treated shares `theta`): the TOTAL excess risk is strictly decreasing in `s` at every tested configuration, and its global maximum over the supercritical side sits exactly AT the edge `s = sqrt(c)` where it connects continuously to the truncation value `(alpha^2/sigma^2)` of an excluded spike. The preregistered kink signature therefore holds exactly in the ansatz: flat below the edge (up to negligible subcritical loading terms), continuous peak at `m = 1`, steepest descent just above it (script-measured |slope| ratio `~ 2200x` versus the far-above-edge slope at `c = 0.5`). Secondary structure worth recording: the two variance channels individually are non-monotone just above the edge (`zeta/lambda` rises from 0 before decaying; sup between `0.05` and `0.26` as `c` goes from 4 down to 0.25), but their contribution is `O(1/T0)`-suppressed and never overturns the bias-channel decline in the total.

## 5. Ingredient-to-source map (also mirrored in theory_targets.md)

| Ansatz ingredient | Source result | Register ID | Adaptation gap |
|---|---|---|---|
| Bulk edge `sigma^2(1+sqrt(c))^2`; existence threshold `s > sqrt(c)` | Baik-Ben Arous-Peché (2005), real spiked covariance | T1 (register C) | None for existence; our use is standard |
| Outlier location `lambda(s,c)` | Benaych-Georges-Nadakuditi (2011), rectangular low-rank deformations | T2 | Direct |
| Overlap `zeta(s,c)` | BGN (2011), eigenvector overlap formulas | T2 | Direct (single-spike); multi-spike perturbation folklore |
| Finite-n fluctuation band around the edge (Witness 1 onset `< m=1`) | Johnstone (2001) TW law; Onatski (2010) ratio tests | T3, T5a | Band width under weak dependence = C3's open problem |
| Ridge-family companion risks (WP-C2 engine, not used in F) | Dobriban-Wager (2018); Hastie et al. (2022) | T6, T7 | Row-targeted (single unit) risk instead of Frobenius; deferred to T1 |
| Lower-bound side ("no estimator below floor when all spikes subcritical or alpha = 0") | NOT COVERED BY THIS ANSATZ | T2 target | Le Cam two-point template; Phase E |

Positioning obligation (watch item W-2, canonical text frozen in `evidence_register.md` S07): the frontier varies signal strength/alignment at fixed estimator class; Spiess et al.'s double/single-descent varies complexity (donor count) at fixed signal via model averaging. Orthogonal axes; joint test = descent curves flattening at the `sigma` floor exactly when the diagnostic reports sub-frontier distance.

## 6. Derivation sketch (what Phase E T1 must make rigorous)

With `v_hat` the sample top right singular vector and `u_hat` its left partner, BGN gives `cos^2(theta_v) -> zeta`, `cos^2(theta_u) -> zeta`, `d_hat/sigma -> sqrt(T0 * lambda)` (all in probability). Substituting,

```
beta      = <r, v_hat>/d_hat,       r = alpha f + eps_1
sh_t      = u_hat^T y_.t^post = sigma sqrt(s) cos(theta_u) f_t + eta_t,
Var(eps_1^T v_hat) = sigma^2,  Var(eta_t) = sigma^2,  eps_1 independent of post noise.
```

Prediction error `yhat*_t - alpha f_t - E_1t = alpha f_t (G - 1) + beta_noise sh_signal_t + beta_signal eta_t + beta_noise eta_t` with `G = zeta tau` (leading order; the mixed term `eps_1`-through-`v_hat`-rotation is second order because the coefficient noise enters already divided by `d_hat`). Taking variances/squares and dividing by `sigma^2` yields the two-line formula in Section 2. The `o(1)` corrections are `O(T0^{-1/2})` fluctuation terms and `O(zeta (1-zeta))` rotation-mixing terms; both are part of the T1 work package, not of this ansatz.

Numerical calibration point (non-decisive, recorded for Phase C overlay tolerance): at the Witness 2 ALIGNED cell (`s = 6, c = 0.5, T0 = 240, alpha^2/sigma^2 = 9*6/119 = 0.4538`), F predicts `rho_excess = 0.0751` (`bias 0.0167 + zeta/lambda channel 0.0548 + T0-channels 0.0036`), i.e. `RMSE = 1.0369 sigma` versus the measured `1.0395 sigma` (delta `0.25%`). Agreement at this cell is encouraging but is one point, not evidence.

## 7. Stated limitations (tagged on the record)

1. iid Gaussian noise only (A2); AR(1)/heteroskedastic drift is WP-C2(iii)'s question.
2. Well-separated spikes; near-degenerate pairs invalidate the per-direction decomposition.
3. Post-window factor transport assumed consistent on the sample basis (A3/A4); finite `T_post` adds a variance channel measured in Witness 2's battery and folded into Phase C cells, not into F.
4. Rank selection treated as oracle-given `K`; selector errors belong to C3 metrics.
5. The formula is a conjecture with verified algebra and one calibration point; its falsifier is the WP-C1 kink criterion (location within 15% of `m = 1` on >= 80% of the c-grid), plus the inherited onset-convergence check as `n` grows.

## 8. Symbolic verification log

Produced by `code/check_frontier_ansatz.py` (sympy 1.14.0 + numpy grid scans), final run 2026-08-24:

```
PASS  SC2 s->inf: bias->0
PASS  SC2 s->inf: zeta/lambda->0
PASS  SC2 s->inf: (s*zeta+1)/lambda->1 (so /T0 -> 0)
PASS  SC4 edge: zeta->0
PASS  SC4 edge: included-signal-bias -> 1 (= truncation value)
PASS  SC4 edge: zeta/lambda -> 0
PASS  edge tau^2 identity
PASS  bias channel strictly decreasing on supercritical grid
PASS  TOTAL excess risk strictly decreasing (all c, theta)
PASS  steepest descent near edge (kink signature)
|slope| near edge 0.4310 vs far 0.002335
ALL PASS
```

SC1/SC3 are structural (empty sums, vanishing alpha) and hold exactly in formula F. Development note for the record: two intermediate monotonicity claims drafted during the derivation were corrected before freeze because a test-harness scaling slip dropped the `1/T0` factor on channel 3; with the formula as specified, total risk is strictly decreasing above the edge at every tested configuration and its maximum sits exactly at the boundary. The frozen document states only verified claims.
