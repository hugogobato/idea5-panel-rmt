# Model Card (WP-A1): Spectral Recoverability Frontiers for Panel Causal Inference

**Project:** Idea 5, phase transitions in panel causal inference.
**Work package:** WP-A1 of the research plan (`Idea5_Panel_RMT_Research_Plan.md`, Section 7, Phase A).
**Status:** Frozen before witness execution. The pass rules in Section 7 were written down here before either witness notebook was run.
**Date:** 2026-08-24.

---

## 1. Purpose and scope

This card fixes the formal model, the estimand, the calibration conventions, and the assumption ledger that all Phase A witnesses (and later Phase B formalization) must respect. It exists so that the two numerical witnesses test a *pre-specified* mechanism rather than one reverse-engineered from simulation output. No claim in this card is proved; every asymptotic statement is tagged CONJECTURE until WP-B1/WP-E1 derive it. The card also states the nonidentifiability region explicitly, per Gate G0.

## 2. Data regime and notation

Units `i = 1..n` are observed over periods `t = 1..T`. Unit `i = 1` is treated; units `i = 2..n` are donors (`n_d = n - 1`). Let `T0` be the number of pre-treatment periods and `T_post` the post window. Both `n` and `T0` grow with fixed aspect ratio

```
c = n / T0  in (0, inf).
```

The pre-period outcome matrix is `Y_pre in R^{n x T0}`; donor outcomes continue into the post window. All spectral diagnostics use only `Y_pre` (leakage rule from G0).

## 3. Data-generating model

Potential-outcome decomposition:

```
Y_it = L_it + tau_it Z_it + E_it
```

with latent trajectory `L` of rank `r` written in factor form

```
L_it = sum_{k=1}^r a_ik f_kt ,      A = [a^(1) ... a^(r)] in R^{n x r}
```

where the factor scores satisfy `E[f_kt] = 0`, `E[f_kt^2] = sigma_f,k^2 = 1` after standardization, and loadings are deterministic. Noise `E_it` is iid `N(0, sigma^2)` in Phase A (Assumption A2 below governs relaxations later). Treatment effects `tau_it Z_it` occupy the post period of unit 1 and play no role in Phase A.

Equivalent singular-value view: `Y_pre = U V^T + E` with `U = A F_pre^{1/2}-scaled`; both views are used interchangeably.

### 3.1 Spiked-covariance calibration (the conventions all code must use)

Row scatter of the full panel:

```
C_hat = (1/T0) Y_pre Y_pre^T   ->   Sigma = A A^T + sigma^2 I_n .
```

Define the **spike excess** of factor k:

```
s_k := ||a^(k)||^2 / sigma^2 .
```

With `gamma := c = n/T0`, the Marchenko-Pastur upper edge of the noise bulk (eigenvalue scale) and the BBP detectability condition are:

```
lambda_+ = sigma^2 (1 + sqrt(c))^2          [bulk edge]
s_k > sqrt(c)                              [supercritical: outlier exists]
s_k <= sqrt(c)                             [subcritical: no outlier]
```

For a supercritical spike the limiting top eigenvalue is CONJECTURE-tagged but standard (BBP/BGN):

```
E[eig_1(C_hat)] -> sigma^2 (1 + s_k + c + c / s_k)     when s_k > sqrt(c).
```

The **spike strength multiplier** used on all grids is

```
m := s / (sigma^0 * sqrt(c)) = s / sqrt(c),        critical value m* = 1.
```

The treated unit's leverage on factor k is `alpha_k := a_{1k}`. Two distinct failure channels exist and must never be conflated: subcritical spikes (channel 1, Witness 1) and zero leverage on supercritical spikes (channel 2, Witness 2).

## 4. Estimand

Primary target: the treated unit's untreated potential outcome path in the post window,

```
y*_post = { L_1t + E_1t : t in (T0, T0 + T_post] } ,
```

i.e., the *realized* control trajectory including its idiosyncratic component, matching what SCM practice compares against. The irreducible floor of any predictor that does not know `E_1` is therefore `RMSE = sigma`. ATT aggregation is deferred to later phases. Primary metric: post-period RMSE normalized by `sigma`.

## 5. Assumption ledger

| ID | Assumption | Status | Diagnostic plan |
|---|---|---|---|
| A1 | Additive linear factor structure `Y = AF^T + E`, r fixed as n, T0 grow | Maintained throughout | Scree/rank diagnostics |
| A2 | Noise iid `N(0, sigma^2)` across units and time | Maintained in Phase A; AR(1)/heteroskedastic relaxation scheduled Phase C (WP-C2 iii) | TW size plots under misspecification; block bootstrap fallback |
| A3 | Factor scores standardized, mean-zero, stationary across the cutoff | Maintained | Pre-period vs post-period moment checks |
| A4 | No structural break in V (factor score law) at the cutoff | Stated assumption, first-class | Spectral pre-trends test (C3); break DGP in WP-C2(iv) |
| A5 | Treated loading has nonzero projection onto the spiked subspace (`sum alpha_k^2 > 0`) | Stated assumption | Alignment estimate reported by diagnostic suite |
| A6 | Donors observe the same factors; no interference/spillovers | Maintained | Placebo-in-space batteries (Phase D) |
| A7 | Spike profile bounded away from the edge in either direction when limits quoted | Technical | Fine grids near m=1 quantify finite-size behavior |

Nonidentifiability region (must be surfaced, never hidden): if A4 fails, the post-period counterfactual is genuinely unidentifiable from pre-data; if A5 fails exactly (`alpha = 0` on all supercritical spikes), the target is unidentifiable regardless of spike strength; if all spikes are subcritical (channel 1), no estimator separates the treated row's factor component from noise. These three regions delimit what the paper may claim.

## 6. Conjectured recoverability frontier (informal statement)

CONJECTURE (to be formalized as function F in WP-B1). The normalized counterfactual risk of any panel method behaves as follows.

1. If all spikes are subcritical (`max_k s_k <= sqrt(c)`), every estimator attains risk within Monte Carlo indistinguishability of the pure-noise regression floor; the treated row carries no learnable signal.
2. If some spike is supercritical but the treated unit has zero leverage on it, the same floor holds despite a clearly visible outlier eigenvalue.
3. If a spike is supercritical and leverage is positive, spectral methods attain risk strictly below the floor, improving monotonically in `s_k` and in `alpha_k^2`, with a kink located at the BBP edge `m = 1`.

Witnesses W1 and W2 test observable consequences of regimes 1 to 3; neither proves the conjecture.

## 7. Witness preregistration (pass rules frozen before execution)

### 7.1 Witness 1, sub-threshold invisibility (`notebooks/witness_subthreshold.ipynb`)

Design: single factor, `n = 120`, `T0 = 240` (`c = 0.5`, bulk edge `(1+sqrt(0.5))^2 = 2.9143 sigma^2`). Donor loadings iid `N(0, s_d^2)` with `s_d^2 = s/(n-1)`; treated loading `alpha = 3.0 * s_d`. Grid over strength multiplier `m = s/sqrt(c) in linspace(0.55, 1.65, 23)`, `R = 300` replications per point. Seeds: one master scan stream seeded at 50101 (per grid point, in order: one loadings draw, then `R` replication draws); baseline pool of `R = 600` pure-noise replications seeded at 50102, split into six pools of 100 for the median KS statistic.

Statistic (amended before the definitive run, see deviation log): `q := (beta_hat)^2` with `beta_hat = < y_1 , v_hat_1 >` the first-PC coefficient. Rationale: the SVD sign is arbitrary per replication, so the signed statistic is symmetric about zero whenever signal is present and carries no information; the squared statistic tests the same hypothesis. For reference the null law of `beta_hat` at c = 0.5 has sd approximately `(1 + sqrt(c))/sqrt(c) sigma = 2.414 sigma` (dominated by the noise-fitting channel), so magnitudes, not signs, are informative.

Pass rules (amended once before the definitive run; original drafting and rationale for every amendment recorded in the deviation log of `gate_g0_g1_decision.md`):
1. P1 invisibility: median two-sample KS p-value (q statistic vs baseline pool) `>= 0.05` at every grid point with `m <= 0.70`.
2. P2 power: median KS p-value `<= 0.01` at every grid point with `m >= 1.30`.
3. P3 transition location: the smallest grid `m` with median KS p `< 0.01` lies in `[0.60, 1.45]`, and Spearman rank correlation between median p and `m` over the full grid is `< -0.7`.
4. P4 outlier location sanity: mean top eigenvalue of `(1/T0) Y Y^T` at `m = 1.65` within 15 percent of the predicted `sigma^2 (1 + s + c + c/s)`, `s = 1.65 sqrt(c)`.
5. P5 eigenvalue silence: mean top eigenvalue `<= 1.05 lambda_+` at every grid point with `m <= 0.95`.

Finite-size note (added after the diagnostic pass, before the definitive run): at `n = 120`, `T0 = 240` the coefficient-based detectability onset is expected somewhat BELOW the asymptotic edge, because the eigenvector overlap turns on continuously through the TW fluctuation band while the eigenvalue itself remains inside the bulk until further above threshold. The witness therefore tests localization up to the wide P3 window; sharp verification that the onset converges to `m = 1` as `n` grows is a Phase C falsifier (WP-C1), not a Phase A one.

Fail interpretation: P1/P2 failing symmetric to prediction (transition far from `m=1`) falsifies the ansatz location; P4 failing while P1-P3 pass indicates a calibration-formula error to be investigated, not a mechanism failure.

### 7.2 Witness 2, visible-but-useless spike (`notebooks/witness_misalignment.ipynb`)

Design: `n = 120`, `T0 = 240`, `T_post = 100`. One donor-carried factor with `s = 6` (`||a_donors||^2 = 6 sigma^2`, predicted outlier at `1 + 6 + 0.5 + 0.5/6 = 7.583 sigma^2`). Three arms sharing identical noise and factor draws per replication index (paired by construction, base seed `50201 + j`, `R = 400`):

1. NULL arm: all loadings zero (`r = 0`).
2. MISALIGNED arm: treated loading exactly `alpha = 0`; donors carry the factor.
3. ALIGNED comparator: same spike strength, treated loading `alpha = 3 s_d = 3 sqrt(6/119)`.

Estimators (weights fit on pre-periods only, predict `T_post` ahead): donor-mean (uniform weights); ridge-SC (donor features plus intercept, penalty by 4-fold contiguous-block CV); simplex SCM (Abadie weights via SLSQP, `w >= 0`, `sum w = 1`); hard-threshold spectral SC (regression on top-k PC scores of the donor matrix, `k` by largest successive eigenvalue-gap ratio, `k <= 4`). Metric: post-period RMSE against realized `y*_post`, normalized by `sigma`.

Pass rules (all required):
1. Q1 coexistence: MISALIGNED mean top donor-scatter eigenvalue `> 6.8` (outlier visibly present) while NULL arm mean `< 3.3` (at the finite-size edge).
2. Q2 uselessness: for every method, `|mean paired RMSE difference (MISALIGNED minus NULL)| < 1 sd of the per-replication paired difference` and the MISALIGNED mean does not exceed the NULL mean by more than 5 percent relative.
3. Q3 sensitivity control: because the target is the realized outcome (Section 4), the oracle-aligned floor is exactly `sigma` (a perfect factor predictor still leaves `E_1`); no estimator can go below it. The battery is therefore sensitive iff ALIGNED spectral-SC RMSE `<= 1.05 sigma` AND at least 10 percent below (paired) the ALIGNED donor-mean RMSE, which cannot extract the factor. A flat instrument would make Q2 vacuous.

Fail interpretation: Q2 failing because misaligned spikes HELP would contradict channel-2 logic and force a PIVOT; Q3 failing means the estimators are broken, not the theory.

## 8. What this card does not establish

No lower bound is proved; "no estimator" statements inside witnesses are operationalized as "none of the implemented estimators," which is the strongest claim simulations can support at Phase A. The deterministic-equivalent risk formula F, the matching impossibility theorem, and calibrated inference are WP-B1 and Phase E objects. Serially correlated or heteroskedastic noise is out of scope here (A2).

## 9. Reproducibility pointers

Seeds: see `seeds.yaml` (registry entries W1_SCAN_SEED = 50101, W1_NULL_SEED = 50102, W2_BASE_SEED = 50201; derivation formulas inside each notebook). Environment: Python 3.12, numpy 2.4.3, scipy 1.17.1, matplotlib 3.10.8. Execution: notebooks run top-to-bottom in a fresh kernel; figures land in `figures/`. Source anchors for all cited results live in `evidence_register.md`.

## 10. Key references (verification level E2-E3; anchors consolidated in evidence_register.md)

1. Baik, Ben Arous, Peche (2005), Ann. Probab. 33(3), 1643-1697. BBP transition. DOI 10.1214/009117904000000923.
2. Benaych-Georges, Nadakuditi (2011), Adv. Math. 227(1), 494-521. Outlier locations and overlaps. DOI 10.1016/j.aim.2011.02.011 (DOI re-check scheduled WP-A3).
3. Johnstone (2001), Ann. Statist. 29(2), 295-327. Tracy-Widom null for largest eigenvalue. DOI 10.1214/aos/1013203451.
4. Onatski (2010), Rev. Econ. Stat. 92(4), 818-835. Eigenvalue-ratio factor testing. DOI 10.1162/REST_a_00037 (re-check scheduled).
5. Athey, Bayati, Doudchenko, Imbens, Khosravi (2021), JASA 116(536). arXiv:1710.10251, DOI 10.1080/01621459.2021.1891924.
