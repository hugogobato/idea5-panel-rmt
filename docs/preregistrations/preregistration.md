# Preregistration (Phase C) — frozen before any decisive run

**Project:** Idea 5, spectral recoverability frontiers for panel causal inference.
**Date frozen:** 2026-08-24 (before any Phase C decisive cell was executed).
**Supersedes:** plan Section 7 Phase C defaults where explicitly deviated below; every deviation is itemized in Section 9.
**Status:** FROZEN. No rule, grid point, metric, or pass rule in this document may change after the first decisive replication is generated. Failures are reported as failures.

---

## 1. Primary metric

Per replication: `rmse_rep = sqrt(mean_{t in post}(yhat_t - y*_t)^2) / sigma`, where `y*_t = L_1t + E_1t` is the realized untreated trajectory of the treated unit (model_card.md Section 4) and `sigma = 1` by construction. Cell statistic: mean over replications; median and IQR reported alongside; per-method paired differences against oracle recorded for all comparisons.

## 2. Secondary metrics (frozen definitions)

1. ATT bias: `mean_t(yhat_t) - mean_t(y*_t)` per replication; reported as mean and IQR over reps.
2. Rank-selection accuracy: gated selector `k` (silence gate + largest gap ratio, k_max = 4) vs true `r`; ungated ablation (`k` by gap ratio only); CV-rank comparator (Section 5.3).
3. Empirical size/power of the spectral pre-trends test at nominal 5%: rejection rate of (a) the simulated-TW-null statistic `Z_tw` and (b) the circular block-bootstrap statistic `Z_boot`; classical comparator: OLS t-test of `y1_pre` on `(1, t)` trend coefficient.
4. Coverage of the 95% spectral prediction interval: per rep, fraction of post periods with `y*_t` inside `yhat_t +/- 1.96 * sd_hat * sqrt(1 + sh_t' (S'S)^{-1} sh_t)` (k = 0 case: `ybar1 +/- 1.96 * sd_hat * sqrt(1 + 1/T0)`), `sd_hat^2 = RSS_pre/(T0 - k)` from the treated-on-scores pre regression.
5. Runtime and peak RSS per method per rep (cost audit, feeds any re-shard decision).

## 3. Decisive grid (WP-C1), approved trade-off bundle

Geometry: **equal-entries cells**, `n * T0 ~= 25,600` with `c = n/T0` exact:

| c | n | T0 | T_post = floor(T0/2) | n*T0 |
|---|---|---|---|---|
| 0.25 | 80 | 320 | 160 | 25,600 |
| 0.50 | 113 | 226 | 113 | 25,538 |
| 1.00 | 160 | 160 | 80 | 25,600 |
| 2.00 | 226 | 113 | 56 | 25,538 |
| 4.00 | 320 | 80 | 40 | 25,600 |

Spike multiplier grid: `m = s/sqrt(c) in linspace(0.2, 3.0, 15) UNION {0.9, 1.1}` (17 points, sorted; m = 1 exactly on-grid, +/- 0.1 spacing at the edge). Spikes equal-strength across factors: `s_j = m * sqrt(c)` for all j.

Alignment arms (G2 obligation: leverage parameterized as `theta_j = alpha_j^2/sigma^2` directly):

| arm | alignment | theta |
|---|---|---|
| orthogonal | "none" | 0 |
| partial | "first" | 0.25 |
| full | "first" | 1.0 |

Factors: r in {1, 3} (r = 3 uses three equal spikes `s_j = m*sqrt(c)`).

Method set, full sweep: donor_mean, scm_simplex, ridge_sc (CV over logspace(-1,4,11), 4-fold), spectral_sc GATED (sigma, c supplied; silence gate at 1.05x MP edge), spectral_sc UNGATED (ablation, same SVD), mc_nn_cv (soft-impute, held-out 10% CV over logspace(0,3,11)). SDID restricted to the preregistered subgrid: c in {0.5, 2} x full arm x r = 1 x all 17 m points (34 cells).

Cell count: 5 c x 17 m x 3 arms x 2 r = 510 cells x 500 reps (+ SDID on 34 of them, + C2(i) battery folded into the same fleet, Section 4).

Seeds: replication i in a cell uses seed `10000 + i`, i = 0..499, identical across methods (all estimators run on the same generated panel within a rep). Seed reuse across cells is by design (cell config differs); derivation formulas live in each notebook header and `seeds.yaml`.

Equal tuning budget convention (frozen): adaptive methods are ridge_sc and mc_nn_cv, each with an 11-point penalty grid and 4-fold/held-out CV respectively as implemented in WP-B2; spectral selectors are rule-based (no tuning); SCM/SDID fixed algorithms. This is declared the equal-budget convention; no method gets extra tuning passes.

## 4. Honesty battery (WP-C2)

1. **C2(i) null battery** (folded into the C1 fleet): r = 0 panels (alignment "none") at the five grid geometries plus one production-size iid null (n=T0=250, T_post=125, full method set incl. SDID). 500 reps each. Pass: pairwise method-vs-oracle excess RMSE within MC bands (no incumbent beats oracle, none trails donor_mean by > 0.05 sigma systematically); gated rank selector returns k = 0 in >= 94% of reps (size <= 6% at nominal 5%); Z_tw and Z_boot reject <= 6% on null panels.
2. **C2(ii) baseline-favorable**: dense weak factors favoring nuclear-norm MC. Geometry (n,T0,T_post) = (160,160,80), alignment "all", theta_total = 1.0 spread equally; four cells: (r=8, m=0.6), (r=8, m=0.8), (r=16, m=0.6), (r=16, m=0.8); each s_j individually subcritical. 500 reps, full method set. Pass: whichever family wins is reported honestly; the battery FAILS the phase only if spectral_sc is dominated by donor_mean here (mechanism inverted, not merely beaten).
3. **C2(iii) calibration under weak dependence**: noise laws {gaussian control, AR(1) rho=0.3, AR(1) rho=0.7, heteroskedastic het_ratio=4} x {null, spiked m=1.6 theta=1 aligned}, geometry (160,160,80), 500 reps, diagnostics-only battery (gated/ungated k, Z_tw, Z_boot, classical t-test, plus rmse for spectral/ridge as robustness riders). The notebook FIRST re-measures the bootstrap cost multiplier (10 timed reps with and without B=200 bootstrap) and records it before any calibration output (G2 obligation 5). Pass: Z_boot size in [3%, 8%] under AR(1)/het nulls AND power >= 80% in spiked cells; Z_tw size drift documented; if Z_tw fails but Z_boot passes, the diagnostic ships bootstrap-native (recorded scope change, not failure).
4. **C2(iv) structural break** (A4 violation detectability): geometry (160,160,80), m=2.0, four cells: (delta=1, aligned theta=1), (delta=2, aligned theta=1), (delta=2, orthogonal theta=0), (delta=0 control, aligned theta=1). 500 reps, full method set + diagnostics. Pass: Z_boot power >= 80% in delta >= 1 break cells, control size <= 8%, and every estimator's RMSE degrades vs its delta=0 pair.

## 5. Diagnostic head-to-head (WP-C3)

1. Metric hierarchy (fixed): primary = gated-selector rank accuracy minus CV-rank comparator accuracy on C1 cells with true r in {1,3}, m >= 1.2; secondary = pre-trends size/power comparisons from C2(iii)/(iv).
2. CV-rank comparator (frozen): choose k in 0..4 minimizing 4-fold contiguous-block CV MSE of the k-PC regression of y1_pre on donor scores (same folds construction as ridge_cv_seed = 1234 stream); computed inside every C1 rep at negligible cost.
3. Ablation isolating the claimed mechanism: gated vs ungated selection (TW-style silence gate on/off), and Z_tw vs Z_boot (iid-parametric vs dependence-robust null).
4. Pass: gated beats CV-rank accuracy by >= 10 percentage points in at least the mid-aspect columns (c in {0.5, 1, 2}) on supercritical cells (m >= 1.2); Z_boot strictly out-sizes the classical t-test under rho >= 0.3. Fail: no regime where the diagnostic wins => INCREMENTAL-ONLY pressure on C3 per plan.

## 6. Practical-bite criterion and falsifiers (frozen)

1. **Kink criterion (primary falsifier of the ansatz).** For each c column, primary curve = mean rmse of GATED spectral_sc, full arm, r = 1, over m in [0.6, 1.6]. Kink estimate = m-grid point maximizing the discrete second difference `RMSE[m+1] - 2 RMSE[m] + RMSE[m-1]`. PASS iff `|m_kink - 1| <= 0.15` in >= 80% of the five c columns (i.e., >= 4 of 5). Sensitivity (non-decisive): two-segment piecewise-linear breakpoint fit reported alongside.
2. **DE overlay (descriptive, T1 support):** closed-form F (frontier_ansatz.md Section 3, evaluated at each geometry's T0, c, theta) plotted atop simulated gated curves for the full arm; systematic deviation beyond MC bands at multiple c triggers the ansatz-revision loop (max two documented revisions naming the corrected ingredient; third failure kills C1 per plan).
3. **Onset-convergence check (inherited Witness-1 falsifier, dedicated slice):** c = 1, sizes (n,T0) in {(81,81),(121,121),(161,161),(241,241),(361,361),(541,541)}, T_post = T0/2, m grid linspace(0.55, 1.65, 23) (Witness-1 grid for comparability), 300 reps/point, full arm theta = 1. Detectability onset per size = smallest m with median KS p-value < 0.01 using the W1 squared-coefficient statistic q = beta_hat^2 against a 600-draw pure-noise pool per size (seeds 51101 master / 51102 pools). PASS iff onset at the largest size lies in [0.90, 1.10]. Diagnostics-only compute (seconds per rep).
4. **Practical bite:** frontier has bite iff some incumbent {scm_simplex, ridge_sc, mc_nn_cv} attains mean rmse >= 2.0 sigma in at least one preregistered substantive region (defined as any (c column, arm, r) cell family with c >= 1) while the gated selector flags that region with rate >= 80%, and flags no more than 20% of clearly-recoverable cells (m >= 1.5, theta >= 0.25, aligned). Flag event := gated k = 0.
5. **WP-C5 decision rules:** unchanged from plan Section 7 (GO / PIVOT / INCREMENTAL-ONLY / KILL), applied to these frozen criteria.

## 7. Execution vehicle (user-approved deviation)

All Phase C packages run as self-contained Google Colab notebooks generated from this repo (user decision 2026-08-24; supersedes the pilot's LOCAL classifications). Notebook families: `nb_c1_shardNN_of40.ipynb` (main grid + folded null battery, cost-balanced via measured per-rep model), `nb_c2ii_baseline_favorable.ipynb`, `nb_c2iii_calibration.ipynb`, `nb_c2iv_break.ipynb`, `nb_c4_scaling.ipynb`, `nb_onset_slice.ipynb` (45 notebooks total). Every notebook: generates all data internally from registered seeds, appends checkpoint CSV rows every 25 reps, writes a JSON metadata header (grid coordinates, library versions, timestamps), ends with the mandatory download fallback, targets <= 4 h runtime (hard cap well under Colab's session limit). Aggregation happens locally only after `merge_shards.py` verifies completeness and checksums against `shard_manifest.yaml`; no figure is ever produced from a partial grid. WP-C4 scaling runs on Colab too (spectral-family timing up to n = T0 = 5000, memory-light; simplex-family costs extrapolated from WP-B3 exponents, stated as such).

## 8. Compute and determinism notes

Panels are bitwise-reproducible across machines (numpy PCG64 seeded draws); floating-point results may differ at LAPACK level between Colab BLAS and local BLAS, which is acceptable for Monte Carlo statements and recorded in shard metadata. Checkpoints make a crashed shard lose at most one 25-rep chunk.

## 9. Deviation log relative to plan draft (all made BEFORE first decisive run)

1. Equal-entries geometry replaces square n = T0 cells (cost; consequence: the 1/T0 variance channels vary across c columns, which the DE overlay tests rather than assumes away).
2. m-grid refined by {0.9, 1.1}; everything else identical to linspace(0.2, 3, 15).
3. Alignment axis replaced by explicit theta parameterization (G2 obligation 3): orthogonal/partial/full = theta 0 / 0.25 / 1.0, direction "first"; channel-2 logic already witness-tested and re-covered by orthogonal arms + C2(iv) orthogonal cell.
4. SDID moved to a preregistered subgrid (pilot Section 5 lever b); baselines SCM + CV-MC remain in every cell.
5. Reps kept at 500 (pilot Section 5 lever c NOT used); precision preserved for the +/- 15% kink criterion.
6. All packages routed to Colab notebooks (user decision), including pilot-LOCAL-classified batteries.
7. Onset-convergence obligation served by a dedicated cheap slice instead of the main grid (cleaner: joint (n, T0) growth at fixed c).
8. Pre-trends test operationalized as post-residual TW statistic with simulated finite-n null (`Z_tw`) plus circular block-bootstrap null (block length 10, B = 200, pseudo-cutoff within observed pre window) (`Z_boot`); classical comparator = linear-trend t-test. Definitions frozen above; implementation must match them exactly.
9. The manual Scholar/RePEc novelty pass (W-3) remains open on the user's side; per plan it does not gate Phase C execution and stays logged as a residual formality.
10. CORRECTION (pre-run, caught during generator construction 2026-08-24): the initially drafted geometry table listed (n, T0) pairs violating the stated aspect ratios; replaced with the exact-c table above before any replication existed. No other section affected.

Nothing beyond Phase C starts before Gate G3.
