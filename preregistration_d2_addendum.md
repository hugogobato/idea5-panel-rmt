# Preregistration D2 Addendum — distance-to-frontier analysis (WP-D2)

**Date frozen:** 2026-08-26, BEFORE any WP-D2 decisive computation on either application panel. **Parent documents:** `preregistration.md` (2026-08-24), `preregistration_c5_addendum.md` (2026-08-25), `gate_g3_memo.md` Section 7 (restricted claim set + GO to Phase D), `preprocessing_frozen.md` (WP-D1 freeze; its Section 8 leakage rules bind here). **Authority:** plan Section 7 Phase D (WP-D2 actions, pass/fail, give-up rules). Every threshold below is inherited from already-frozen Phase C instruments or fixed from theory alone; none is tuned toward any real-panel number, and no real-panel spectral quantity had been computed at freeze time. Failures are reported as failures.

---

## 1. Restricted claim under test

Per gate_g3_memo Section 7: C1 ships as "measured distance-to-frontier explains pre-treatment fit quality" on canonical SCM panels; C3 ships as the iid-calibrated diagnostic suite with documented scopes. The G3 condition discharged here: channel-2 (misalignment) flagging MUST integrate `alignment_energy` before any application claim (Section 4).

## 2. Analysis space and leakage

Primary space: **unit-centered pre-period levels**. For each unit row i of the donor matrix, subtract that row's own pre-period mean (leakage-safe; uses no post data and no other unit's data). Rationale: the spiked-model calibration (MP bulk, BBP inversion) presumes fluctuations around a stable mean; raw macro levels are dominated by cross-unit level heterogeneity that is not part of the noise sea. Common trends are NOT removed: they are factors (spikes) by the model's own definition, and spanning them is exactly how donors carry counterfactual information.

Secondary space (descriptive decomposition only, never gated): **first differences** of the pre window (T0-1 columns, each row self-centered). It answers a different question ("high-frequency comovement") than the primary space ("level-path span"); disagreement between spaces is reported as an interpretation datum, not a sensitivity failure.

Leakage: every spectral quantity is computed from donor PRE data plus the treated PRE row only. Donor POST windows may enter only through the two instruments explicitly designed for them (`pre_trends_post_test`, descriptive control in Section 8); treated POST outcomes never enter any D2 statistic. All panels enter exclusively through `applications.loaders`.

## 3. Channel 1: spike spectrum and distance d

For a donor pre matrix Y_D (n_d x T0, unit-centered), c = n_d/T0:

1. evals := descending eigenvalues of (1/T0) Y_D Y_D' (`scree`).
2. sigma2_0 := median(evals) / mu_MP(c), where mu_MP(c) is the MEDIAN of the unit-mean Marchenko-Pastur law with parameter c (computed by numerical CDF inversion; deterministic).
3. k_0 := gate(evals, sigma2_0): 0 if evals[0] <= 1.05 * sigma2_0 * (1+sqrt(c))^2, else largest-gap rank among top 4.
4. One bulk-refinement step (frozen): if n_d - k_0 >= max(10, ceil(n_d/4)), sigma2_1 := mean(evals[k_0:]) else sigma2_1 := sigma2_0. Justification: the MP trace identity makes the non-spike MEAN exactly sigma^2, so this removes mild upward bias from sub-k_0 contamination without iteration.
5. Final rank k := gate(evals, sigma2_1) with the same 1.05 rule, k_max = 4.
6. Spikes: s_j := invert_bbp(evals[:k], c) (BBP/BGN law lambda = 1+s+c+c/s); m_j := s_j/sqrt(c).

**Distance definition (frozen):** d := max_{j<=k} m_j when k >= 1; d := 0 when k = 0 (convention: no supercritical signal strength exists above the frontier). The frontier sits at d = 1 by construction; d < 1 can therefore only arise through bootstrap uncertainty around a silent or marginal spectrum.

## 4. Channel 2: alignment_energy integration (G3 condition)

e_obs := alignment_energy(Y_D_centered, y1_centered, k) = sum_{j<=k} <y1, v_j>^2 / ||y1||^2, with y1 the treated pre row centered by its OWN pre mean, projected on the OBSERVED top-k sample right singular directions.

Null distribution (frozen): because e is self-normalized by ||y1||^2, the pure-noise null is scale-free: draw G = 2000 Gaussian rows g ~ N(0, I_T0), recompute e(g) against the SAME fixed basis. p_align := (#\{e(g) >= e_obs\} + 1)/(G+1). This is the calibrated Witness-2 test: it detects whether the treated row holds statistically demonstrable leverage on the retained spike subspace beyond what an exchangeable noise row achieves. Channel 2 is CONFIRMED iff p_align < 0.05 (one-sided; the same 5% nominal discipline as every shipped instrument). r_align := e_obs / median(null) reported alongside.

## 5. Classification (frozen)

| Label | Rule |
|---|---|
| FRAGILE-INVISIBLE | k = 0 (silence gate fires; nothing supercritical exists) |
| RECOVERABLE | k >= 1 AND p_align < 0.05 AND donor-bootstrap q05(d) >= 1 |
| INCONCLUSIVE-SUBEDGE | k >= 1 AND p_align < 0.05 AND bootstrap CI for d straddles 1 |
| FRAGILE-SUBEDGE | k >= 1 AND p_align < 0.05 AND bootstrap q95(d) < 1 |
| FRAGILE-MISALIGNED | k >= 1 AND d point estimate >= 1 AND p_align >= 0.05 (visible spikes, no demonstrated treated leverage: channel-2 dead) |

## 6. Uncertainty (bootstrap, preregistered layers)

B1 (primary for classification): donor resampling, B = 500 draws of n_d rows with replacement from the unit-centered donor matrix; full pipeline of Sections 3-4 rerun per draw (k, d, p_align recomputed); percentile 95% interval for d. B2 (robustness column): circular moving-block resampling of the TIME index, block length 4, B = 200; same pipeline. Both layers are computed for both panels and every placebo.

## 7. Placebo battery and separation criteria

In-space placebos: for each donor i, placebo panel := treated unit i, donor pool := all units except {i, true treated} (true treated dropped, standard contamination guard). Placebos inherit the frozen D1 specifications verbatim: classic ADH predictor spec, windows, and V-search budget; the V multistart seed 60001 is IDENTICAL for every placebo fit (search structure held constant across placebos). Recorded per placebo: pre-fit RMSE in native units (smoking: packs; germany: % of treated-unit mean pre level), d_i, p_align_i, class_i.

Separation criteria (all frozen, one-sided):
- P1 (fit-ordering): Spearman rho(d_i, rmse_i) across placebos <= -0.30 (smoking, n=38) / <= -0.42 (germany, n=16; the ~5% one-sided critical value at that n), permutation p <= 0.05 (999 reps).
- P2 (treated separation): d_treated >= Q90({d_i}).
- P3 (date alarms): z_shift (shipped C5c instrument, verbatim defaults) run on donors-pre with pseudo-cutoffs smoking {1978, 1980, 1982}, germany {1972, 1976, 1980, 1984}; alarm rate (p < 0.05) must be <= 0.20.

## 8. Incumbent comparison and controls (descriptive, not gated)

Published-spec incumbents from WP-D1 (classic SCM weights/effects, synthdid SC/SDID/DID) plus mc_nn_cv (CV-rank nuclear-norm MC) are tabulated against the classification: their realized pre-fit versus what the frontier certifies as recoverable. Descriptive donor-post factor-law control: `pre_trends_post_test` on each panel (donor-pre basis vs REAL donor-post window, simulated iid finite-n null), reported with its scope caveat (Gaussian-noise calibration).

## 9. Sensitivity (gated stability)

Within the primary space, classification stability required in >= 3 of 4 combos: selector {primary gate, gate_lrv} x window {full pre, trimmed (drop 2 leading + 2 trailing pre-years)}. The ungated gap-ratio selector and CV-rank comparator are reported alongside. Difference-space results are excluded from this criterion (Section 2).

## 10. Predeclared novel-finding arms and verdict rules

NF arms (any one suffices for a novel finding, evaluated per panel):
- N1: the treated unit classifies into a FRAGILE label with the frontier explaining observed pre-fit weakness that informal inspection could not quantify.
- N2: a poor-looking pre-fit certified RECOVERABLE (provably frontier-safe despite appearance).
- N3: certification asymmetry: placebo units achieving pre-fit comparable to the treated unit classify differently (treated RECOVERABLE while comparable-fit placebos FRAGILE, or vice versa), demonstrating the frontier measures something equal-quality pre-fit cannot separate. Comparable := placebo pre-fit RMSE within +/-20% of the treated unit's pre-fit RMSE (frozen before any run).

WP-D2 verdict rules (frozen):
- Identification control: neither canonical treated unit may classify FRAGILE-INVISIBLE (both effects are real and precisely estimated; flagging them as unrecoverable artifacts would be positive-control failure => investigate before any reframing, KILL candidate per plan Phase D give-up rule 1).
- Applied-value PASS requires: NF on >= 1 panel AND P2 PASS on that panel AND P3 PASS on both panels AND Section 9 stability on the NF panel. If additionally P1 PASSes on the NF panel, verdict FULL PASS; otherwise PASS-WITH-LIMITATION (P1 failure documented).
- If all panels classify comfortably RECOVERABLE and no NF arm fires: applied value = certification-only; INCREMENTAL-ONLY pressure per plan give-up rule 4, honestly recorded.
- Placebo-separation failure everywhere (no P2/P3 separation on any panel): plan give-up rule 2 territory.

## 11. Seeds, environment, compute

Seeds registered in `seeds.yaml` (phase_d.wp_d2): alignment_null 60101 (shared stream offset per panel/placebo: base + index*10^6), donor_bootstrap 60102 (+index*10^6), time_bootstrap 60103 (+index*10^6), spearman_permutation 60104, z_shift 8880001 (C5 convention), post-test null pools 7770001 (C5 convention), placebo V-search 60001 for every fit. Environment unchanged from preprocessing_frozen.md Section 9. Compute: local only, minutes (no Colab trigger: everything << 2 h, << 4 GB).

## 12. Deviation log

D0 (operationalization pinned PRE-run, 2026-08-26, before any decisive computation): the phrase "poor-looking pre-fit" in arm N2 is operationalized as: the treated unit's classic-SCM pre-fit RMSE exceeds the 75th percentile of its own panel's placebo pre-fit RMSE distribution (worse-than-typical fit by the application's internal standard). No other amendment.

D1 (numerical safeguard, logged 2026-08-26 during the FIRST decisive execution, which crashed before writing any artifact; no decisive output existed yet): under time-block bootstrap draws with heavy column duplication the unit-space scatter is rank-deficient and `eigvalsh` returns ~-1e-17 dust on the zero half, which can push the median-based sigma^2 estimate numerically negative (math-domain crash). Safeguards: eigenvalues clamped at 0 before estimation, and sigma^2_0/sigma^2_1 floored at 1e-8 x lambda_1. These affect no non-degenerate configuration (dust clamping changes nothing above machine precision; the floor binds only in degenerate draws). First run aborted; the analysis was restarted from scratch afterward.
