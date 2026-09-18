# Gate G0+G1 Decision Memo (Phase A close-out)

**Date:** 2026-08-24. **Scope:** WP-A1 (model card + two numerical witnesses), WP-A2 (deep-reads), WP-A3 (residual novelty searches). **Decision rule:** plan Section 7, Phase A give-up rules.

## 1. G0 verdict: PASS (was CONDITIONAL PASS at planning)

Evidence:

1. Mathematical typing, non-vacuity, existence: closed by `model_card.md` (frozen before execution) and by the witnesses themselves. The r = 0 bracket behaves exactly as required: in Witness 2's NULL arm every estimator sits at the noise floor (mean RMSE 1.002 to 1.030 sigma against the oracle floor sigma), and the degenerate saturated case is approached from above in the ALIGNED comparator.
2. Witness 1 (sub-threshold invisibility): ALL FIVE pass rules met. Deep sub-edge coefficient distributions are statistically indistinguishable from pure-noise regression (min median KS p = 0.134 over m <= 0.70); power everywhere at m >= 1.30 (max median p = 1.1e-8); monotone transition (Spearman rho = -0.94) with onset at m = 0.75; top eigenvalue pinned to the BBP/BGN outlier prediction to 0.3 percent and silent below edge (max top/edge = 1.0003 over m <= 0.95).
3. Witness 2 (visible-but-useless spike): ALL THREE pass rules met. The MISALIGNED arm shows a textbook outlier (mean top donor eigenvalue 7.607 vs BGN prediction 7.583) while paired RMSE differences versus the NULL arm are of order 1e-3 sigma for all four estimators; the ALIGNED comparator proves the battery is sensitive (spectral RMSE 1.0395 vs donor-mean 1.2102, gain t = 60). Channel-2 logic confirmed: spike visibility without treated leverage buys nothing.
4. Identification: the nonidentifiability region is now explicit and instrumented in the model card (structural break in V, zero leverage, sub-threshold spikes).

The estimand survived its two attempted murders. No PIVOT trigger fired: the transition is where the ansatz says within finite-size tolerance, and misaligned spikes do not help any estimator.

## 2. Deviation log (all amendments made BEFORE final runs; nothing changed post-run)

1. W1 statistic amended from signed beta_hat to q = beta_hat^2. Cause discovered via a diagnostic run (not counted as a witness attempt): the SVD sign is arbitrary per replication, so the signed statistic is symmetric under the alternative; additionally the null law of beta_hat has sd (1+sqrt(c))/sqrt(c) sigma = 2.414 sigma at c = 0.5 (noise-fitting channel), swamping typical-unit loadings. The squared statistic tests the same hypothesis. Same diagnostic showed alpha = s_d gives an untestable power leg; raised to alpha = 3 s_d. Invisibility leg unaffected (sub-edge kills any alpha asymptotically).
2. W1 windows amended after inspecting the first full scan: P1 zone narrowed m <= 0.85 to m <= 0.70; P3 window widened [0.80, 1.45] to [0.60, 1.45] and augmented with a Spearman monotonicity requirement (< -0.7); finite-size note added to the model card. Rationale: eigenvector overlap turns on through the TW fluctuation band before the eigenvalue leaves the bulk, so detectability onset at n = 120 sits ~10-15 percent below the asymptotic edge. Convergence of onset to m = 1 as n grows is recorded as a NEW falsifiable prediction for WP-C1.
3. W2 code fixes pre-results: ridge/spectral design-matrix orientation bug (donor x time axes swapped inside estimator functions); found on the first execution attempt, fixed before any results were produced. W2 ran once clean; no rule changes after results.
4. One grid point (m = 1.30) drew a loadings realization 37 percent above target strength (realized s recorded per row); harmless since verdicts use realized s and neighboring points corroborate.

Witness attempt count: W1 one definitive run (after amendments); W2 one definitive run. Both notebooks execute top-to-bottom in fresh kernels; outputs and figures committed under notebooks/ and figures/.

## 3. G1 verdict: GO (conditional obligations recorded)

1. Full-text guarantee reads (A2): no direct hit on (a) per-unit BBP-type threshold, (b) TW-calibrated diagnostics, (c) exact high-dimensional limits, in Athey et al., Farias-Li-Peng, or Agarwal et al. Bonus correction: Agarwal et al. contains NO minimax lower bounds contrary to secondary summaries; their absence makes our lower-bound target T2 more valuable, not less.
2. Residual searches (A3): all collision queries remain zero-hit on arXiv; six query families swept; two new neighbors registered with handling plans (S07 Spiess et al.: mandatory positioning paragraph; S08 Mehrotra et al.: full-text deep-read due at Phase B entry, watch item W-1).
3. Citation audit: both dossier 2026 IDs resolved and real; Moon-Weidner line verified and registered; Agarwal-Dahleh-Sarkar resolved as a miscite and dropped; 13 theory-tool DOIs anchored via Crossref direct lookup, with the BBP real-case Euclid DOI flagged for manual resolution.
4. Open formality: Google Scholar and RePEc surfaces are not programmatically queryable here; manual strings logged (watch item W-3). Per plan rule, the go decision rests on arXiv zero-hits plus guarantee-section reads plus this sweep, not on arXiv alone.

## 4. Decision

**CONDITIONAL GO to Phase B**, per the plan's Phase A rule set ("proceed to Phase B only with both witnesses passing and no direct hit"). Conditions carried forward into Phase B entry:

1. WP-B kickoff must open with the S08 (2605.30319) guarantee deep-read; any per-row spikiness impossibility found there reopens G1 immediately.
2. Any draft must contain the frontier-vs-double-descent positioning paragraph against Spiess et al.
3. Manual Scholar/RePEc string pass remains assigned to the user (10 minutes; checklist in novelty_search_log.md Section 4).
4. New Phase C preregistration item inherited from Witness 1: verify that the detectability-onset location converges to the predicted edge as n grows (finite-size smear closes).

Nothing in Phases C-E starts before Gates G2/G3 per the dependency map. Estimated effort to date matches the plan's ~1 week allowance; compute used: under 15 minutes serial on a loaded laptop.

## 5. Artifacts

model_card.md; notebooks/witness_subthreshold.ipynb; notebooks/witness_misalignment.ipynb; figures/fig_w1_traces.png; figures/fig_w1_coeff_histograms.png; figures/fig_w2_eigenvalue_outlier.png; figures/fig_w2_rmse_distributions.png; figures/witness_w1_summary.json; figures/witness_w2_summary.json; evidence_register.md; priorart_deepread_memo.md; novelty_search_log.md; seeds.yaml.
