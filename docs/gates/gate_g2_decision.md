# Gate G2 Decision Memo (Phase B close-out)

**Date:** 2026-08-24. **Scope:** Phase B entry conditions (S08 deep-read, Spiess positioning), WP-B1 (frontier ansatz + symbolic checks), WP-B2 (reference implementation library + tests), WP-B3 (pilot cost model). **Decision rule:** plan Section 7, Phase B give-up rules.

## 1. Entry conditions carried from G0+G1: both discharged

1. S08 (Mehrotra-Tran-Vu-Zampetakis, arXiv 2605.30319v1) full-text guarantee read completed at kickoff (E3, anchored in `evidence_register.md`). No impossibility or lower-bound statement implying a per-row spikiness threshold; their necessity remarks concern incoherence/SNR as regularity for THEIR upper bound. Answers: (a) NO, (b) NO, (c) NO. G1 not reopened; W-1 CLOSED. Bonus tool registered: their Theorem B.2 row-wise truncated-SVD perturbation bound, candidate input for T1/T2.
2. Spiess-Imbens-Venugopal full text read (E3). Purely mechanical model-averaging theory; zero spectral content. Canonical positioning paragraph frozen in the register (W-2 CLOSED): frontier varies signal strength/alignment at fixed estimator class; descent curves vary complexity at fixed signal; orthogonal axes with a joint empirical test (descent flattening at the sigma floor).

## 2. WP-B1: formula standing

`frontier_ansatz.md` states the conjectured DE risk for hard-threshold spectral SC:

```
F = 1
  + sum_{j in K} (alpha_j^2/sigma^2) [ (1 - zeta_j tau_j)^2 + zeta_j/lambda_j ]
  + sum_{j in K} (s_j zeta_j + 1) / (T0 lambda_j)
  + sum_{j not in K} (alpha_j^2/sigma^2)
```

with `lambda` (BBP/BGN outlier location), `zeta` (BGN overlap), `tau = sqrt(s/lambda)`; every ingredient mapped to a verified source (`theory_targets.md` stub created; T4 marked cut-by-default). Special cases verified symbolically (`code/check_frontier_ansatz.py`, ALL PASS): r=0 gives exactly the noise floor; s -> inf collapses to the floor; alpha=0 kills all spike terms exactly; the frontier is continuous at m=1 with a KINK (included-spike value connects to truncation value) and total risk is strictly decreasing on the supercritical side with maximum exactly at the edge. Calibration point: at Witness 2's ALIGNED cell the formula predicts RMSE 1.0369 sigma vs measured 1.0395 (delta 0.25 percent); recorded as one point, not evidence.

Deviation log for WP-B1: two intermediate monotonicity statements were corrected before freeze (a harness scaling slip dropped the 1/T0 factor on channel 3 in test code, not in the formula); the final document claims only what the checks verify. No rule changes after results; no revision cycle consumed.

## 3. WP-B2: green tests, honest battery

Package `code/scm_frontier/` (dgps, estimators, diagnostics) implementing the DGP ladder knobs and all seven estimators plus diagnostics (scree, gap-ratio rank selector with TW-style silence gate, TW statistic, BBP inversion, alignment energy). `tests/test_estimators.py`: 10 tests, ALL PASS (final run 390 s under load; leakage enforced structurally by API plus behavioral determinism tests). Scientific rules verified: oracle floor exact; r=0 at floor; infinite-spike near-oracle for tuned methods (donor-mean excluded by design); spectral beats donor-mean paired t >> 3; trivial baseline never usually beats tuned methods; simplex solver stable (feasibility-based convergence criterion after diagnosing SLSQP status-8 stalls-at-optimum).

Smoke run: 20 reps x 1 cell x all methods in 35.5 s (PASS vs 60 s target; machine under load ~41).

Deviation log for WP-B2 (all fixes pre-results, i.e., before any Phase C run):
1. ridge_sc broadcasting bug (donor-feature axis), fixed on first execution.
2. SDID rewritten twice: (i) missing analytic Jacobian made SLSQP numerically differentiate ~160 dims and stall far from optimum; (ii) the weighted-least-squares design-scaling formulation silently rescaled the ridge penalty relative to the data term by ~T0; replaced with direct row-weighted residuals. Final reduced port behaves correctly in null, supercritical, and infinite-spike cells (RMSE 1.02 / 1.26 / 1.02 respectively at the probe cells) and is documented as a reduced port without inference.
3. mc_nn_cv warm-started down the lambda grid; tolerance loosened to 1e-4 (CV selection does not need tighter).
4. sdid default sweeps set to 1 after verifying the second alternating pass changes nothing while doubling cost.
5. Test bounds recalibrated to PAIRED method-minus-oracle differences after discovering that absolute RMSE levels near the floor are dominated by the realized post-noise draw of each seed batch (oracle itself averaged 0.9406 sigma on seeds 10000-10009 versus 0.9953 over 200 seeds).

No PIVOT trigger: the simplex solver is stable in its home turf (>99% feasibility-converged), so the ridge-only fallback scope was not needed.

## 4. WP-B3: sane pilot, classification recorded

`pilot_cost_report.md` (raw JSON in `figures/pilot_costs.json`). Production decisive cell (n=T0=250, T_post=125): spectral 28 ms/rep, ridge 216 ms, scm 4.35 s, mc 4.79 s, sdid 11.0 s; peak RSS 93 MiB. Scaling exponents vs n*T0: spectral 1.31, mc 1.47, sdid 1.67, scm 1.79. Classification under the plan's rule (LOCAL iff < 2 h @ <= 8 workers AND < 4 GiB): C2(ii)/C2(iv)/C2(iii) LOCAL; C1 grid, C1 slices, and C2(i) COLAB as measured on a heavily loaded machine (conservative). Per-rep stability within the 30% model rule for all expensive methods. Inputs for preregistration (cell size, method subset for the full sweep, rep counts) are listed neutrally in the report; no decision taken in Phase B. Seeds registered in `seeds.yaml`.

## 5. Decision

**GO to Phase C**, per the plan's Phase B rule set ("GO: green tests, sane pilot, formula standing"). All three legs satisfied:
1. Formula standing: ansatz internally consistent, special cases symbolic-PASS, kink signature verified, one calibration point at 0.25 percent.
2. Green tests: 10/10 including the scientific ordering rules and the leakage guard.
3. Sane pilot: cost model accepted within its pass rule; LOCAL/COLAB table recorded.

Phase C obligations inherited (to be honored in `preregistration.md` BEFORE any decisive run):
1. Freeze the cell-size / method-subset / rep-count trade-off using Section 5 of the pilot report, then freeze `preregistration.md`.
2. Include the inherited Witness-1 falsifier: detectability-onset location converges to m=1 as n grows.
3. Parameterize treated leverage as theta_j = alpha_j^2/sigma^2 directly (ansatz scaling caveat): non-vanishing frontier requires non-vanishing treated share.
4. Carry the positioning paragraph (register S07) into any draft; W-2 text is frozen.
5. Re-measure the bootstrap-cost multiplier before C2(iii) instead of relying on the assumed 50x.

Nothing beyond Phase C starts before Gate G3 per the dependency map.

## 6. Artifacts added this phase

frontier_ansatz.md; theory_targets.md; code/scm_frontier/ (dgps.py, estimators.py, diagnostics.py, __init__.py); code/tests/ (conftest.py, test_estimators.py); code/check_frontier_ansatz.py; code/run_smoke.py; code/run_pilot.py; pilot_cost_report.md; figures/pilot_costs.json; seeds.yaml updated; evidence_register.md updated (S07 E3 + positioning paragraph, S08 E3, W-1/W-2 closed).
