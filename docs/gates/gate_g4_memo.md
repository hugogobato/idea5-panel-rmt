# Gate G4 Memo (Phase D close-out)

**Date:** 2026-08-26. **Preregistration:** `preregistration_d2_addendum.md`, frozen BEFORE any decisive computation (deviations D0-D1 logged there; D1 was a mid-run numerical safeguard applied after an abort that wrote no artifacts). **Execution:** single local run of `scripts/wp_d2_analysis.py` (~10 min): treated-unit pipelines, 38 + 16 in-space placebos (frozen ADH specifications, V seed 60001 throughout), bootstrap layers B1 = 500 donor / B2 = 200 block(4) time resamples, alignment nulls G = 2000, shipped z_shift/gate_lrv/post-test instruments verbatim. Artifacts: `results_d2/{distance_*,placebos_*,summary}.json|csv`, `figures/fig_d2_distance_to_frontier.png`; tests 40 green (32 prior + 8 new WP-D2 regressions).

## 1. Preregistered criteria vs outcomes

| # | Frozen criterion | Outcome | Verdict |
|---|---|---|---|
| 1 | Identification control: neither canonical treated unit FRAGILE-INVISIBLE | Smoking: RECOVERABLE (k=2, d=2856, donor-CI [1865, 4490], p_align 0.0005 = test floor, r_align 12.3). Germany: RECOVERABLE (k=1, d ~= 4.8e8, p floor, r_align 63.3) | PASS |
| 2 | NF arm N1 (fragile treated) | Neither treated unit fragile | FAIL |
| 3 | NF arm N2 (poor-looking yet certified safe) | Predicate never engaged: treated pre-fits are good (smoking RMSE 1.75 packs vs placebo Q75 4.00; Germany 1.48% of level, BETTER than all 16 placebos) | FAIL |
| 4 | NF arm N3 (certification asymmetry) | 53/54 placebos also RECOVERABLE (sole exception Connecticut, FRAGILE-MISALIGNED); comparable-fit placebo sets uniformly RECOVERABLE | FAIL |
| 5 | P1 fit-ordering rho(d_i, rmse_i) <= -0.30/-0.42 | Smoking: rho = -0.53, perm-p <= 0.001 (n=38). Germany: rho = +0.26, perm-p 0.82 | PASS / FAIL |
| 6 | P2 treated separation: d_treated >= Q90(placebo d) | Smoking: 2856 vs Q90 2893 (42nd percentile). Germany: 4.755e8 vs Q90 4.756e8 (indistinguishable) | FAIL / FAIL |
| 7 | P3 z_shift date-alarm rate <= 20% | Smoking 0/3. Germany 3/4 (p = 0.045 at pseudo-cutoffs 1976/1980/1984) | PASS / FAIL |
| 8 | Sensitivity stability >= 3/4 combos | Both panels 2/4: labels stable across trimmed windows under the primary gate, but gate_lrv silences (k=0) both windows on both panels | FAIL / FAIL |
| 9 | Composite applied-value PASS (NF AND P2 AND P3 AND stability on >= 1 panel) | Not met anywhere | FAIL |

## 2. What the lens did and did not see (substance, not just gates)

F1 (saturation): in unit-centered LEVELS both panels sit astronomically above the frontier (lambda_1/edge ~= 51 smoking, ~= 2100 Germany; m-hat up to 2856 and ~10^8). The BBP inversion is invoked orders of magnitude outside its calibrated regime, so ABSOLUTE d carries no unit-discriminating information and the classification saturates at RECOVERABLE essentially everywhere. This is itself the central empirical answer: shared stochastic trends make canonical SCM panels deep-supercritical BY CONSTRUCTION; the recoverability frontier never binds on them, which is precisely why SCM can work there at all.

F2 (relative lens validated once): within the smoking panel, frontier position orders pre-fit quality across 38 placebos exactly as theory predicts (rho = -0.53, permutation p <= 0.001). The German panel fails this (rho > 0): with one dominant world-growth spike and n = 16, cross-country fit variation is not spectrum-ranked.

F3 (channel-2 integration, discharging the G3 condition): the alignment test hits its evidence floor on both treated units and surfaces exactly one real Witness-2 phenotype among 54 placebo units: Connecticut shows a supercritical k=2 spectrum with NO demonstrable treated leverage (p = 0.105) and a worst-quintile pre-fit. Visible-spike-without-leverage exists in real data and the integrated flag catches it.

F4 (transport risk made concrete): the donor-post factor-law control rejects post-treatment donor-law stability on BOTH panels at the simulated-null floor (p ~= 0.003). Pre-period certification does not certify post-window transport; Assumption A4 is visibly strained on the canonical applications. This feeds T3 directly.

F5 (robust-gate scope boundary confirmed out of the box): gate_lrv reads trend-dominated rows as near-unit-root noise and silences both panels (the documented serial-correlation limitation manifesting exactly as diagnosed at G3). No repair attempted; scope stands.

F6 (difference-space decomposition): differencing does NOT rescue variation in d (smoking k=1, d ~= 360, aligned; Germany k=1, d ~= 1.8e6, aligned): even growth-rate comovement is massively supercritical. Fragility in SCM practice therefore cannot be attributed to signal scarcity on these panels; where it occurs it runs through break and misalignment channels (F3, F4), not through sub-frontier spectra.

Incumbent descriptives: mc_nn_cv reproduces the smoking effect (-18.5 vs published -19.6) but fails badly on Germany (+2293, loses the level path); gated spectral SC lands at -29.6 (smoking) and -1580 (Germany ~= ADH -1588) - descriptive only, no gates attached.

## 3. Strongest baseline's best case and the proposed lens's failure regions

The strongest incumbent analysis (published SCM/SDID specs plus CV-rank MC plus eyeballed pre-fit) is unbeaten here: every preregistered claim the frontier lens could exclusively own on these two panels failed (N1-N3, P2, German P1/P3). Failure regions of the lens: any panel whose scatter is dominated by shared stochastic trends (absolute d uncalibrated); small panels with one dominant common factor (fit variation not spectrum-ranked); instruments whose nulls assume stationarity (gate_lrv, z_shift under trending regimes - the German date-alarm failures are exactly this).

## 4. Plan give-up-rule mapping

Rule 1 (KILL: diagnostic flags true effects as artifacts): NOT FIRED - both positive controls classified correctly (RECOVERABLE + maximally aligned).
Rule 2 (KILL: distances uniformly uninformative AND placebos inseparable): NOT FIRED - the first conjunct fails (F2 smoking ordering; F3 misalignment catch; F4 transport detection are genuine signal).
Rule 3 (PIVOT: lens works only for a different class than advertised): SUBSTANTIALLY INDICATED - the working residue is relative/within-panel ordering and channel-2 flagging, not absolute distance classification.
Rule 4 (INCREMENTAL-ONLY: purely descriptive additions): APPLIES TO THE C5 CLAIM - the applied yield is certification and scope-definition, not a changed reading of either case study.

## 5. Decision

**WP-D2 formal outcome: the preregistered novel-finding arms FAILED on both panels; the applied-value PASS was not met. Per plan Section 7 Phase D, C5 is demoted: no claim that the frontier lens changes the interpretation of the California smoking or German reunification applications survives this gate.**

**Gate G4 disposition: PROCEED TO PHASE E UNDER RESTRICTED CLAIMS, no reruns required.** Specifically:

1. C1/C3 stand on Phase C evidence; nothing in WP-D2 contradicts them, and the identification controls PASSED. The paper's spine remains the threshold/calibration layer (T1-T3).
2. C5 ships as a certification case study ONLY: "on two canonical panels the frontier certifies recoverability (both treated units deep-supercritical and maximally aligned), quantifies why it never binds under shared stochastic trends (F1/F6), demonstrates the channel-2 misalignment phenotype in the wild (Connecticut, F3), and exposes post-window donor-law drift as the operative transport risk (F4)." The distance-to-frontier reporting standard is withdrawn as a headline contribution; the smoking cross-sectional ordering (F2) is reported as supporting evidence, not a claim.
3. T3 gains two mandatory scoped limitations from D2: donor-law drift detectability (F4) and the trend-saturation regime boundary of all shipped calibrations (F5/F1).
4. User decision point preserved by the plan: if the bar excludes incremental applied contributions entirely, the C5 thread is cut at no loss to C1/C3; terminating the applied thread does NOT trigger project termination.

Every criterion above traces to `results_d2/`; figures regenerate via `scripts/wp_d2_analysis.py` (deterministic seeds registered in seeds.yaml, phase_d.wp_d2).
