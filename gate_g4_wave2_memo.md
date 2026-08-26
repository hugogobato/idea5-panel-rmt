# Gate G4 Wave-2 Memo (hypothesis-driven panel extension)

**Date:** 2026-08-26. **Authorization:** preregistration_d2b_addendum.md (frozen BEFORE any Wave-2 computation; declared selection bias: panels chosen because the frontier was EXPECTED to bind). **Execution:** `scripts/wp_d2b_wave2.py`, single-threaded BLAS under heavy external load; artifacts `results_d2/wave2_{basque,basque_placebos,degradation}.*`, `figures/fig_d2_wave2.png`. Basque data: R `Synth` package extract (CRAN mirror), sha256-pinned `data/raw/basque_synth_ag2003.csv`, provenance E2 (register D06).

## 1. Frozen gates vs outcomes

| # | Frozen item | Outcome | Verdict |
|---|---|---|---|
| 1 | W1 treated fragile (FRAGILE-MISALIGNED / SUBEDGE family) with CI support | RECOVERABLE: k=1, d=5.83, donor-CI [4.46, 7.53], p_align at floor (r_align high) | FAIL (hypothesis rejected) |
| 2 | W1 P3 date alarms <= 0.20 (pseudo-cutoffs 1969/1971) | 0/2 alarmed | PASS |
| 3 | W1 P2 separation d_treated >= Q90(placebo d) | PASS: 5.833 at the top of the placebo range (placebos 4.68-5.88) | PASS |
| 4 | W1 P1 fit-ordering rho <= -0.45 | rho = -0.24, perm-p n.s.; pre-fit variation across Spanish regions not spectrum-ranked | FAIL |
| 5 | W1 sensitivity stability >= 3/4 | gate_lrv again silences trend-dominated rows; 2/4 | FAIL |
| 6 | N3 comparable-fit asymmetry | sole comparable placebo (Castilla-La Mancha) also RECOVERABLE | FAIL |
| 7 | **W1_PASS composite** | | **FALSE** |
| 8 | W2 control: full panel loud AND some size >= 4 reaches >= 50% silent draws | Full loud (k >= 1 everywhere); silence only at n_d = 3 (10%); median d-hat declines 2856 -> ~700 but never near the edge | FAIL as designed transition; informative |

Deviation D-W3 (logged in the addendum): the first W1 run used the true-treated row as every placebo's regression target (leave-one-donor-out refits of Basque; detected via degenerate duplicate pre-fits). No conclusion was drawn from it; the corrected run restarts from scratch and is the only decisive W1 record.

## 2. Findings

F7 (moderate regime EXISTS on real macro data): Basque lands at d ~= 5.8 - three orders of magnitude below smoking/Germany saturation and INSIDE the regime where BBP inversion is roughly calibrated. Non-saturated distances are therefore attainable without simulation; Wave-1's degeneracy was a property of those two panels' shared-stochastic-trend dominance, not of the lens.

F8 (first natural P2 separation): Basque certifiably separates from Spanish placebo regions on frontier distance (P2 PASS), something neither canonical panel could deliver. The lens discriminates when the spectrum lives at a discriminating scale.

F9 (the fragility hypothesis for "bad pre-fit" applications rejected here): Basque's mediocre OUTCOME-ONLY pre-fit (RMSE ~= 1.46 thousand USD, ~26% of level - the published AG2003 fit uses predictor-based V machinery we deliberately did not replicate, deviation D-W1) is NOT explained by frontier proximity or misalignment: the panel is supercritical and maximally aligned. Where classic outcome-only SCM fits badly on an aligned supercritical panel, the failure is estimator-side (level/growth mismatch that predictor reweighting fixes), not information-side.

F10 (pool thinning does not approach the frontier): removing up to 34 of 38 smoking donors leaves median d-hat in the hundreds; silence appears only at n_d = 3. Frontier distance is governed by the factor structure (a shared trend survives any subset), not by donor count - consistent with F1 and directly relevant to how the next dataset hunt should screen candidates.

## 3. Disposition

W1_PASS = FALSE per frozen rule: no natural fragile SCM panel has been found in the candidate set {smoking, Germany, Basque}; C5 remains demoted exactly as gate_g4_memo.md Section 5 recorded. No KILL trigger fired anywhere; positive controls keep passing.

What Wave-2 adds to the program: (i) proof of attainability (F7) and separability (F8) on natural panels, which upgrades the lens from "saturated certification" to "functional instrument awaiting a binding application"; (ii) two concrete screening criteria for the dataset hunt now underway: prefer panels WITHOUT a dominant shared stochastic trend (else d saturates) and do NOT expect pool size alone to create fragility (F10); target small-n/high-noise/weak-factor settings with documented poor pre-fit.

Next step owned by the user: the external deep-research sweep over published SCM/SDID replications for a genuinely sub-frontier application.
