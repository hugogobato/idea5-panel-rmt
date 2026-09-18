# Prior-Art Deep-Read Memo (WP-A2)

**Date:** 2026-08-24. **Inputs:** full-text reads of arXiv 1710.10251 (v5), 2106.02780 (v2), 2109.15154 (v1); structured extraction with theorem-level anchors (see evidence_register.md S01-S03).

## The three kill questions

| Question | Athey et al. (S01) | Farias-Li-Peng (S02) | Agarwal et al. (S03) |
|---|---|---|---|
| (a) Per-unit recoverability threshold of BBP/spiked type | No. Signal enters via L_max and p_c in an order-rate Frobenius bound (Thm 2); no computable spikiness boundary; no lower bound at all | No. Their impossibility (Prop. 1) is tangent-space geometric identifiability, never tied to an eigenvalue spectrum; their minimax bound (Prop. 2) is for ATT, not per-unit trajectories | No. Entry-wise max-norm geometry; Assumption A6 is a balance regularity condition for the upper bound; zero lower-bound content in v1 |
| (b) Tracy-Widom / eigenvalue-calibrated diagnostics | No occurrences | No occurrences | No; Gavish-Donoho cited as related work only |
| (c) Exact high-dimensional limits (MP/TW/deterministic equivalents) | None; Matrix Bernstein + Massart concentration | None; operator-norm sub-Gaussian concentration only | None; non-asymptotic O_p conditioning on instance |

**Verdict: PASS. No direct hit on any of (a), (b), or (c) in the three closest papers' guarantee sections.**

## What each paper actually owns

1. Athey et al. own the unification of panel estimators as matrix-completion objectives (their Thm 1) and the staggered-adoption sampling-density story (p_c). Any paper on panel counterfactuals must position against them; none of their results constrain spectral detectability questions.
2. Farias-Li-Peng own rate-optimal ATT under general (even adaptive) intervention patterns. Notably useful to us: their tangent-space impossibility is a different (non-spectral) axis, so our frontier claim does not collide; and their minimax template is a candidate scaffold for our T2 lower-bound target.
3. Agarwal et al. own MNAR entry-wise imputation with arbitrary missingness and a CLT for the estimator under anchor-biclique structure. Their well-balanced-spectrum assumption is precisely the regime our spiked analysis replaces with explicit spike/edge calculus.

## Two additional neighbors found during A3 (abstract level, E2)

1. Spiess-Imbens-Venugopal (2305.00700): double-descent view of high-dimensional SC; "more donors help even past perfect pre-fit." Orthogonal axis to ours: they characterize weight-estimation variance in over-parameterized regressions, we characterize detectability of the treated unit's factor component against the MP sea. Mandatory positioning paragraph; watch item W-2.
2. Mehrotra-Tran-Vu-Zampetakis (2605.30319): improved PER-ROW matrix-completion guarantees for treatment effects. Abstract shows upper-bound improvements, not thresholds; scheduled for full-text guarantee deep-read at Phase B entry (watch item W-1). If it hides a per-row spikiness impossibility, G1 reopens.

## Implication for claim ledger

C1's remaining delta survives contact with the closest literature: nobody supplies a per-unit recoverability boundary as a function of spike strengths, aspect ratio, and treated-unit leverage, nor TW-calibrated rank/pre-trends diagnostics in the panel-causal setting. C3's incumbents are now precisely identified: CV-rank MC (no type-I control anywhere in S01-S03) and classical pre-trend t-tests plus the Roth-line pretest critiques (S11). The empirical lens C5 gains a new named incumbent analysis (S09 square-root lasso SC).

No change to the plan's architecture is required. Phase B may proceed once the G0+G1 decision memo records the two open watch items.
