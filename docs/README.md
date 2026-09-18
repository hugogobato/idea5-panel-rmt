# Documentation index

Research trail for **spectral-frontier**, internally "Idea 5" in the parent research program.
The code lives in `code/`, experiment outputs in `results_*/`, and figures in `figures/`.
Every document below is part of the provenance chain: the plan froze the gates, the
preregistrations froze the analyses, and the gate memos record what actually happened,
including negative results.

## Where to start

1. `plan/Idea5_Panel_RMT_Research_Plan.md`: the full phased program, claim ledger, simulation
   design, risk register, and give-up rules. Start here for the big picture.
2. `model/model_card.md`: the formal model, estimand, calibration conventions, and assumption
   ledger.
3. `model/frontier_ansatz.md`: the conjectured deterministic-equivalent risk formula, its
   symbolically verified special cases, and its stated limitations.
4. `gates/gate_g3_memo.md`: the decisive simulation gate, with the restricted claim set that
   survived.
5. `applied/applied_study_closure_memo.md`: the applied casebook scope and its negative results.
6. `../paper/powell_summary/powell_summary.pdf`: a four-page intuitive summary with the main
   figures.

## Folders

### `plan/`
- `Idea5_Panel_RMT_Research_Plan.md`: master plan (gates G0 to G6, work packages, compute policy,
  Colab sharding rules, artifact map).

### `model/`
- `model_card.md`: model, estimand, spiked-covariance conventions, assumption ledger A1 to A7,
  witness preregistrations.
- `frontier_ansatz.md`: the DE risk formula F, special-case algebra, and source map (tagged
  CONJECTURE; proof is Phase E target T1).
- `theory_targets.md`: T1 to T4 proof packages mapped to verified sources and numerical
  falsifiers.

### `preregistrations/`
- `preregistration.md`: frozen Phase C simulation design, metrics, pass rules, and seed policy.
- `preregistration_c5_addendum.md`: repair-and-confirm battery after the first kink estimator
  defect.
- `preregistration_d2_addendum.md`, `preregistration_d2b_addendum.md`: canonical-panel and wave-2
  applied analyses.
- `preregistration_bonander_addendum.md`: Florida stand-your-ground certification asymmetry.
- `preregistration_bi63_wisconsin_addendum.md`, `preregistration_fbi_rtc_nebraska_addendum.md`:
  disclosed-proxy and FBI Nebraska screens.

### `gates/`
- `gate_g0_g1_decision.md`: model validity and prior-art verdicts (Phase A).
- `gate_g2_decision.md`: prototype and enabling-formalization verdict (Phase B).
- `gate_g3_memo.md`: the decisive simulation gate: PIVOT to the restricted claim set, plus the
  C5 repair-and-confirm outcomes (Phase C).
- `gate_g4_memo.md`, `gate_g4_wave2_memo.md`: canonical-panel applied gate; novel-finding arms
  failed, C5 demoted.
- `gate_gBon_memo.md`, `gate_gBi63_memo.md`, `gate_gFBI_RTC_memo.md`: later applied screens and
  their dispositions.

### `applied/`
- `preprocessing_frozen.md`: frozen preprocessing and 21/21 trusted-result reproduction gates.
- `applied_study_closure_memo.md`: which applied paths are closed, why, and the only condition
  for reopening the applied gate.
- `options.md`: one-line decision note kept for provenance.
- `../Applied_Study/`: local-only replication archives (gitignored).

### `evidence/`
- `evidence_register.md`: source ledger with anchors, verification levels, and positioning
  obligations.
- `priorart_deepread_memo.md`: guarantee-section deep reads of the closest prior art.
- `novelty_search_log.md`: collision-search queries, dates, and hit counts.
- `pilot_cost_report.md`: LOCAL vs Colab classification from the pilot cost model.

### `literature/`
Candidate-panel and literature-search notes produced during the applied hunt (multiple model
assistants, kept for audit). These informed, but do not carry, the go decisions.

### `colab/`
- `README_COLAB.md`: fleet manual for the 47 self-contained Phase C notebooks, merge discipline,
  and determinism rules.

## Conventions

- Paths in the documents are repository-relative; the historical layout used `research/idea5/`
  and the references were rewritten during the 2026-09-18 reorganization.
- Result files (`results_*/`, `figures/`) are experimental outputs and are never edited by hand
  or gitignored; corrections are made by re-running the frozen drivers.
- Seeds are registered in `../config/seeds.yaml`; the row schema is
  `../config/results_schema.yaml`.
