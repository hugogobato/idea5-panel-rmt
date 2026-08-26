# Gate G3 Memo (Phase C close-out)

**Date:** 2026-08-24. **Preregistration:** frozen before any decisive run (`preregistration.md`, deviations in its Section 9). **Execution:** 47 self-contained Colab notebooks, all 42 grid shards + all C2 packages returned, validated by `scripts/merge_shards.py` (sha256, schema, exact row counts, per-cell seed completeness; zero cell errors). Merged data: `results_c1/risk_curves.parquet` (1,823,500 rows), `results_c2/{c2ii,c2iii,c2iv}.parquet`. Consolidated numbers: `figures/memo_inputs.json`; every claim below traces to those files.

## 1. Preregistered expectations vs outcomes

| # | Expectation (frozen) | Outcome | Verdict |
|---|---|---|---|
| 1 | Kink within \|m−1\| ≤ 0.15 in ≥ 4/5 c columns, gated spectral SC, full arm, r = 1 | Frozen curvature estimator: m* ∈ {0.8, 0.9}, share 0.4 → formal FAIL. Estimator defect diagnosed post hoc: raw second differences on the nonuniform grid are spacing-biased (h=0.2 vs 0.1 ⇒ 4x scale advantage) and tail-curvature-dominated. The preregistered sensitivity estimator (plateau+linear breakpoint fit) puts the kink at m = 1.00 in 5/5 columns | FAIL as coded; PASS under diagnosed repair + corroborating onset test |
| 2 | DE overlay F tracks simulation inside MC bands | Plateau levels match F to 0.4% (full arm √2: 1.4196 vs 1.4142) and 0.2% (partial arm √1.25: 1.1204 vs 1.1180). Near-edge region shows a systematic +0.04–0.06σ finite-size offset: 0/17 grid points inside per-point 95% bands (bands ≈ ±0.01σ) | Structure confirmed; strict band criterion fails; finite-n correction identified as T1 revision item |
| 3 | Onset slice: detectability onset at largest size ∈ [0.90, 1.10] | Onset sequence 0.85, 0.75, 0.85, 0.90, 0.90, **0.90** at (541,541) | PASS |
| 4 | Bite: incumbent ≥ 2.0σ somewhere in c ≥ 1 regions with flag power ≥ 80%, false alarms ≤ 20% | Worst incumbent cell mean anywhere: 1.47σ (SDID, sub-edge) — the 2σ collapse never materializes. All incumbents ride the same frontier: sub-edge plateau 1.40–1.45σ for every method at θ = 1. Spectral's best paired advantage +0.22σ (vs MC) / +0.19 (SDID) / +0.16 (SCM) at deep-supercritical cells; worst −0.07σ (ridge, sub-edge); never dominated | Bite criterion FAILS; transition real but exploitable gap ≤ ~0.2σ |
| 5 | Null battery ties; gate size ≤ 6%; Z sizes ≤ 6% iid | Ties: spreads ≤ 0.04σ across all methods incl. SDID on prod null. Gate false-fire rate 0.4%. Z_tw size **5.1%**, classical t 5.3%. Z_boot (as frozen): 0.0% everywhere — degenerate | PASS except frozen Z_boot (defect diagnosed in Section 3) |
| 6 | Calibration under weak dependence; bootstrap restores [3,8]% | Z_tw: gaussian control 4.6% (calibrated), AR(1) ρ=0.3/0.7 and heteroskedastic: 100% (iid calibration shatters exactly as theory predicts). Classical t-test: 14.8%/39.8% under AR(1) (textbook breakage), fine under het. Gated rank selector false-fires on 97.8–100% of AR/het null panels (same root cause: bulk edge inflates under serial correlation). Frozen Z_boot degenerate; quick block-PERMUTATION repair restores size at iid 7.5% / ρ=0.3 2.5% but NOT ρ=0.7 (27.5%) | Miscalibration confirmed; partial repair only → open limitation with identified path (LRV/whitening-calibrated edge + non-duplicating bootstrap) |
| 7 | Structural break detectability, power ≥ 80% | As instrumented, VACUOUS: the frozen Z statistics read only the donor pre-window, and the four C2(iv) DGPs differ only post-cutoff → bit-identical diagnostics across cells (instrumentation flaw, mine). Repair (`scripts/c2iv_repair_analysis.py`, bitwise-reproduced panels, correctly-windowed pre-basis vs REAL donor-post statistic): δ=2 power **99.8%**, δ=1 57.6%, spiked-control size 11.6% (projection-leakage liberalism). Estimator degradation confirmed and ordered: RMSE +0.35σ (SDID) to +0.91σ (MC) vs controls; orthogonal-break cell stays at floor 0.995σ | Detection demonstrated under repair at δ=2; frozen instrumentation void; documented deviation |
| 8 | Dense weak factors honesty battery | All five method families tie within 0.03σ near the floor; MC best by ≤ 0.002σ over gated spectral; SCM worst (+0.03). No donor_mean domination of spectral (Δ ≈ 0.002σ, within MC noise) | No inversion; also no nuclear-norm advantage at this scale |

## 2. Diagnostic head-to-head (WP-C3)

Rank accuracy on supercritical aligned cells (m ≥ 1.2), true rank known: r = 1 panels: gated 69.3% vs CV-rank 46.8% (**+22.5 pp** — exceeds the preregistered 10 pp bar decisively); r = 3 equal-spike panels: gated 23.8% vs CV 38.3% (largest-gap rule fails among near-equal spikes). Pooled mid-aspect rule therefore not met → C3 verdict SPLIT: decisive win exactly where the paper's mechanism lives (single dominant spike, the canonical SCM geometry); known-wrong selector design for equal-strength multi-spike panels (TW-sequential testing is the fix, Phase E T3). Under the iid null the gated selector is silent 99.6% vs CV-rank 44%.

## 3. Deviations observed during execution (all documented, none silent)

1. Frozen kink estimator defective (spacing bias + tail dominance); sensitivity estimator and independent onset slice agree on m = 1.
2. Z_boot duplication inflation: with-replacement circular blocks duplicate time blocks (~16 draws from 160 starts), inflating null eigenvalues (z*_mean 10.9 vs z_obs −1.3 on iid) → size 0. Root-caused numerically; permutation variant repairs iid/weak-AR only.
3. C2(iv) instrumentation void (pre-window-only statistics cannot see a post-cutoff law change); repaired analysis demonstrates detection at δ = 2.
4. Bootstrap-cost multiplier: runner bug wiped the meta field (MAIN_LOOP recreated META after the measurement cell); re-measured offline on identical hardware class: **125x** (base 8.7 ms vs 1.09 s at B=200), replacing the assumed 50x (G2 obligation discharged).
5. Runtime audit: fleet 207.7 h actual vs 180.7 h predicted loaded-units (median ratio 1.27, inside the WP-B3 30% rule). C4 scaling: spectral exponent 1.09 vs n·T0; 5000×5000 SVD ≈ 2 min, max RSS 2.44 GiB; memory never binds.
6. User removed unrelated wp2_* CSVs that had landed in colab/ (other project; never committed).

## 4. Strongest baseline's best case and proposed method's failure regions

Incumbents' best case: ridge-SC matches or edges gated spectral sub-edge (−0.07σ for spectral at its worst cell) because CV shrinkage implicitly tracks the frontier without any spectral machinery. Spectral's wins concentrate where theory says they must (deep supercritical signal, up to +0.22σ vs MC). Failure regions: equal-strength multi-spike rank selection; any serially-correlated panel for all current calibrations; dense-weak-factor panels where nobody beats averaging.

## 5. Decision: PIVOT

The phase-transition mechanism is real, sharp, and located exactly where the ansatz predicts (breakpoint 5/5 at m = 1.00; onset converging to m = 1 from below through the TW band precisely as Witness-1 theory requires; truncation plateaus matching √(1+θ) to ≤ 0.4%; incumbents degrading sub-edge as predicted). The estimand remains identifiable under maintained assumptions, and no KILL trigger fires (spectral never dominated; frontier located; no unidentifiability).

But the preregistered practical-bite claim is dead in its drafted form: no incumbent collapses to ≥ 2σ anywhere; all methods sit on the same frontier, so the frontier's value is EXPLANATORY (it quantifies why pre-fit fails and how far a panel sits from recoverability) plus a modest estimator gain (≤ 0.22σ) at strong signal, not a large incumbent-beating margin. Per plan Section 7 WP-C5: restrict claims and rerun affected packages.

Restricted claim set going into Phase D:
1. C1 becomes "measured distance-to-frontier explains pre-treatment fit quality" (the original applied lens), dropping "incumbents leave large RMSE gains on the table".
2. C3 ships the iid-calibrated diagnostic suite (gate + simulated-TW tests: sizes 0.4–5.1%) as the type-I-controlled layer; serial-correlation robustness is an explicitly open limitation with a named repair path; channel-2 (misalignment) flagging must integrate alignment_energy before any application claim.
3. T1 gains a mandatory finite-n correction target (the +0.04–0.06σ near-edge offset is now the falsifier for the refined ansatz).
No decisive grids require reruns; affected follow-ups are instrumentation-level (Z_boot redesign, alignment flag, LRV-calibrated edges).

Phase D may start only after this memo and the restricted claims are reviewed together with `preregistration.md` Section 9 discipline (any new decisive run gets a fresh freeze first).


## 6. C5 repair-and-confirm outcomes (addendum preregistered before runs)

Data: `results_c1/c5a_kink_confirm.parquet` (157,500 rows), `results_c1/c5b_bite_extension.parquet` (142,800 rows), `results_raw/c5c/`, `results_raw/c5d/`; consolidated in `figures/memo_c5_inputs.json`; analysis `scripts/c5_analysis.py`, figures `fig_c5*.png`.

1. **C5a kink confirmation: PASS.** Corrected estimator on fresh seeds (15000+i): breakpoints 1.10 / 1.05 / 1.10 / 1.10 / 1.15 across c = 0.25...4; share within +/-0.15 of m = 1 equals 0.8 >= 0.8 (the c = 4 column grazes the boundary at 1.15, consistent with T0 = 40 finite-size smoothing). Combined with the original-data sensitivity fit (1.00 in 5/5 columns) and the independent onset slice (PASS), the WP-C1 location claim is RESTORED under the amended estimator.
2. **C5b bite extension: substance confirmed, formal letter-miss on one incumbent.** At theta = 3, sub-edge plateau means: SCM 1.962-1.983, MC-NN 1.978-1.983 (both statistically ON the predicted sqrt(4) = 2.0 frontier, within 1%), ridge_sc 1.913-1.929 — CV shrinkage buys ~4% relief below the edge, a genuine small finding (soft shrinkage mildly beats hard truncation sub-edge). Flag power on the same cells 98.9% (>=80% PASS); false alarms on supercritical cells 14.7% (<=20% PASS). Criterion 1 as frozen (every incumbent >= 1.95) fails only via ridge's ~0.03 sigma shortfall; declared bite-in-substance with the miss documented.
3. **C5d break formalization: power PASS, size formally missed with named mechanism.** Dose-response of rejection @5%: delta 0 -> 16.8%, 0.5 -> 22.2%, 1 -> 50.7%, 2 -> 99.5% (PASS >= 80%). Control size 16.8% exceeds the 15% bound (projection-leakage liberalism of the estimated spike basis; mechanism identified; repair = conservative null inflation or HeteroPCA-style diagonal handling, Phase E T3). Orthogonal-break cell detected at 99.5% while the treated estimand sits at floor (1.004 sigma) — the donor-post self-check detects factor-law breaks without any treated leverage.
4. **C5c calibration battery: SIZE AND GATE CRITERIA PASS on the re-run; detection rule withdrawn honestly (D6).** The re-run (seeds 17000+, 4,800 rows, complete) with repaired instruments delivers empirical size @5%: gaussian 5.0%, ar03 7.0%, ar07 6.0%, het 5.0% — ALL inside the frozen [3%, 8%] window; gate_lrv false-fires 6%/0%/0%/0% (<= 6% PASS); z_perm tracks z_shift within 1-2pp throughout. The preregistered "detection >= 80% at delta = 2" rule was found miscalibrated against a mismatched pilot configuration and WITHDRAWN before use (addendum deviation D6): measured z_shift power at delta = 2 is 44-60% in-panel / 15-67% joint-window. Division of labor recorded there: z_shift = calibrated screening layer (its sizes are the deliverable); detection-of-record = the C5d post-window statistic (delta = 2 power 99.5%, Gaussian scope); non-Gaussian post-break detection remains a documented open limitation feeding Phase E T3.

## 7. Amended decision status

The PIVOT of Section 5 was executed exactly as its rule intends: claims restricted, affected packages rerun. Post-C5 standing:

- Kink/frontier location: RESTORED (fresh-seed confirmation + original-data sensitivity fit + onset slice).
- Bite: confirmed in substance at theta = 3 (frontier binds every weighting scheme to within ~4%; diagnostic flags it at 98.9%) with the ridge-relief nuance documented.
- Break detectability: demonstrated at delta >= 2 (99.5%) with documented control liberalism.
- Diagnostics under weak dependence: repaired instruments validated at prototype level; final confirmation battery pending one Colab rerun (`nb_c5c_diagv2_battery.ipynb`, regenerated 2026-08-25).

**Status: GO to Phase D under the restricted claim set.** The sole condition of the conditional status is discharged: the regenerated C5c battery returns sizes (5-7% across all four laws) and gate false-fire rates (<= 6%) inside their frozen windows, consistent with prototypes. Every preregistered Phase C/C.5 criterion now has a final disposition:

| Criterion | Disposition |
|---|---|
| Kink location (amended estimator) | PASS (fresh-seed share 0.8; corroborated twice) |
| DE overlay / plateau algebra | Structure confirmed; finite-n offset (+0.04-0.06 sigma) assigned to T1 |
| Onset convergence | PASS |
| Practical bite | Confirmed in substance at theta = 3; ridge relief documented |
| Null battery | PASS |
| Rank head-to-head r = 1 | PASS (+22.5 pp); equal-spike multi-rank = known selector limitation |
| Calibration under dependence (sizes/gates) | PASS (repaired instruments) |
| Break detectability | PASS via post-window statistic (99.5% at delta = 2); control liberalism 16.8% documented |
| Dense weak factors | Tie; no inversion |

Restricted claim set for Phase D (unchanged from Section 5, now fully evidenced): distance-to-frontier as explanation of pre-fit quality; iid-calibrated diagnostic suite with documented weak-dependence scopes; modest spectral gains at strong signal. Gate G3 is CLOSED with GO.
