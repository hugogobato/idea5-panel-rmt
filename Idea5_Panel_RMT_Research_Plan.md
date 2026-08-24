# Research Plan: Phase Transitions in Panel Causal Inference (Idea 5)

**Project:** Spectral recoverability frontiers for synthetic control, matrix completion, and staggered-adoption estimators on spiked factor panels.

**Derived from:** `Random_Matrix_Research_Ideas.md`, Idea 5 (2026-08-23).

**Current classification:** Promising but unproven. Prior-art novelty verified at the abstract/search level (E2-E3); all empirical scores UNTESTED until Gate G3/G4.

**Date:** 2026-08-23

---

## 1. Executive verdict

| Item | Statement |
|---|---|
| Classification | Promising but unproven |
| Confidence | Medium-high on novelty, low-medium on mechanism (untested numerically) |
| Proposed contribution | A BBP-type spikiness threshold below which the treated unit's counterfactual trajectory is unrecoverable by *any* panel method, plus TW-calibrated rank selection and a spectral pre-trends test, turning "pre-treatment fit quality" from folklore into a measured distance-to-frontier |
| Real contribution | The threshold/diagnostic layer (Claims C1, C3) |
| Engine | Deterministic-equivalent risk machinery (Claim C2); simulation infrastructure |
| Application | Distance-to-frontier re-analysis of canonical SCM panels + staggered panels (Claim C5) |
| Decoration | Staggered-adoption block-spiked generalization (Claim C4): cut unless it becomes the strongest application |
| Strongest reason it could become a strong field paper | Zero-hit collision searches verified on arXiv (2026-08-23); the literature has norm-bound guarantees only, and every applied SCM user already worries about pre-fit, so a calibrated frontier diagnostic has immediate uptake |
| Strongest reason it could fail | Panel noise is serially correlated and heteroskedastic, so iid Marchenko-Pastur/Tracy-Widom calibration may break exactly where it matters; and the threshold may exist but carry no practical bite if incumbents already behave near-oracle everywhere |
| Next unresolved gate | G0+G1 (validity + novelty deep-read), bundled into Phase A below |
| Single cheapest decisive next action | Work package A2: deep-read the two closest papers' guarantee sections and run the remaining alternate-vocabulary searches |

---

## 2. Idea reconstruction and claim decomposition

**Scientific problem.** Panel causal estimators assume outcomes follow Y = L + tau*Z + E with L low-rank (interactive fixed effects). Existing guarantees (Abadie-style SCM, nuclear-norm matrix completion, SDID, rate-optimal ATT results) are operator/Frobenius norm bounds under strong identifiability conditions. None characterizes *how spiky* the factor structure must be, relative to the MP noise sea of the donor-pool Gram matrix, before the treated unit's counterfactual is determined by donors. Poor pre-treatment fit, the dominant practical failure mode of SCM, is hypothesized to be precisely the sub-threshold regime.

**Unit of analysis:** a single treated unit observed over T₀ pre-periods against n₀ untreated donor units, extended to cohorts under staggered adoption.

**Observable data:** an n₀ x T₀ pre-period outcome matrix M = UVᵀ + E; donor outcomes continue into post-periods; the treated unit's post-period counterfactual is the target.

**Target estimand:** the treated unit's untreated potential outcome trajectory in post-periods, aggregated to ATT where relevant.

**Method:** (i) derive the BBP-type recoverability threshold as a functional of spike strengths, aspect ratio c = n₀/T₀, and the treated unit's leverage/alignment with recovered factor directions; (ii) deterministic-equivalent risks for SCM variants; (iii) Tracy-Widom/Jacobi-calibrated rank selection and pre-trends tests; (iv) release as software.

**Intended audience:** applied econometrics (program evaluation), panel causality methods literature, JASA/AEJ/Econometrics Journal / JRSS-B tier.

### Claim ledger

| ID | Claim | Type | Feasibility | Novelty | Importance | Evidence status | Keep / cut / pivot |
|---|---|---|---|---|---|---|---|
| C1 | Unrecoverability threshold: below a computable spikiness frontier no estimator attains o(1) counterfactual risk; matching construction above it | Contribution | Medium-high (BBP + spiked-MC machinery exists) | High (zero hits verified on arXiv) | High (explains the #1 practical failure mode) | Untested; ansatz not yet simulated | Keep; load-bearing |
| C2 | Exact limiting risk of constrained SC, SC+intercept, nuclear-norm MC, SDID via deterministic equivalents | Engine (partially contribution) | Medium for ridge/unconstrained; hard for simplex constraints | Medium-high | Medium (enables apples-to-apples comparison) | Untested | Keep ridge family first; simplex deferred |
| C3 | Calibrated rank selector (type-I control) and spectral pre-trends test valid under serial correlation | Contribution | Medium (Onatski 2010 template exists; serial correlation is the gap) | High | High (directly addresses known weakness of classical pre-trend testing) | Untested | Keep; co-load-bearing |
| C4 | Block-spiked staggered-adoption extension | Application/enabling | Medium | Medium | Medium (leverages DiD-BCF/CUPED infrastructure) | Untested | Cut unless Phase D shows it is the best application |
| C5 | Distance-to-frontier empirical lens on California smoking, German reunification (+ one staggered panel) | Application | High | Medium-high (new empirical reporting standard) | High for adoption | Untested | Keep |

**Load-bearing contribution:** C1 + C3. **Load-bearing assumption:** the post-period factor space is spanned by pre-period factor directions (no structural break in V), and the treated unit's loading has nonzero projection onto the spiked space. If both fail, the target is genuinely unidentifiable and the theory must say so explicitly rather than hide it.

**Strongest simple baseline:** plain Abadie SCM with informal pre-treatment-fit inspection, and nuclear-norm matrix completion with cross-validated rank (Athey-Bayati-Doudchenko-Imbens-Khosravi). Both must appear in the first decisive experiment.

**Cheapest decisive experiment:** sweep spike strength across the conjectured threshold in a Gaussian factor panel; plot counterfactual RMSE of all estimators against the predicted frontier. A sharp transition aligned with the prediction validates the ansatz; a smooth curve with no exploitable boundary kills the phase-transition story.

**Hardest credible referee objection:** "MP/TW laws require iid noise; real panels have serial correlation, heteroskedasticity, and time-varying loadings, so your calibrated tests miscalibrate exactly in the applications you advertise; moreover Athey et al. (JASA 2021) already unified these estimators, so the remaining delta looks like new notation for old objects." Answer required by end of Phase C/D: size control under weak dependence (or honest bootstrap fallback), plus an application finding incumbents cannot see.

---

## 3. Fatal-flaw certificate (Gate G0)

Checks performed at planning time (symbolic reasoning level; numerical witnesses scheduled in Phase A):

| Check | Result | Detail |
|---|---|---|
| Mathematical typing | Pass | M = UVᵀ + E with U in R^{n0 x r}, V in R^{T0 x r}, E iid-scale entries; estimand is a held-out row-block of UVᵀ; well-defined |
| Non-vacuity | Pass, with witness scheduled | r = 0 (pure noise) and r >= min(n0,T0) degenerate cases must bracket behavior; numeric check in A1 |
| Identification | Conditional pass | Counterfactual identifiable iff (i) post-period factor values are learnable from donor rows (donors observe the same V) and (ii) treated loading u* has nonzero projection onto spiked directions. Structural breaks in V across the cutoff are a genuine nonidentifiability region; the plan treats this as a stated assumption with a diagnostic, not a hidden condition |
| Existence/well-definedness of method | Pass | All six planned estimators (oracle, donor-mean, simplex SCM, ridge-SC, hard-threshold spectral SC, nuclear-norm MC, SDID) exist and are implementable |
| Information leakage | Watch item | Rank selection and threshold tuning must use pre-period data only; leakage rules written into WP B2 |
| Method-target mismatch | None detected | Spectral tools match the low-rank-plus-noise geometry |
| Smallest counterexamples | Scheduled | Two witnesses to build in A1: (a) sub-threshold treated unit indistinguishable from noise; (b) visible spike contributing zero bias (misaligned u*) |
| Vacuity-by-equivalence | Partially cleared | Athey et al. (2021) show SCM/MC/unconfoundedness optimize related objectives; this unifies *estimators*, not their *spectral limits or failure boundaries*. Remaining delta is the frontier + calibration layer. Deep-read in A2 to confirm |

**G0 verdict: CONDITIONAL PASS.** No material defect found at planning resolution. Two named witnesses must be produced in Phase A before any grid investment.

---

## 4. Verified prior art and nearest-neighbor map (Gate G1)

Searches executed 2026-08-23 via the arXiv API:

| Query | Hits |
|---|---|
| `"synthetic control" AND "eigenvalue"` | 0 |
| `"panel data" AND "Marchenko-Pastur"` | 0 |
| `"synthetic control" AND "random matrix"` | 0 |
| `"synthetic control" AND ("BBP threshold" OR "Tracy-Widom")` | 0 |

Google Scholar and RePEc/SMR remain unchecked (bounded task A3); arXiv zero-hits alone cannot clear G1.

### Nearest neighbors

| Source (status) | Same problem? | Same target? | Same method? | Same evidence? | Remaining gap | Direct-hit risk |
|---|---|---|---|---|---|---|
| Athey, Bayati, Doudchenko, Imbens, Khosravi, "Matrix Completion Methods for Causal Panel Data Models," JASA 2021 (verified, arXiv:1710.10251, DOI 10.1080/01621459.2021.1891924) | Yes | Counterfactual imputation in panels | Nuclear-norm MC; unifies SCM/MC objectives conceptually | Simulations + theory under missingness | No spectral thresholds, no RMT limits, no calibrated rank/pre-trends tests | Medium: they own the unification story; we must position as the *quantitative frontier* layer |
| Agarwal, Dahleh, Shah, Shen, "Causal Matrix Completion," arXiv:2109.15154 (verified abstract) | Adjacent (MNAR missingness) | Entry-wise potential outcomes | SNN; max-norm bounds | Finite-sample consistency/asymptotic normality | No BBP-type frontier, no spectral diagnostics | Low-medium |
| Farias, Li, Peng, arXiv:2106.02780 (verified abstract) | Yes | Rate-optimal ATT under general intervention patterns | Regularized estimation | Norm rates | Rates are not sharp constants; no detectability frontier; no tests | Medium |
| Agarwal, Han, Saha, Syrgkanis, Yoon, "Synthetic Blips," arXiv:2210.11003 (verified abstract) | Adjacent (dynamic treatments) | Unit-specific dynamic effects | Backward-induction weighting | Identification + algorithms | No spectral analysis | Low |
| Moon-Weidner line on interactive fixed effects / PC estimators in dynamic panels (UNVERIFIED, to inspect in A3) | Yes (estimation) | Factor estimates, ATT | Principal components | Consistency/rates under many factors | Not RMT-exact; no individual-unit recoverability question | Unknown until inspected; flagged |
| Spiked matrix completion / sparse PCA minimax literature (Birnbaum et al.; Koltchinskii-Lounici; Perry-Wein style) (UNVERIFIED at detail level) | Adjacent (statistics) | Low-rank recovery under missingness | Spectral methods + lower bounds | Minimax rates | Missing-data spectrum distortion studied, but never mapped to panel causal estimands or SCM practice | Low-medium: supplies our proof tools rather than competing |
| Chen-Ma alignability tests, arXiv:2511.21074 (cited in dossier, UNVERIFIED) | Adjacent (two-sample) | Distribution alignment | Principal variances | Asymptotic tests | Not panel/counterfactual | Low |
| Onatski (2010 REStat) factor-count testing; Ait-Sahalia-Xiu (2019) | Adjacent | Number of factors | Eigenvalue ratios | Asymptotic null distributions | Never used for treated-unit recoverability | Complementary, must cite |

**What is novel after verification:** the *question* (per-treated-unit recoverability frontier as a function of spikiness and leverage), the *calibration layer* (rank selection and pre-trends tests with controlled size in the q/n -> c regime), and the *empirical lens* (distance-to-frontier as reported quantity). What is NOT novel: low-rank panel estimators themselves, unification arguments, and generic spiked-model tools.

**Strongest incumbent:** Athey et al. nuclear-norm MC with CV rank choice, and plain SCM with pre-fit inspection.

**G1 verdict: CONDITIONAL GO.** Proceed to Phase A deep-reads; a direct hit in the remaining vocabularies downgrades to PIVOT or INCREMENTAL-ONLY immediately.

---

## 5. Impact thesis and skeptical-referee test

1. **Why the problem matters.** SCM-family methods are the default for major policy evaluations; their most common diagnostic (pre-treatment fit) is currently eyeballed, with no statistical meaning.
2. **What changes if this works.** Practitioners get a calibrated number: "this treated unit sits 0.3x below the recoverability frontier; no weighting scheme can fix it," plus honest rank selection and a pre-trends test that survives serial correlation.
3. **Who uses/cites it.** Applied microeconomists running SCM/SDID; methods authors needing a baseline theory for why fit fails; the DiD-BCF/DiD-CUPED projects in this repository.
4. **Why simpler incumbents are insufficient.** Pre-fit RMSE conflates noise level, panel size, and factor strength into one ad hoc glance; CV rank selection has no type-I control; classical pre-trend tests lose size under serial correlation.
5. **Why not merely a combination.** The threshold is a new *limit object* (not a re-weighted incumbent), and it produces predictions incumbents do not make (e.g., misaligned-but-visible spikes contribute nothing; invisible-below-threshold units are unrecoverable).
6. **Most damaging referee paragraph.** As in Section 2: iid-noise fragility + "Athey et al. already unified this."
7. **Evidence needed to answer it.** Size plots of the spectral tests under AR(1)/heteroskedastic errors (Phase C), and an application finding invisible to standard pipelines (Phase D).

### Impact dimension scores (initial triage; empirical scores UNTESTED)

| Dimension | Score | Note |
|---|---|---|
| Problem importance | 3 | Core applied-econometrics methodology |
| Novelty after prior art | 2 | Verified zero-hit at interface; adjacent literatures strong |
| Mechanism or insight | 1 (UNTESTED) | BBP physics plausible; not yet simulated |
| Empirical advantage | 0 (UNTESTED) | Pending Gate G3 |
| Applied value | 0 (UNTESTED) | Pending Gate G4 |
| Generality | 2 | Whole SCM/MC/SDID family |
| Credibility | 2 | Standard assumptions, clearly stated limits |
| Paper coherence | 2 | One story: how much signal does SCM need? |

---

## 6. Dependency graph and gate map

```text
G0 validity + G1 novelty  (Phase A)
    -> G2 enabling formalization + prototype ladder  (Phase B)
        -> G3 simulation-first falsification  (Phase C)
            -> G4 applied study  (Phase D)
                -> G5 evidence-earned theory + G6 submission case  (Phase E)
```

Everything downstream of a pending gate is DORMANT UNTIL that gate. Read-only data feasibility checks (file existence, schema) may run during Phase B; no outcome inspection before Gate G3 passes.

---

## 7. Phase-by-phase execution program

### Phase A: Validity and novelty preflight (Gates G0 + G1)  [ACTIVE]

**Purpose.** Prove the idea is worth two weeks before touching grids: produce the two mathematical witnesses, verify no direct hit hides in adjacent vocabulary, and pin the exact model assumptions.

**Prerequisites.** None beyond the dossier and this plan.

#### WP-A1: Model statement + numerical witnesses
- Status: ACTIVE
- Gate served: G0
- Objective: Fix the formal model, estimand, and assumption ledger; build two toy witnesses.
- Why this changes a decision: a nonidentifiable target or vacuous threshold kills the project for the cost of one notebook.
- Inputs: Section 2 of this plan; dossier Idea 5.
- Actions:
  1. Write `model_card.md`: data regime (n0, T0 -> infinity, c = n0/T0 fixed), model Y = L + tau Z + E, L = UVᵀ with r fixed spikes of strengths l_1 >= ... >= l_r > sqrt(c)-scale MP edge, E with variance sigma^2; estimand = treated-row counterfactual in post-periods; assumption ledger including the no-structural-break-in-V assumption and its diagnostic.
  2. Witness 1 (sub-threshold invisibility): generate a panel where the treated unit's loading component along one factor falls below the BBP edge; verify empirically that its estimated coefficient distribution is asymptotically indistinguishable from pure-noise regression.
  3. Witness 2 (visible-but-useless spike): construct a strong spike orthogonal to the treated unit's loading; verify large outlier eigenvalue coexists with zero reduction in counterfactual error.
- Outputs: `research/idea5/model_card.md`; `research/idea5/code/witness_subthreshold.ipynb`; `witness_misalignment.ipynb`.
- Verification (mechanical): both notebooks run clean top-to-bottom in a fresh kernel.
- Verification (scientific): Witness 1 shows overlapping coefficient/noise distributions within Monte Carlo bands; Witness 2 shows outlier eigenvalue present while RMSE unchanged vs. spike-free panel within 1 MC sd.
- Pass rule: both witnesses behave as described.
- Fail rule: either witness contradicts the mechanism (e.g., sub-threshold units are recoverable by plain OLS weights).
- Gate consequence: fail => PIVOT (re-scope the threshold claim) or KILL; pass => proceed to A2.
- Compute: laptop-scale, minutes.
- Likely trap: accidentally giving Witness 1 a spike just above the edge; scan a fine grid of spike strengths around the conjectured edge.
- Recovery: rerun with denser grid and larger T0.

#### WP-A2: Deep-read of closest prior art
- Status: ACTIVE
- Gate served: G1
- Objective: Inspect guarantee sections, not abstracts, of the three closest sources.
- Inputs: arXiv PDFs 1710.10251, 2106.02780, 2109.15154.
- Actions:
  1. Extract from each: exact estimand, assumptions (noise model, missingness, factor scaling), guarantee form (norm/rate/constant), inference offered, limitations admitted.
  2. Fill the evidence-register rows with page/theorem anchors and verification level E3.
  3. Write a half-page memo: does any of them contain (a) a per-unit recoverability threshold, (b) TW-calibrated diagnostics, (c) exact high-dimensional limits? Expected answer: no; record otherwise if found.
- Outputs: `research/idea5/evidence_register.md` (updated), `priorart_deepread_memo.md`.
- Verification: every row carries an anchor (section/theorem/page).
- Pass rule: no direct hit on (a)-(c).
- Fail rule: direct hit found.
- Gate consequence: fail => KILL original framing, convene pivot decision.
- Dependencies: none; can run parallel with A1.
- Compute: reading only.

#### WP-A3: Residual novelty searches + unverified-citation audit
- Status: ACTIVE
- Gate served: G1
- Objective: Close the remaining search surface.
- Actions:
  1. Query families (Google Scholar, RePEc, arXiv full-text): "factor augmented" + "prediction" + "many units"; "principal components" + "panel" + "consistency" (Moon-Weidner line); "spiked covariance" + "missing data"; "denoising" + "spiked"; "synthetic control" + "high dimensional asymptotics"; "pre-treatment fit" + "test".
  2. Verify the dossier's uncited/2026-ID citations (2603.24833, 2605.30319) and the Agarwal-Dahleh-Sarkar 2019 reference; replace or annotate anything unresolved as UNVERIFIED.
  3. Record Moon-Weidner conclusions in the register.
- Outputs: updated `evidence_register.md`; `novelty_search_log.md` with queries, dates, hit counts.
- Verification: log contains query strings and counts; no UNVERIFIED entry supports any go decision.
- Pass/fail: same as A2.
- Gate consequence: completes G1 decision.

**Phase A give-up rules.**
- KILL: a verified direct hit (someone already derives the per-unit recoverability threshold or TW-calibrated SCM diagnostics); or WP-A1 shows the counterfactual estimand is unidentifiable even under the favorable regime beyond the acknowledged structural-break caveat.
- PIVOT: witnesses hold but reveal the interesting object is something else (e.g., transition is in coverage, not RMSE); rewrite C1 accordingly and rerun G0.
- CONDITIONAL GO: proceed to Phase B only with both witnesses passing and no direct hit.

**Estimated effort:** ~1 week calendar. Compute trivial.

---

### Phase B: Minimum enabling formalization + prototype correctness ladder (Gate G2)  [ACTIVE after Phase A]

**Purpose.** Build the smallest honest end-to-end pipeline and the *heuristic* threshold formula (tagged conjecture). No proofs in this phase beyond what makes the code meaningful.

**Prerequisites.** Phase A passed.

#### WP-B1: Conjectured frontier formula (paper-and-pencil + symbolic check only)
- Status: BLOCKED until Phase A
- Objective: State the deterministic-equivalent ansatz for the counterfactual risk of ridge/hard-threshold SC and the location of the recoverability boundary as an explicit function F(l_j, c, sigma^2, alignment(u*, V)).
- Actions:
  1. Derive the regression-of-treated-row-on-top-r-right-singular-vectors risk decomposition: bias from truncating spikes below the BBP edge + variance sigma^2 * c-trace terms.
  2. Special cases to satisfy algebraically: r=0; l_j -> infinity (risk -> interpolation floor); u* exactly orthogonal to factor space (risk equals noise floor regardless of spikes).
  3. Tag the output CONJECTURE. Map each ingredient to its source result (BBP outlier locations; Benaych-Georges-Nadakuditi overlaps; Dobriban-Wager ridge risk) in `theory_targets.md` stub.
- Outputs: `frontier_ansatz.md` with formula F and the three special-case derivations.
- Verification: special-case reductions hold symbolically.
- Pass rule: consistent special cases; formula predicts monotone risk decreasing in spike strength with a kink near the BBP edge.
- Fail rule: internal contradiction in special cases.
- Gate consequence: fail => revise once; second failure => downgrade C1 from sharp threshold to empirical transition map and note scope change at G2 review.
- Compute: none.

#### WP-B2: Reference implementation library
- Status: BLOCKED until WP-B1
- Objective: One Python module implementing DGPs and all estimators with unit tests.
- Actions:
  1. `dgps.py`: factor-panel generator with knobs (n0, T0, r, {l_j}, alignment profile, noise law: gaussian / ar(1) / heteroskedastic; optional structural break).
  2. Estimators: oracle (uses true L), donor-mean (trivial baseline), simplex-constrained SCM (Abadie), ridge-SC, hard-threshold spectral SC, nuclear-norm MC with CV rank, SDID (simple implementation or vetted package port), plus the diagnostic suite (scree, eigenvalue-ratio rank selector, TW statistic).
  3. Unit tests: shapes; r=0 case returns noise-floor predictions; infinite-spike case returns near-oracle; leakage guard asserting rank selection touches pre-period data only; determinism under fixed seeds.
  4. Smoke run: 20 reps, one cell, all methods, < 1 minute.
- Outputs: `research/idea5/code/scm_frontier/` package; `tests/test_estimators.py`; CI-style pytest run log.
- Verification (mechanical): `pytest` green.
- Verification (scientific): oracle beats everything everywhere by construction; trivial baseline never beats tuned methods by more than MC noise in favorable cells.
- Pass rule: all above.
- Fail rule: any estimator unstable in its home turf (e.g., simplex solver non-convergence in >1% of runs).
- Gate consequence: instability => PIVOT within Phase B to ridge-only scope (documented), since simplex DEs were already flagged as the risky branch.
- Compute: minutes on laptop.

#### WP-B3: One-seed pilot + cost model
- Status: BLOCKED until WP-B2
- Objective: Measure wall time/RAM per replication; decide what stays local vs. goes to Colab.
- Actions:
  1. Run the decisive-cell pilot (see Section 8) with 10 seeds on the local machine; record wall time, peak RSS per rep.
  2. Fit linear cost model; classify each planned experiment in Section 8 as LOCAL (predicted < 2 h total AND < 4 GB peak) or COLAB (otherwise).
  3. Check current machine load before choosing worker counts; cap workers at 8 physical cores, set thread limits (OMP_NUM_THREADS=1 inside workers) to avoid nested parallelism, leave headroom for other jobs.
- Outputs: `pilot_cost_report.md` with the LOCAL/COLAB classification table.
- Verification: numbers logged, not guessed.
- Pass rule: cost model explains pilot within 30%.
- Gate consequence: feeds Phase C sharding decisions; no scientific decision.

**Phase B give-up rules.**
- PIVOT: simplex SCM cannot be made stable (restrict scope to ridge/hard-threshold family; document as deliberate scope change, revisit simplex via approximation later).
- KILL: the ansatz contradicts its own algebra twice after revision (downgrade path already defined) OR the reference implementation cannot reproduce the oracle ordering even in noiseless DGPs (implementation logic broken beyond repair).
- GO: green tests, sane pilot, formula standing.

**Estimated effort:** ~1 week. Compute: < 1 h local.

---

### Phase C: Simulation-first falsification (Gate G3)  [DORMANT UNTIL GATE B]

**Purpose.** Decide whether the phase-transition mechanism is real, sharp, and exploitable, and whether the diagnostics calibrate. This is the decisive gate of the whole project.

**Prerequisites.** Phase B passed; preregistration below frozen BEFORE running decisive cells.

#### Preregistration (freeze in `preregistration.md` before any decisive run)
- Primary metric: post-period counterfactual RMSE normalized by sigma, averaged over treated rows, per (c, spike-profile, alignment) cell.
- Secondary: ATT bias; rank-selector accuracy; empirical size of TW pre-trends test at nominal 5%; coverage of spectral CI.
- Practical-effect threshold: the frontier "has bite" if some incumbent suffers RMSE >= 2x oracle-normalized floor in at least one substantive region while the spectral diagnostic flags it with power >= 80% at 5% size.
- Seeds: 500 replications per cell, seeds 10000+i, identical across methods; equal tuning budget (same inner CV folds) for all adaptive methods.
- Failure rule (falsifier): if the empirical risk-vs-spike curve shows no kink aligned with F within finite-size tolerance (kink located within +/- 15% of predicted edge in >= 80% of c-grid points), the ansatz is falsified.

#### WP-C1: Decisive threshold sweep
- Objective: risk curves vs. predicted frontier across the crossover grid.
- Grid: c in {0.25, 0.5, 1, 2, 4}; spike strength multiplier s in linspace(0.2, 3, 15) around the BBP edge; alignment in {full, partial, orthogonal}; r in {1, 3}.
- Methods: all seven from WP-B2. Baselines present from day one: plain SCM and CV-rank nuclear-norm MC.
- Outputs: `results_c1/` parquet shards + merged `risk_curves.parquet`, figure `fig_threshold_kink.png`.
- Verification: all seeds present, schema matches `results_schema.yaml`.
- Scientific pass: preregistered kink criterion met; incumbents degrade sub-edge as predicted.
- Fail: smooth curves or kink far from prediction.
- Gate consequence: fail triggers the ansatz-revision loop (max two documented revisions, each naming the corrected ingredient); third failure => KILL C1.

#### WP-C2: Null, baseline-favorable, and misspecification cells
- Objective: honesty battery.
- Cells: (i) r=0 null: all estimators statistically tied; spectral tests size <= 6% nominal 5%; (ii) dense weak-factor DGP favoring nuclear-norm MC (baseline-favorable); (iii) AR(1) errors rho in {0.3, 0.7} and heteroskedastic noise: report calibration drift of TW tests and whether block-bootstrap correction restores size; (iv) structural-break DGP: confirm all methods fail and the diagnostic fires (this cell demonstrates the assumed-away region is detectable).
- Outputs: `results_c2/`, `fig_size_calibration.png`, `fig_baseline_favorable.png`.
- Pass: null ties hold; size restored within [3%, 8%] under bootstrap correction or documented as open limitation; structural-break cell shows diagnostic power >= 80%.
- Fail: miscalibration unfixable under weak dependence AND diagnostic fires on benign panels (false alarms dominate).
- Gate consequence: miscalibration-only failure => PIVOT diagnostics toward bootstrap-native versions (scope change recorded).

#### WP-C3: Diagnostic head-to-head
- Objective: spectral rank selector + pre-trends test vs. incumbents (CV rank selection; classical linear pre-trend t-test under serial correlation).
- Metric hierarchy fixed in preregistration; ablation: diagnostic with and without TW correction isolates the claimed mechanism.
- Pass: spectral selector dominates CV-rank accuracy by >= 10 percentage points in at least the mid-aspect regimes, and pre-trends test strictly out-sizes the classical test under rho >= 0.3.
- Fail: no regime where the diagnostic wins => INCREMENTAL-ONLY pressure on C3.

#### WP-C4: Scaling study
- Objective: runtime/memory scaling to n0 = T0 = 2000-5000; identifies which production sweeps need Colab.
- Outputs: `scaling_report.md`.

#### WP-C5: Gate memo and decision
- Standalone memo `gate_g3_memo.md` containing preregistered expectations, result tables/plots, deviations, strongest baseline's best case, proposed method's failure regions, ablation evidence, practical effect sizes, and ONE decision.
- Decision rules:
  - GO to Phase D: threshold has bite (preregistered criterion), diagnostics calibrate or have a working fallback, and at least one incumbent-beating region exists with practical size.
  - PIVOT: transition real but bite restricted (e.g., only high-c or only ridge variants); restrict claims and rerun affected C-packages.
  - INCREMENTAL-ONLY: threshold exists but incumbents already near-oracle everywhere, or diagnostics never beat CV/pre-fit heuristics; if user bar excludes incremental work, terminate.
  - KILL: no sharp transition (story dead), or methods worse than benchmark models across all substantive regimes (the specific user-stated give-up condition), or diagnostics miscalibrated with no repair.

**Phase C give-up rules summary (explicit).**
1. KILL if the proposed spectral methods are worse than benchmark models (plain SCM, CV-rank MC, SDID) in every substantive regime, not just contrived cells.
2. KILL if the wanted object (recoverability frontier) turns out empirically nonexistent or unlocatable (no kink, no exploitable boundary).
3. KILL if the estimand proves unidentifiable in the favorable regime (structural-break-like failure even under maintained assumptions).
4. PIVOT/INCREMENTAL-ONLY per WP-C5 rules above.

**Estimated effort:** 2-3 weeks. Compute: decisive grid roughly 5 c-values x 15 spikes x 3 alignments x 2 r-values x 500 reps x 7 methods; per-rep SVDs on matrices up to 2000 x 500 (~seconds). Pilot will refine, but expect several LOCAL experiments in the 1-2 h range (parallelize over 8 workers) and the AR(1)+bootstrap battery and the 2000+ scaling cells flagged COLAB (each > 2 h projected). Sharding per Section 12.

---

### Phase D: Applied study (Gate G4)  [DORMANT UNTIL GATE C]

**Purpose.** Show the frontier lens reveals something incumbents cannot see on real panels.

**Prerequisites.** G3 GO. Early read-only feasibility (data existence/schema) permitted since Phase B.

#### WP-D1: Data acquisition and trusted-result reproduction
- Panels: California smoking (Prop 99), German reunification; optionally Basque conflict; one staggered-adoption dataset aligned with the DiD-BCF project if C4 survived.
- Actions: ingest, clean, reproduce the canonical published point estimates with a standard pipeline within tolerance; freeze preprocessing (`preprocessing_frozen.md`) before comparative runs.
- Pass: reproduction matches published estimates to reporting precision.
- Fail: cannot reproduce trusted benchmark => fix pipeline before any novel analysis; unresolved mismatch => treat as bug, not discovery.

#### WP-D2: Distance-to-frontier analysis
- Actions: for each dataset/unit, compute estimated spike spectrum, treated-unit alignment, and distance d to the predicted frontier; report d with uncertainty (bootstrap over donors/time blocks).
- Predeclared novel finding: at least one canonical application sits in a quantifiably fragile region (d < 1) where the frontier explains observed pre-fit weakness that informal inspection could not quantify, or conversely a case where poor-looking pre-fit is provably frontier-safe.
- Placebos/negative controls: in-space placebos (untreated units assigned pseudo-treatment) must yield d-distributions consistent with their healthy pre-fit; placebo treatment dates must not trigger the pre-trends alarm at nominal rates.
- Comparison: strongest incumbent analysis = published SCM/SDID specifications plus CV-rank MC.
- Pass: novel finding survives placebos and sensitivity (alternative rank selectors, trimmed windows); changes interpretation of at least one application.
- Fail: every canonical panel sits comfortably far from the frontier AND the diagnostic adds nothing beyond existing pre-fit inspection => applied value = 0.

**Phase D give-up rules.**
1. KILL if identification fails on real data (diagnostic fires on the canonical positive controls, i.e., flags true effects as unrecoverable artifacts) with no repair.
2. KILL if no useful finding: frontier distances are uniformly uninformative and placebos cannot separate treated from untreated units.
3. PIVOT if the lens works only for a different class than advertised (e.g., staggered panels rather than single-treated-unit SCM); retarget the paper accordingly and rerun G3-affected cells.
4. INCREMENTAL-ONLY if findings are purely descriptive additions to known case narratives.

**Estimated effort:** 3-4 weeks including package polish. Compute: light (< 1 h per panel); bootstrap batteries may hit the COLAB trigger.

---

### Phase E: Evidence-earned theory + paper consolidation (Gates G5 + G6)  [DORMANT UNTIL GATE D]

**Purpose.** Formalize only what surviving evidence makes load-bearing, then assemble the paper.

#### Theory target table (all conditional on G3/G4 outcomes)

| Target | Why it matters (evidence link) | Sketch | Tag | Source result to lean on | Adaptation gap | Numerical falsifier | Stop rule |
|---|---|---|---|---|---|---|---|
| T1 Frontier formula (C1 upper side) | Explains the kink WP-C1 demonstrated | Ridge-SC counterfactual risk -> deterministic equivalent; boundary where risk stops improving | Adaptation | Dobriban-Wager (2018) ridge risk; Benaych-Georges-Nadakuditi (2011/2012) overlaps; BBP (2005) outliers | Row-targeted (single unit) risk instead of Frobenius; tractable resolvent variant | Overlay DE curve on WP-C1 simulations to plotting accuracy | Two failed revision cycles => demote to conjecture with empirical support |
| T2 Matching lower bound (C1 lower side) | Turns folklore into theorem; the paper's spine | Le Cam two-point between treated-loading configurations separated below noise floor | Adaptation->Conjecture | Spiked matrix completion minimax templates (Koltchinskii-Lounici line; Birnbaum et al.); missing-noise spiked models | Missing pattern is structured (block row), not MCAR | Sub-threshold witness (WP-A1) at scale: no estimator crosses floor | If resisted after bounded effort: publish GO-WRITE version with empirical frontier only |
| T3 TW/Jacobi calibration of diagnostics | Validates C3's type-I claims | Largest-eigenvalue null after whitening; universality under weak dependence | Direct (iid) / Conjecture (AR) | Johnstone (2001, 2008); Onatski (2010) ratio statistics; Bao-Pan-Zhou universality | Serial correlation breaks exact nulls | WP-C2(iii) size plots | Bootstrap fallback stands as the shipped method; theory limited to iid case |
| T4 SDID deterministic equivalents | Apples-to-apples zoo comparison | Resolvent calculus for SDID's two-weight system | Conjecture | Arkhangelsky et al. (2021) definition; Dobriban-Wager machinery | Coupled weight systems lack closed forms | DE vs simulation overlay | Cut T4 entirely if it resists; it is decoration unless reviewers demand it |

Each target maps to a verified source entry in `evidence_register.md`; no target proceeds from an E0/E1 source.

#### WP-E1: Proof packages for T1-T3 in priority order (one package per target, each with special-case reductions and a counterexample search before the long proof).
#### WP-E2: Software release (`spectral-frontier` R/Python package): diagnostic suite + vignettes reproducing every paper figure.
#### WP-E3: Skeptical-referee pass + reproducibility audit (rerun-from-clean-env of every figure/table).
#### WP-E4: Venue stress test and submission case (target tier reasoned from current comparably-scoped papers: JASA/Econometrics J./JRSS-B for full program; arXiv + applied seminar circuit for the GO-WRITE variant).

**Phase E give-up rules.**
1. T2 resisted after bounded effort => switch to GO-WRITE (empirical + DE-risk paper); do not manufacture decorative theorems.
2. Any theorem that neither explains a demonstrated phenomenon nor validates reported inference gets cut (anti-decoration rule).
3. If referee simulation shows the story depends on promised-but-unproven theory, return to G5 and shrink claims to proven scope.

**Estimated effort:** 6-10 weeks. Compute: negligible beyond reruns.

---

## 8. Simulation study specification

### 8.1 Claims-to-experiments matrix

| Claim | Mechanism | DGP | Metric | Baselines | Ablation | Threshold | Falsifier | Output |
|---|---|---|---|---|---|---|---|---|
| C1 threshold | BBP-edge truncation vs. variance tradeoff | Gaussian factor panel, spike sweep (C1) | normalized counterfactual RMSE | oracle, donor-mean, SCM, CV-MC, SDID | ridge vs hard-threshold | preregistered kink within 15% of edge | smooth risk curve | fig_threshold_kink |
| C3 rank selector | TW gap separates spikes from bulk | same, rank unknown | rank accuracy | CV rank; Onatski ratio | with/without TW correction | +10pp over CV mid-regime | tie everywhere | fig_rank_accuracy |
| C3 pre-trends test | residual projection onto bulk under H0 | AR(1) null panels | empirical size/power | classical t-test | with/without bootstrap | size in [3,8]%; power >= 80% | size > 8% unfixable | fig_size_calibration |
| C2 DE risks | resolvent calculus | ridge family grid | DE-vs-empirical overlay error | n/a | n/a | overlay within finite-size tol | systematic deviation | fig_de_overlay |
| C5 lens | frontier distance ordering matches fragility | real panels + placebos | d with CI | published analyses | alternative rank inputs | placebo separation | no separation | fig_distance_to_frontier |

### 8.2 DGP ladder
Deterministic near-noiseless -> correctly specified Gaussian -> r=0 null -> dense-weak-factor (baseline-favorable) -> strong-alignment (mechanism-favorable) -> crossover grid -> AR(1)/heteroskedastic misspecification -> structural-break violation -> runtime scaling. Each rung answers one preregistered question; no decorative scenarios.

### 8.3 Fair comparison protocol
Shared seeds (10000+i), identical preprocessing, equal CV budgets, same hardware class per comparison, oracle variants labeled diagnostic-only, failed runs logged with reason, uncertainty across 500 replications reported (median + IQR), worst-regime highlighted, runtime and peak memory per method recorded.

---

## 9. Applied study specification

See Phase D. Identification considerations specific to panels: treatment timing and anticipation effects (exclude anticipation windows; sensitivity to window choice), composition of donor pool (placebo-in-space), serial-correlation structure (block bootstrap by time blocks), interference/spillovers (state-level spillover sensitivity note for Prop 99), and the structural-break caveat surfaced as a first-class diagnostic rather than buried in assumptions.

Application gate memo requirements: trusted-benchmark reproduction, acceptable identification diagnostics, survival of predeclared sensitivity checks, demonstrated change in understanding versus incumbent analysis, honest uncertainty, anomalies investigated before being framed as discoveries.

---

## 10. Deferred theory program

Contained in Phase E's target table (Section 7). Principle: no substantial proof work before Gates G3 and G4 pass; enabling items in Phase B are limited to the tagged conjecture formula and its special-case sanity algebra.

---

## 11. Compute policy and Google Colab sharding

**Local machine:** 13th Gen Intel i9-13900H (10 physical cores, 20 threads). Rules: process-level parallelism across replications only; cap workers at 8 to leave headroom; set BLAS/OMP threads to 1 inside workers (no nested parallelism); check current load before launching; RAM ceiling respected (matrices up to 5000 x 5000 float64 = 200 MB are safe; avoid materializing many such simultaneously).

**Colab trigger.** Any experiment whose one-seed-pilot-based projection exceeds EITHER 2 hours wall time on the local machine OR 4 GB peak RAM is routed to Google Colab notebooks instead of being run locally.

**Sharding budget.** Up to 40 independent, fully self-contained notebooks. Properties of every notebook:
1. Generates all data internally from explicit seeds (no external state, no cross-notebook dependencies);
2. Runs to completion within Colab limits (well under the ~10 h session cap; target <= 4 h per shard);
3. Writes its outputs (parquet/csv + a JSON metadata header with grid coordinates, git commit, library versions);
4. Ends with the mandatory download fallback:
```python
try:
    from google.colab import files
    files.download(output_file)
    print("Downloaded:", output_file)
except Exception as e:
    print("(Not on Colab / download skipped):", e)
```
5. Named `nb_<experiment>_shard<NN>_of<NN>.ipynb` with a deterministic shard-to-grid-cell map stored in `shard_manifest.yaml`.

**Merge discipline.** After download, a local `merge_shards.py` verifies completeness against the manifest (every expected grid cell x seed block present, checksums match) before any figure is drawn. Missing shards are re-run individually; figures are never produced from partial grids.

**Expected split under current design:** WP-C1 decisive grid and WP-C2 null/baseline cells stay LOCAL (projected < 2 h each with 8 workers); WP-C2(iii) AR(1)+bootstrap battery and WP-C4 large-n scaling cells are the primary COLAB candidates (bootstrap resampling multiplies runtime; 2000+ sized SVD batteries multiply memory). Final classification happens in WP-B3's cost report.

---

## 12. Risk register

| Risk | Probability | Damage | Earliest detector | Prevention | Recovery | Terminal? | Owner |
|---|---|---|---|---|---|---|---|
| Direct prior-art hit in remaining vocabularies | Low-med | High | WP-A2/A3 logs | Broad query families incl. econometrics jargon | Pivot to the untouched sub-question | Terminal if direct hit on C1+C3 | A3 |
| Weak-baseline illusion (win only vs untuned baselines) | Med | High | WP-C2(ii) baseline-favorable cell | Equal tuning budgets; strongest incumbents from day one | Add tuning; re-run | No | C1/C3 |
| Mechanism absent (smooth risk, no kink) | Med | High | WP-C1 kink criterion | Preregistered falsifier | Max two ansatz revisions | Terminal on third failure | C1 |
| TW miscalibration under serial correlation | High | Med-High | WP-C2(iii) size plots | Bootstrap fallback designed upfront | Ship bootstrap-native tests; limit theory to iid | No (scope change) | C2/C3 |
| Simplex SCM instability | Med | Low | WP-B2 tests | Ridge-first scope | Restrict scope, document | No | B2 |
| Tuning leakage (rank chosen using post-period info) | Low | High | WP-B2 leakage guard test | Hard code-path separation | Fix and rerun affected | No | B2 |
| Real-panel identification failure (diagnostic fires on true effects) | Low-Med | High | WP-D1/D2 controls | Positive-control battery | Investigate before reframing; repair or retreat | Terminal if unrepairable | D2 |
| Data access friction | Low | Med | Phase B read-only audit | Public datasets chosen (Replogle-style scale not needed here) | Swap canonical panel | No | D1 |
| Decorative-theory drift | Med | Med | G5 review against target table | Anti-decoration rule in Phase E | Cut T4-type targets | No | E1 |
| Compute overrun / crashed local runs | Med | Low | WP-B3 cost model | Colab routing per Section 11 | Re-shard | No | B3/C* |

---

## 13. Reproducibility and artifact map

```text
research/idea5/
  model_card.md                  # WP-A1: formal model, assumptions, estimands
  evidence_register.md           # WP-A2/A3: source ledger with anchors + E-levels
  priorart_deepread_memo.md      # WP-A2
  novelty_search_log.md          # WP-A3
  frontier_ansatz.md             # WP-B1: conjecture + special cases
  preregistration.md             # frozen before decisive Phase C runs
  pilot_cost_report.md           # WP-B3: LOCAL/COLAB classification
  shard_manifest.yaml            # Colab shard map (when triggered)
  code/scm_frontier/             # WP-B2 package (dgps, estimators, diagnostics)
  code/tests/test_estimators.py
  notebooks/                     # witnesses + Colab shards nb_*_shardNN_ofNN.ipynb
  results_c1/, results_c2/       # parquet shards + merged tables (schema: results_schema.yaml)
  figures/                       # every figure regenerated by scripts/make_figures.py
  gate_g3_memo.md                # WP-C5
  preprocessing_frozen.md        # WP-D1
  gate_g4_memo.md                # Phase D decision
  theory_targets.md              # WP-E1 proof packages
```

Seed policy: global seed registry `seeds.yaml`; every table/figure states its seed range and commit hash. Environment lock: `environment.yml` + pip freeze snapshot. Raw vs processed boundary: raw panel extracts immutable under `data/raw/`; all cleaning scripted into `data/processed/`. Checkpointing: every Colab shard writes incrementally (partial CSV append every 25 replications) so a crashed session loses at most one chunk.

---

## 14. Immediate actions (stop at next unresolved gate)

Only Phase A work is authorized now:

1. WP-A1: write `model_card.md`, then build and run the two witness notebooks (`witness_subthreshold.ipynb`, `witness_misalignment.ipynb`).
2. WP-A2 (parallel): deep-read 1710.10251, 2106.02780, 2109.15154 guarantee sections; update `evidence_register.md` with anchors.
3. WP-A3 (parallel): run the six residual query families; verify the dossier's 2026-ID citations; log everything in `novelty_search_log.md`.
4. Convene the G0+G1 decision using the Phase A give-up rules. Nothing in Phases B-E starts before that decision.

---

## 15. References (verification level at planning date 2026-08-23)

Verified today via arXiv API (abstract inspected):
- Athey, Bayati, Doudchenko, Imbens, Khosravi (2021), "Matrix Completion Methods for Causal Panel Data Models," JASA 116(536). https://arxiv.org/abs/1710.10251 , https://doi.org/10.1080/01621459.2021.1891924
- Farias, Li, Peng (2021+), "Learning Treatment Effects in Panels with General Intervention Patterns." https://arxiv.org/abs/2106.02780
- Agarwal, Han, Saha, Syrgkanis, Yoon (2022+), "Synthetic Blips." https://arxiv.org/abs/2210.11003
- Agarwal, Dahleh, Shah, Shen (2021), "Causal Matrix Completion." https://arxiv.org/abs/2109.15154

Classic/theory sources cited from the dossier, to be anchored (exact theorem + DOI) during WP-A2/A3 before any go decision relies on them:
- Baik, Ben Arous, Peche (2005), Ann. Probab. [BBP transition]
- Benaych-Georges, Nadakuditi (2011 Adv. Math.; 2012 JMVA) [outlier locations and overlaps]
- Johnstone (2001, Ann. Statist.) [TW largest eigenvalue]; Johnstone (2008, Ann. Statist.) [Jacobi ensembles]
- Onatski (2010, Rev. Econ. Stat.) [eigenvalue-ratio factor testing]; Onatski, Moreira, Hallin (2013, Ann. Statist.)
- Dobriban, Wager (2018, Ann. Statist.) [ridge risk]; Hastie et al. (2022) [ridgeless limits]
- Abadie, Diamond, Hainmueller (2010, JASA); Abadie (2021, JEL); Arkhangelsky et al. (2021, AER) [SDID]; Li (2020, JASA) [SCM inference]
- Bai (2009, Econometrica) [interactive fixed effects]; Ait-Sahalia, Xiu (2019, J. Econometrics)
- Moon-Weidner dynamic-panel line [UNVERIFIED, inspect in A3]
- Dossier 2026 IDs 2603.24833, 2605.30319 [UNVERIFIED, verify or drop in A3]

Zero-hit collision searches (arXiv API, 2026-08-23): listed in Section 4. Google Scholar/RePEc surfaces pending (WP-A3); no go decision may rest on the arXiv zero-hits alone.
