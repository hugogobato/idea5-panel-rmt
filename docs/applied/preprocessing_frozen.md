# Preprocessing Freeze (WP-D1): Phase D application panels

**Project:** Idea 5, spectral recoverability frontiers.
**Work package:** WP-D1 (`docs/plan/Idea5_Panel_RMT_Research_Plan.md`, Section 7, Phase D).
**Status:** FROZEN 2026-08-26, after trusted-benchmark reproduction succeeded and BEFORE any comparative or distance-to-frontier analysis (WP-D2). No rule in Sections 2 to 5 may change after this date; any need to touch ingestion, cleaning, estimator specs, anchors or tolerances requires a new dated addendum first. Failures discovered later are reported as failures, never tuned toward.

---

## 1. Purpose and scope

This document fixes how the two Phase D panels enter the project and which published numbers count as the trusted benchmarks. WP-D1 demonstrated that the frozen pipeline below reproduces those benchmarks to reporting precision (Section 6); WP-D2 consumes its outputs under the leakage rules of Section 7 without re-touching anything defined here.

## 2. Data provenance (immutable raw layer)

Raw extracts live under `data/raw/`, checksums in `data/raw/SHA256SUMS`. Download date: 2026-08-26.

| File | Source | sha256 (prefix) |
|---|---|---|
| `repgermany_updated.tab` | Abadie, Diamond, Hainmueller (2015) replication archive, updated erratum version, Harvard Dataverse doi:10.7910/DVN/24714, file id 13454770 | `0ec596dc` |
| `codebook_repgermany_updated.md` | same archive, file id 13578067 | `96fd3e4f` |
| `README_ADH2015_dataverse.md` | same archive, file id 13578064 (erratum note) | `a57f5928` |
| `rep_updated.r` | same archive, file id 13595688 (author replication code) | `9e58393e` |
| `california_prop99_synthdid.csv` | synthdid reference implementation repo, `data/california_prop99.csv`, github.com/synth-inference/synthdid@master | `550d24c2` |
| `california_prop99_processing.r` | same repo, processing script proving the extract's construction | `06b1307c` |
| `smoking_adh_covariates_facure.csv` | community mirror of the ADH (2010) smoking panel with covariates (python-causality-handbook), cross-validated against the synthdid extract at load time | `bae31804` |

Primary-source anchors were pinned from the journals themselves: ADH (2010), JASA 105(2), 493-505 (full text inspected; abstract, Tables 1-2, Figures 3-8 statements); ADH (2015), AJPS 59(2), 495-510 (full text inspected; Tables 1-2, Figure 3 statements); Arkhangelsky, Athey, Hirshberg, Imbens, Wager (2021), arXiv:1812.09970v4, Table 1 and Eqs. (2.1)-(2.4), Algorithm 1; the synthdid R reference implementation (`R/synthdid.R`, `R/solver.R`, `R/utils.R` at master, 2026-08-26) was used to mirror the estimator semantics exactly.

## 3. Frozen cleaning contract

### 3.1 California smoking (Prop 99)

Grid: 39 states x 31 years, 1970-2000. Treated unit: California, from 1989 (T0 = 19, T_post = 12). Donors: the 38 remaining states of the ADH donor list (Alabama ... Wyoming as coded in `loaders.SMOKING_DONOR_ORDER`). Outcome: packs per capita (`cigsale`). No rows are dropped, imputed, smoothed, or transformed; the two independent mirrors must agree on the full state-year outcome grid to <= 1e-6 (enforced at load time; measured max diff 0).

ADH (2010) predictor block X (paper Table 1 order): ln(per-capita personal income, 1997 dollars) mean 1980-1988; percent aged 15-24 mean 1980-1988 (mirror stores a fraction, scaled by 100); average retail price per pack mean 1980-1988; beer consumption per capita mean 1984-1988; cigarette sales levels for 1988, 1980, 1975. Treated-row means reproduce paper Table 1 to <= 0.01 on every row except percent aged 15-24, where the mirror gives 17.353 vs the printed 17.40 (documented source nuance, tolerated).

### 3.2 German reunification

Grid: 17 countries x 44 years, 1960-2003, complete. Treated unit: West Germany (index 7), from 1990 (T0 = 30, T_post = 14). Donors: the other 16 OECD countries of the author archive. Outcome: GDP per capita, PPP current international USD (erratum labeling; not deflated). Predictors follow `rep_updated.r` verbatim, two stages:

1. Training stage X: gdp, trade, infrate each averaged 1971-1980; industry share averaged 1971-1980; schooling averaged over {1970, 1975}; investment rate (invest70) at 1980. V is fit against the outcome window 1981-1990 inclusive (the authors' `time.optimize.ssr`, which intentionally spans the treatment year).
2. Main stage X: gdp, trade, infrate averaged 1981-1990; industry averaged 1981-1990; schooling averaged over {1980, 1985}; investment rate (invest80) at 1980. W is re-solved at the training-stage V, frozen.

## 4. Estimator specifications

### 4.1 Classic SCM (ADH mode), instrument A1/B1

Synth conventions replicated exactly: predictor rows of [X0 | X1] are divided by their combined sample SD (ddof = 1); unit weights solve min (X1 - X0'w)' V (X1 - X0'w) subject to w >= 0, sum(w) = 1; predictor weights solve min_v MSPE(Z1 - Z0'w(v)) over the preregistered outcome window Z (smoking: full pre path 1970-1988, single stage, per paper Section 2.3; Germany: the two-stage procedure above). V search: deterministic multistart (uniform, all one-hots, 16 seeded Dirichlet draws, seed 60001) Nelder-Mead on a softmax parameterization, plus restart polish at the incumbent optimum. Smoking uses the default search depth (16 starts, 250 iterations, 1 polish round); see deviation D3 on why a deeper search is NOT shipped.

Inner simplex QP solvers (deviation D1): SLSQP terminates at its initial point reporting success on this problem class and is therefore banned; exact Lawson-Hanson NNLS on an augmented system (equality encoded as a penalty row at 1e6 x design scale, renormalized afterwards) is used everywhere a plain simplex QP is needed. The synthdid estimators use the reference implementation's own sparsified Frank-Wolfe port instead (below).

### 4.2 synthdid instruments, A2

Ported 1:1 from `synthdid` master (`sc.weight.fw`, `fw.step`, `sparsify_function`, estimate formula): Frank-Wolfe with exact line search from uniform start, phase 1 of 100 iterations, then the reference sparsify tie-break (zero weights <= max/4, renormalize), then up to 10,000 iterations; stopping when the objective decrease falls below (1e-5 sigma_hat)^2 with sigma_hat = sd of pooled donor pre-period first differences (ddof = 1). SC: omega intercept OFF, zeta_omega = 1e-6 sigma_hat, lambda identically zero, tau = post-window mean of the level gap. SDID: both intercepts ON, zeta_lambda = 1e-6 sigma_hat, zeta_omega = (N_tr T_post)^{1/4} sigma_hat, tau via the reference closed form c(-omega, 1)' Y c(-lambda, 1/T_post). DID: uniform weights (sanity instrument only).

## 5. Anchors, tolerances, gating

Published values live in `scripts/wp_d1_reproduce.py` (`ANCHORS`) next to the tolerances. Gating rules:

1. Donor weights: tolerance +/- 0.01 (smoking) and +/- 0.02 (Germany). Papers print two decimals; our solver swap (ipop vs NNLS/FW) leaves third-decimal residue, observed <= 0.0045 on every gated weight.
2. Effect anchors: windows proportional to how the papers phrase them: smoking average gap "almost 20" (+/- 1.5), year-2000 gap "about 26" (+/- 3), pre-MSPE "about 3" (+/- 50%), post/pre MSPE ratio "~130" (+/- 60%); Germany average gap "-1600 USD/yr" (+/- 250), relative 2003 gap "+12%" (+/- 4 pp, sign fixed positive: synthetic ABOVE actual).
3. synthdid Table 1: SC -19.6 (+/- 0.4), SDID -15.6 (+/- 0.5), DID -27.3 (+/- 0.5). These tolerances absorb the papers' own tie-break ambiguity on rank-deficient weight fits (19 pre periods, 38 donors); the reference FW tie-break is mirrored precisely and lands within 0.02 and 0.004 respectively.
4. REPORTED ONLY, never gated: predictor-weight vectors V (outputs of the authors' ipop runs, not independent claims); the ADH Table 1 synthetic balance column and the ADH 2015 Table 2 synthetic column. Justification (measured, not asserted): the smoking pre-MSPE objective admits near-equivalent optima whose donor weights differ beyond published rounding while effects barely move; reallocating 0.002-0.005 mass among the five named donors moves the cigsale-1975 balance row across 126.5-127.7, which brackets both our value (127.10) and the printed 126.99. Gating derived columns of this kind would gate solver noise, not substance.

## 6. Achieved reproduction (WP-D1 verdict inputs)

Driver: `python3 scripts/wp_d1_reproduce.py` (writes `results_d1/reproduction_{smoking,germany}.json`, series CSVs, figures `fig_d1_repro_{smoking,germany}.png`). Exit code 0 iff every gate passes. Result at freeze time: 21/21 gates PASS, 0 FAIL.

Smoking: weights CO 0.1595 / CT 0.0679 / MT 0.2019 / NV 0.2356 / UT 0.3351 (targets .164/.069/.199/.234/.334, all within 0.0045, zeros elsewhere exact); pre-MSPE 3.077; average gap 1989-2000 -18.98; gap 2000 -25.73; ratio 129.0; synthdid SC -19.62, SDID -15.60, DID -27.35.

Germany: weights AT 0.4184 / US 0.2205 / JP 0.1577 / CH 0.1090 / NL 0.0944 (targets .42/.22/.16/.11/.09, all within 0.0046, zeros elsewhere exact); average gap 1990-2003 -1587.6 USD/yr; relative gap 2003 +11.75%. Pass-through diagnostic at the paper's own V reproduces the same five donors and gives average gap -1585.7 and relative gap +11.75%, confirming pipeline correctness independently of the V search.

Per plan WP-D1: PASS. The canonical pipelines reproduce the published point estimates to reporting precision.

## 7. Deviations and findings recorded at freeze time

1. D1 (solver swap, before any anchor was matched): SLSQP is unreliable as a simplex-QP solver here (terminates at the initial point reporting success); replaced by augmented-NNLS for classic-SCM inner problems and by the faithful Frank-Wolfe port for the synthdid instruments. Both choices precede and are independent of the anchor outcomes they were tested against.
2. D2 (report-only derived columns): as justified in Section 5.4, with the measured sensitivity experiment recorded in the repository log of 2026-08-26.
3. D3 (search depth basin): a deliberately deeper V search (96 starts, 1200 iterations) finds a strictly better pre-MSPE optimum (2.867 vs 3.077) in a DIFFERENT basin (V collapses onto the 1975 lag; Colorado drops out of the weight support), which no longer matches the published weight table. The shipped instrument therefore uses Synth-default search depth targeting the published basin; the deeper-basin numbers are retained in the repository log as evidence that the ADH point estimates are one representative of a flat objective rather than its unique argmax. This is itself a small frontier-relevant observation: classic SCM's published solutions can be V-flat, which the WP-D2 lens will quantify.
4. D4 (Germany SSR window): the authors' updated script fits the training stage on outcomes 1981-1990 INCLUSIVE, i.e., spanning the treatment year; replicated verbatim rather than sanitized.
5. D5 (anchor sign fix, pre-freeze): the 2003 relative-gap anchor was initially drafted with the wrong sign during implementation; corrected against the paper sentence ("in 2003, per capita GDP in the synthetic West Germany is estimated to be about 12% higher than in the actual West Germany") before any freeze-time conclusion relied on it. Recorded for transparency.

## 8. Rules inherited by WP-D2

Leakage: every spectral quantity (scree, spike estimates m-hat, alignment energy, gate/test statistics) is computed from PRE-treatment data only; treated post-period outcomes enter solely through the estimator comparison objects already defined in Phases B/C. Distance-to-frontier d, placebo batteries, and incumbent comparisons must consume the panels exclusively through `applications.loaders` so the contract above cannot drift. Any new decisive D2 criterion gets its own dated preregistration addendum before running, per gate_g3_memo Section 7 discipline.

## 9. Reproducibility

Seeds: `config/seeds.yaml` `phase_d.wp_d1.v_multistart_seed = 60001` (numpy PCG64 via `default_rng`; the only stochastic input in WP-D1). Environment: Python 3.12.3, numpy 2.4.3, pandas 3.0.1, scipy 1.17.1, matplotlib 3.10.8 (matches `config/seeds.yaml` environment block). Tests: `pytest code/tests/test_wp_d1_reproduction.py` (6 tests) re-runs both reproductions end-to-end and asserts the Section 5 gates; full suite `pytest code/tests` (32 tests) green at freeze time.
