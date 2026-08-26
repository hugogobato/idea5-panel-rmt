# Preregistration D2b Addendum — Wave-2 panels (hypothesis-driven extension)

**Date frozen:** 2026-08-26, BEFORE any Wave-2 computation. **Parents:** `preregistration_d2_addendum.md` (machinery, unchanged), `gate_g4_memo.md` Sections 4-5 (rule-3 PIVOT indication; user-approved exploration of alternative panels). **Declared selection bias:** Wave-2 panels are chosen BECAUSE the frontier is expected to bind there; these are hypothesis-driven confirmations, not neutral replications. Any claim produced here carries that label. All Wave-1 thresholds, instruments, labels, and leakage rules apply verbatim unless explicitly amended below; failures are reported as failures.

---

## W1. Basque panel (natural candidate)

Data: Abadie-Gardeazabal (2003) Basque-country panel, outcome real GDP per capita, extracted from the Rdatasets mirror of the R `Synth` package dataset (`basque.csv`); provenance level E2 (community mirror), sha256-pinned under `data/raw/`, structural validation at load time (region set, 1955-1997 grid, no missing outcome entries). If a primary replication archive surfaces later, it supersedes the mirror and triggers re-validation.

Frozen configuration: treated unit Basque Country; treatment year 1975 (AG2003 terrorism-onset convention); pre window 1961-1974 (T0 = 14), post 1975-1997; donors = all other Spanish regions in the extract. Boundary sensitivity covered by the inherited trimmed-window arm.

Wave-2 specification amendment (deviation D-W1, made HERE, before runs): placebo and treated pre-fit RMSE are computed with OUTCOME-ONLY classic SCM (simplex-constrained least squares on the pre window via the frozen augmented-NNLS solver, `synth_adh._nnls_simplex_ridge`, center=True; SLSQP remains banned per preprocessing_frozen.md Section 4.1), because replicating the AG2003 predictor architecture is out of scope for a lens study. No published-effect anchors are gated in Wave-1 style; data integrity is enforced structurally instead. V-search machinery is therefore not used anywhere in Wave-2.

Gates (all inherited): treated pipeline (k, d, CIs, p_align, classification); placebo battery with P1 threshold rho <= -0.42-class value interpolated at the observed donor count via the same one-sided permutation test (threshold fixed at -0.45 for n=17, the conservative neighbor of the frozen table); P2 (d_treated >= Q90); P3 z_shift date-alarm rate <= 0.20 at pseudo-cutoffs {1969, 1971} (indices 8, 10 of the pre window, T_post = 6); sensitivity stability >= 3/4. NF arms N1/N2/N3 evaluated as in Wave 1 (poor-looking predicate: treated pre-fit worse than Q75 of placebo pre-fits).

W1-PASS (novel-finding claim revived for this panel class) iff the treated unit classifies FRAGILE-MISALIGNED or INCONCLUSIVE-SUBEDGE/FRAGILE-SUBEDGE with donor-bootstrap support, AND P3 passes AND stability >= 3/4. FRAGILE-INVISIBLE additionally requires posttest_control reported (scope note: silence + drift would mean the pre basis cannot transport, feeding the certification story instead). Every other combination is reported honestly and does not revive C5.

## W2. Donor-pool degradation battery on smoking (positive-sensitivity control)

Purpose: bridge between simulation and natural panels — does the frozen lens respond when signal mass is genuinely scarce on OTHERWISE-REAL data? This manipulates the pool and is labeled a pool-manipulation control, never a natural finding.

Design: from the unit-centered smoking pre matrix, draw random donor subsets WITHOUT replacement of sizes n_d in {38 (full), 24, 16, 12, 8, 6, 4, 3}, 20 seeded draws per size (seeds 61301 + 100*size_rank + j); treated row untouched. Primary-selector pipeline per draw with bootstrap DISABLED (adaptation D-W2: lite labels k, d, p_align only — CIs are irrelevant to a silence transition and this keeps the battery minutes-scale under current machine load).

Reported per size: median d-hat, share of draws with k = 0, share misaligned. Control criterion: the battery demonstrates sensitivity iff some size >= 4 yields >= 50% silent draws while the full panel is loud (k >= 1); the transition location is reported as a descriptive curve. Failure mode (instrument insensitive until n_d = 3): recorded as such.

## Seeds

`seeds.yaml` phase_d.wp_d2_wave2: baseline instruments inherit Wave-1 offsets (align 60101/donor 60102/time 60103 + idx x 10^6, perm 60104); Basque treated uses idx = 0, its placebos idx = 1..n; degradation battery 61301+; z_shift 8880001; post-test null pools 7770001.

## Compute

Local, minutes per panel; no Colab trigger. Machine currently heavily loaded by unrelated jobs (measured load ~49/20 threads): runs are sequential, single-threaded BLAS, and sized to tolerate it.

## Deviation log

D-W1 and D-W2 above (both pre-run).

D-W3 (implementation defect, logged 2026-08-26 during the FIRST W1 execution): the placebo loop mistakenly kept the TRUE treated row (Basque) as the regression target for every placebo, producing leave-one-donor-out refits of Basque rather than in-space placebos; detected immediately via degenerate duplicate pre-fit values (most dropped donors carry zero weight). No scientific conclusion was drawn from the defective run. Fix restores the frozen Wave-1 placebo definition verbatim (placebo target = unit i's own pre row; pool = all units except {i, true treated}); the run was restarted from scratch afterward.
