# Evidence Register (WP-A2/A3)

**Project:** Idea 5, panel causal inference phase transitions.
**Verification levels:** E0 = unverified mention; E1 = citation located but content uninspected; E2 = abstract/metadata verified (date noted); E3 = full-text guarantee sections inspected with anchors.
**Dates of verification:** 2026-08-24 (arXiv API + Crossref + full-text reads). Machine-assisted full-text extraction; anchors quoted from arXiv versions as noted.

---

## A. Closest prior art (deep-read, E3)

### S01. Athey, Bayati, Doudchenko, Imbens, Khosravi (2021), "Matrix Completion Methods for Causal Panel Data Models," JASA 116(536)
- Links: https://arxiv.org/abs/1710.10251 , DOI 10.1080/01621459.2021.1891924 (verified)
- Version read: arXiv v5 (2022-04-21), JASA-corresponding.
- Estimand: ATT over missing blocks of Y(0); Y = L* + eps (Eq. 4.1); estimator MC-NNM (nuclear norm + two-way fixed effects, Eq. 4.3).
- Guarantees: Theorem 1 (Sec. 5) is an equivalence/representation result among MC-NNM, horizontal/vertical regression, SC, DID; no rates. Theorem 2 (Sec. 6.2): RMSE order-rate bound depending on L_max = ||L*||_max, rank R, noise sigma, control-propensity p_c; universal constants only.
- Assumptions: Assumption 1 (Sec. 4): iid sigma-sub-Gaussian noise independent of L*. Assumption 2 (Sec. 6.1): adoption dates independent of eps given L*; p_c = min_i pi_T^(i). No singular-value lower bound, no incoherence condition in the theorem.
- Lower bounds: none (only informal remark Sec. 6.1 about row-supported rank-1 matrices being unrecoverable when t_i < T).
- Inference: none asymptotic ("beyond the scope", Sec. 4); CV-split ad hoc intervals; no rank-selection guarantee.
- Admitted limitations: Secs. 8.3-8.6 (serial correlation, propensity weighting, p_c conservatism, approximate factors unproven).
- Answers: (a) per-unit BBP-type recoverability threshold: **NO**. (b) TW/eigenvalue-calibrated diagnostics: **NO**. (c) exact high-dimensional limits (MP/TW/DE): **NO**.
- Role: strongest incumbent; owns the unification story. Our delta must be the quantitative frontier + calibration layer.

### S02. Farias, Li, Peng (2021+, v2 2023), "Learning Treatment Effects in Panels with General Intervention Patterns"
- Link: https://arxiv.org/abs/2106.02780
- Version read: arXiv v2 (2023-03).
- Estimand: ATT only (tau_m, Secs. 1.1, 2); unit-specific trajectories are nuisance.
- Guarantees: nuclear-norm regularized LS + de-biasing; Theorem 1 rate depends on kappa = sigma_max/sigma_min, r, ||Z||_F; Proposition 2 minimax lower bound matching up to logs. Identification impossibility (Proposition 1) is TANGENT-SPACE geometric (Assumption 1 conditions on ZV*, Z'U*), not spectral/spikiness-based.
- Inference: Theorem 2 asymptotic normality of de-biased ATT under sigma_min = Omega(n), incoherence max-row norms O(sqrt(r/n)).
- Answers: (a) **NO** (threshold is geometric identifiability, not a computable spikiness boundary tied to a spectrum). (b) **NO**. (c) **NO** (operator-norm concentration only).
- Role: shows rate-optimal ATT work exists; leaves per-unit detectability frontier open. Their lower-bound technology is a template for our T2 target.

### S03. Agarwal, Dahleh, Shah, Shen (2021), "Causal Matrix Completion"
- Link: https://arxiv.org/abs/2109.15154
- Version read: arXiv v1 (2021-09-30, only version).
- Estimand: entry-wise potential outcomes A_ij = E[Y_ij | u_i, v_j] (Eq. 10) under arbitrary missingness D (no model of P at all); panel treated-block is contrasted in Secs. 3.1/3.3.
- Guarantees: SNN estimator (anchor bicliques + PCR); Theorem 2 / Corollary 1 entry-wise O_p rates (order-only, N^{-1/4}-type leading term); Assumptions A1-A7 include well-balanced spectra (A6) and linear span inclusion (A3).
- Lower bounds: **NONE PRESENT** (zero occurrences of lower bound/minimax/impossibility in full text; theorem inventory checked). Note: some secondary summaries claim minimax lower bounds here; that claim is false for v1.
- Inference: Theorem 3 CLT for the estimator itself per entry, conditional on instance (K -> infinity, K = o(N^{1/2}) cap).
- Answers: (a) **NO**. (b) **NO** (Gavish-Donoho cited only as related work). (c) **NO**.
- Role: adjacent (MNAR geometry); their A6 balance assumption is exactly where our spiked analysis will live instead.

## B. Nearest neighbors verified at abstract level (E2, 2026-08-24)

### S07. Spiess, Imbens, Venugopal (2023+, v3 2023-10), "Double and Single Descent in Causal Inference with an Application to High-Dimensional Synthetic Control"
- Links: https://arxiv.org/abs/2305.00700 , DOI 10.3386/w31802 , DOI 10.48550/arXiv.2305.00700
- Version read: arXiv v3 full text (HTML rendering), E3. Deep-read executed in Phase B (medium priority item closed).
- Content: (i) CPS/LaLonde wage imputation: double descent in number of randomly ordered covariates l vs training size n=3000, minimal-norm interpolating LS; loss peaks at interpolation threshold l=n then descends. (ii) California smoking (Abadie et al. 2010 data), T=3 pre-years, varying subsets of N=20 controls: SINGLE descent, RMSE monotonically decreasing in number of control units, no regime change past l=T. (iii) Theory is MECHANICAL: Prop 1 model averaging for interpolating LS; Prop 2 variance reduction beyond interpolation; Prop 3 variation hierarchy; Prop 4 model averaging for SC (driven by convexity of weights); Sec. 4 Jensen-type risk bounds from the (MA) property under high-level assumptions, "largely agnostic about the true data-generating process".
- Answers: (a) per-unit spikiness threshold: **NO** (complexity axis = donor count vs pre-period count; no eigenvalue quantities anywhere; keyword sweep for eigenvalue/singular/spectral/Marchenko/Tracy/low-rank/factor: zero substantive hits). (b) TW diagnostics: **NO**. (c) exact RMT limits: **NO**.
- POSITIONING PARAGRAPH (mandatory in any draft, watch W-2; canonical text recorded here):
  Spiess et al. vary ESTIMATOR COMPLEXITY (number of donors relative to pre-treatment periods) at fixed signal; we vary SIGNAL STRENGTH and treated-unit ALIGNMENT at fixed estimator class. Their descent curves describe the weight-estimation variance of interpolating/convex-weight imputators via model averaging; our frontier describes the INFORMATION CONTENT of the donor panel for the treated row (BBP detectability). The axes are orthogonal and complementary: (1) their single-descent guarantee is consistent with our supercritical regime, where more donors monotonically improve factor recovery toward the sigma floor, a limit invisible to their distribution-free machinery; (2) our channel-2 result (Witness 2) is a prediction their framework cannot make: if the treated loading is orthogonal to the spiked space, NO number of donors helps, and their monotone-improvement curve must flatten at the noise floor; (3) in sub-threshold panels the two theories give disjoint explanations of poor pre-fit (absent signal vs weight variance), yielding a joint empirical test: descent curves flattening AT the sigma floor indicate a binding frontier, not insufficient complexity.
- Watch item W-2 status: CLOSED (canonical positioning text frozen here; to be carried verbatim-adapted into any manuscript).

### S08. Mehrotra, Tran, Vu, Zampetakis (2026), "Improved Guarantees for Heterogeneous Treatment-Effect Estimation via Matrix Completion"
- Link: https://arxiv.org/abs/2605.30319 (dossier ID RESOLVED; exists), DOI 10.48550/arXiv.2605.30319
- Version read: arXiv v1 (2026-05-28), full text (HTML rendering), E3. Deep-read executed at Phase B entry per watch item W-1.
- Setting: panel experiments, D_ij ~ Ber(p_ij) with UNKNOWN non-uniform propensities; signal-plus-noise Y(a) = A(a) + E(a); estimand is the treatment-effect matrix M = A(1) - A(0), controlled in row-wise norm (1/sqrt(m)) ||M - M_hat||_{2,infty}.
- Guarantees: Theorem 3.2 (Sec. 3) row-wise l2 upper bound for a row-scaled truncated-SVD estimator; Corollary 3.3: max_i ||M_hat_i - M_i||_2 / sqrt(m) <~ K r^{3/2} mu log^4(m+n) [ sqrt(r_p/q (1/n + n/m^2)) + max_a ||P(a)||_op / sqrt(m min{m,n}) ]. Assumptions (Assumption 3.1): approximate low rank (sigma_{r+1} <~ K sqrt(m+n)), bounded entries, iid mean-zero noise, row/column incoherence mu, SNR floor sigma_1 >~ K r T(a).
- Answers: (a) per-unit BBP-type recoverability threshold: **NO**. The only necessity language is Sec. 1.1: "The incoherence and signal-to-noise assumptions are regularity assumptions that can also be shown to be necessary", i.e., regularity needed by THEIR upper bound (degenerate concentrated-signal case), not a spectral spikiness impossibility. Zero occurrences of minimax/impossibility/lower-bound theorems (full-text keyword sweep: "lower bound", "impossib*", "converse", "minimax" hit only the SNR assumption and sharpness discussion of their own perturbation bound vs operator-norm alternatives, Appendix B). (b) TW/eigenvalue-calibrated diagnostics: **NO**. (c) exact high-dimensional limits (MP/TW/DE): **NO**; machinery is operator-norm concentration + contour-expansion perturbation theory (Sec. 4, App. B).
- Relation to Idea 5 (positioning, recorded for any draft): complementary axes. S08 gives UNIFORM per-row upper bounds under entrywise random missingness + row-incoherent (spread-out) signal; Idea 5 gives a PER-UNIT spectral frontier whose whole point is the atypical treated row (alignment/leverage), under structured block-row missingness (SCM panels). Their Theorem B.2 (sharp ||.||_{2,infty} truncated-SVD perturbation bound, App. B) is registered as a candidate TOOL for targets T1/T2 adaptations.
- Watch item W-1 status: CLOSED (no G1 reopening trigger).

### S09. Shen, Song, Abadie (2025+), "Efficiently Learning Synthetic Control Models for High-dimensional Disaggregated Data"
- Link: https://arxiv.org/abs/2510.22828
- Multivariate square-root Lasso SC weights; error bounds under time-series dependence; ATT estimation. Norm-rate family; no spectral thresholds. Cite as incumbent-family member (Abadie co-author).

### S10. Ferman, Pinto (2021), "Synthetic Controls with Imperfect Pre-Treatment Fit"
- Links: https://arxiv.org/abs/1911.08521 , DOI 10.3982/qe1596 (Quantitative Economics)
- SC biased under imperfect pre-fit with correlated assignment; demeaned SC; specification test. This is the applied-econometrics owner of "pre-fit has statistical meaning"; our delta: computable spikiness frontier + TW-calibrated tests rather than bias analysis under correlated assignment. Cite as the primary diagnostic incumbent.

### S11. Roth (2022), "Pretest with Caution," AER: Insights, DOI 10.1257/aeri.20210236
- Pre-trend pretesting distortions; context for C3's type-I framing. Related: Roth et al. parallel-trends literature. E2.

### S12. Wang (2024+), "Counterfactual and Synthetic Control Method ... Instrumented Principal Component Analysis," https://arxiv.org/abs/2408.09271
- Method variant (generalized SC with instrumented loadings). No spectral limits. Low direct-hit risk. E2.

### S13. Yan, Chen, Fan (2021+), "Inference for Heteroskedastic PCA with Missing Data," https://arxiv.org/abs/2107.12365
- TOOL SOURCE: non-asymptotic inference for principal subspaces under spiked covariance with missing data + heteroskedastic noise (HeteroPCA base). Directly relevant to C3 robustness roadmap. E2.

### S14. Agarwal, Choi, Yuan (2026), "Robust Matrix Estimation with Side Information," https://arxiv.org/abs/2603.24833 (dossier ID RESOLVED; exists)
- Four-component matrix decomposition with side information; nuclear-norm machinery. Tool/adjacent; no causal panel estimands. E2.

### S04-S06. Moon-Weidner line and IFE panel PC (all E2 abstracts, arXiv)
- S04: Moon, Weidner, "Linear Regression for Panel With Unknown Number of Factors as Interactive Fixed Effects" (Econometrica 2015), https://arxiv.org/abs/2605.00614 . Limit distribution of LS estimator invariant to over-specified factor count; inference without consistent rank estimation.
- S05: Moon, Weidner, "Nuclear Norm Regularized Estimation of Panel Regression Models," https://arxiv.org/abs/1810.10987 . Convex estimators, consistency.
- S06: Peng, Su, Westerlund, Yang, "Interactive Effects Panel Data Models with General Factors and Regressors," https://arxiv.org/abs/2111.11506 .
- Assessment: consistency/rates and limit distributions for regression coefficients under many factors; not RMT-exact, no per-unit recoverability question, no calibrated eigenvalue tests for counterfactual targets. Register conclusion recorded per WP-A3 action item 3.

### S15. Agarwal, Dahleh, Sarkar (2019) dossier reference RESOLUTION
- Resolves to "A Marketplace for Data: An Algorithmic Solution," https://arxiv.org/abs/1805.08125 . Data-pricing mechanism design; NOT a causal-panel or spectral result. Verdict: miscite in dossier for our purposes; DROP from the causal lineage (keep only if data-marketplace angle ever needed).

## C. Theory tool sources (anchors for Phase E targets)

| ID | Source | DOI / URL | Verified |
|---|---|---|---|
| T1 | Baik, Ben Arous, Peche (2005), Ann. Probab. 33(3), 1643-1697 (BBP transition, real case) | arXiv https://arxiv.org/abs/math/0408040 ; DOI UNRESOLVED VIA CROSSREF (Euclid indexing gap). Complex-case companion verified: 10.1214/009117905000000233. Manual Euclid DOI lookup flagged. | E2 |
| T2 | Benaych-Georges, Nadakuditi (2011), Adv. Math. 227(1), 494-521 (outlier locations AND overlaps for low-rank perturbations) | DOI 10.1016/j.aim.2011.02.007 (verified; note: 10.1016/j.aim.2011.02.011 is a DIFFERENT paper) | E2 |
| T3 | Johnstone (2001), Ann. Statist. 29(2), 295-327 | DOI 10.1214/aos/1009210544 | E2 |
| T4 | Johnstone (2008), Ann. Statist. 36(4), Jacobi ensembles/TW limits | DOI 10.1214/08-AOS605 | E2 |
| T5a | Onatski (2010), Rev. Econ. Stat. 92(4), eigenvalue-ratio factor testing | DOI 10.1162/rest_a_00043 | E2 |
| T5b | Onatski (2009), Econometrica 77(4), factor-count testing | DOI 10.3982/ECTA6964 | E2 |
| T6 | Dobriban, Wager (2018), Ann. Statist. 46(6A), ridge risk DEs | DOI 10.1214/17-AOS1549 | E2 |
| T7 | Hastie, Tibshirani, Friedman, Wainwright (2022), Ann. Statist., ridgeless limits | DOI 10.1214/21-AOS2133 | E2 |
| T8 | Abadie, Diamond, Hainmueller (2010), JASA 105(490) | DOI 10.1198/jasa.2009.ap08746 | E2 |
| T9 | Arkhangelsky, Athey, Hirshberg, Imbens, Wager (2021), AER 111(6), SDID | DOI 10.1257/aer.20190159 | E2 |
| T10 | Bai (2009), Econometrica 77(4), interactive fixed effects | DOI 10.3982/ECTA6135 | E2 |
| T11 | Ait-Sahalia, Xiu (2019), J. Econometrics 211(1), PCA of high-frequency factor models | DOI 10.1016/j.jeconom.2017.08.015 | E2 |
| T12 | Li (2020), JASA 115(529), inference for ATT by SCM | DOI 10.1080/01621459.2019.1686986 | E2 |
| T13 | (bonus) "Testing in high-dimensional spiked models," Ann. Statist. 2020 | DOI 10.1214/18-AOS1697 | E2 |

Pending anchors (E1, harmless until Phase E relies on them): Abadie (2021), JEL survey; Onatski-Moreira-Hallin (2013, Ann. Statist.).

## D. Watch items

1. W-1: S08 (Mehrotra et al.) full-text guarantee deep-read due at Phase B entry. [CLOSED 2026-08-24, Phase B kickoff: full text read, E3; no impossibility/lower-bound statement implying a per-row spikiness threshold; G1 not reopened. See S08 entry.]
2. W-2: S07 (Spiess et al.) positioning paragraph mandatory in any draft: frontier vs double-descent axes.
3. W-3: Google Scholar and RePEc surfaces could not be queried programmatically from this environment (see novelty_search_log.md for the exact attempted URLs and the manual checklist). No go decision rests solely on the automated surface; manual pass remains a bounded human task.
