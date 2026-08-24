# Novelty Search Log (WP-A3)

**Date executed:** 2026-08-24. **Tooling:** arXiv API (export.arxiv.org/api/query), Crossref REST API, direct arXiv id_list lookups. Google Scholar offers no programmatic access; RePEc/IDEAS returned a non-parseable page (details below).

## 1. Collision re-checks (plan Section 4 queries re-run)

| Query (arXiv API syntax) | Hits | Assessment |
|---|---|---|
| all:"synthetic control" AND ("Marchenko-Pastur" OR "Marcenko-Pastur") | 0 | Zero-hit maintained |
| all:"synthetic control" AND all:"random matrix" | 0 | Zero-hit maintained |
| all:"synthetic control" AND ("Tracy-Widom" OR "BBP threshold") | 0 | Zero-hit maintained |
| all:"matrix completion" AND all:"Tracy-Widom" | 0 | Zero-hit maintained |
| all:"panel data" AND all:"spiked" | 1 (Bayesian sparse heterogeneity panel; unrelated) | No threat |
| all:"counterfactual" AND all:"spiked" | 8 (all unrelated: spiking neurons, fairness, diabetes) | No threat |

## 2. Six plan-mandated query families

| Family | Query | Hits | Outcome |
|---|---|---|---|
| factor augmented + prediction | "factor augmented" + "prediction" | 17 | Forecasting/methods papers only; no spectral recoverability content |
| PC + panel consistency (Moon-Weidner) | au:Moon AND au:Weidner AND cat:econ.EM; plus "interactive effects"+"panel"+"principal components" | 5 + 1 | Moon-Weidner line pinned: 2605.00614 (Econometrica LS unknown factors), 1810.10987 (nuclear-norm panels), plus Peng-Su-Westerlund-Yang 2111.11506. Recorded in register S04-S06. Not RMT-exact; no per-unit question |
| spiked covariance + missing data | "spiked covariance" + missing | 2 | Yan-Chen-Fan HeteroPCA inference = tool source S13 |
| denoising + spiked | spiked + denoising | 77 (mostly spiking-neuron noise) | One relevant tool hit: spiked F-matrix shrinkage 2211.00986; tool-level only |
| SC + high-dimensional asymptotics | "synthetic control" + high-dimensional | 13 | Spiess et al. double descent (S07), Shen-Song-Abadie sqrt-lasso (S09), Wang instrumented PCA (S12); no thresholds/TW anywhere at abstract level |
| pre-treatment fit + test | "pre-treatment fit"; "pre-trends"+test+econ.EM | 10 + 7 | Ferman-Pinto imperfect-fit line (S10); Roth pretest critique (S11); classical DiD pre-trend testing literature identified as C3 incumbents |

## 3. Citation audit

| Item | Status | Resolution |
|---|---|---|
| Dossier ID 2603.24833 | RESOLVED | Agarwal-Choi-Yuan, robust matrix estimation with side information; register S14 |
| Dossier ID 2605.30319 | RESOLVED | Mehrotra-Tran-Vu-Zampetakis per-row MC HTE guarantees; register S08; deep-read scheduled Phase B entry |
| Agarwal-Dahleh-Sarkar (2019) | RESOLVED | "A Marketplace for Data" (1805.08125); irrelevant to causal panels; dropped from lineage (register S15) |
| Moon-Weidner UNVERIFIED flag | CLEARED | Two entries verified E2 (register S04-S05) |
| Theory-tool DOIs | ANCHORED | 13 of 15 verified via Crossref direct lookup; BBP real-case DOI unresolved by Crossref (flagged in register T1); Abadie 2021 JEL + Onatski-Moreira-Hallin left E1 pending |

## 4. Surfaces not closable programmatically (manual checklist, bounded task)

Google Scholar: blocked to automated query from this environment. Manual strings to run (10 min):
1. "synthetic control" "Marchenko-Pastur"
2. "synthetic control" "Tracy-Widom"
3. "recoverability" "synthetic control" eigenvalue
4. "BBP threshold" panel counterfactual
5. "spiked" "difference-in-differences"
6. "distance to frontier" synthetic control spectral

RePEc/IDEAS: search endpoint fetched but result markers absent (page appears script-rendered). Same six strings on https://ideas.repec.org and EconPapers.

Per plan rule (Section 4): the G1 decision does not rest on arXiv zero-hits alone; it rests on arXiv zero-hits PLUS full-text guarantee reads of the three closest works (A2) PLUS the abstract-level sweep above. The Scholar/RePEc manual pass is a residual formality logged as watch item W-3.

## 5. Log verdict

No direct hit on (a) per-unit BBP-type recoverability threshold, (b) TW-calibrated SCM diagnostics, or (c) exact high-dimensional limits in any queried vocabulary. Two watch items recorded (W-1, W-2). G1 evidence collection complete up to the programmatic limit of this environment.
