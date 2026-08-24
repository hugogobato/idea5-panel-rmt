# Theory Targets (WP-E1 stub, populated incrementally from Phase B)

**Rule (plan Section 10):** no substantial proof work before Gates G3/G4 pass. This file is the mapping layer between surviving evidence and future proof packages. Each target must trace to a verified register entry (E-level >= E2 with anchors).

| ID | Target | Serves | Source results to lean on (register IDs) | Adaptation gap | Numerical falsifier | Status |
|---|---|---|---|---|---|---|
| T1 | DE risk formula for hard-threshold spectral SC + ridge-SC; boundary where risk stops improving (frontier_ansatz.md made rigorous) | C1 upper side, C2 engine | BBP 2005 (T1); BGN 2011 outlier locations AND overlaps (T2); Dobriban-Wager 2018 ridge risk (T6); Hastie et al. 2022 ridgeless (T7) | Row-targeted single-unit risk instead of Frobenius; resolvent variant for treated-row projection; post-window transport term | WP-C1 overlay: F vs simulated RMSE to plotting accuracy across c-grid | ANSATZ STATED (Phase B, conjecture tag); proof dormant until G3 |
| T2 | Matching lower bound: sub-threshold or zero-leverage units unrecoverable by any estimator at o(1) risk | C1 lower side (paper spine) | Spiked MC minimax templates (Koltchinskii-Lounici line; Birnbaum et al.); Farias-Li-Peng Prop 2 tangency template (S02); note Agarwal et al. (S03) contains NO lower bounds (gap confirmed full-text) | Missing pattern is structured block-row, not MCAR; Le Cam two-point between treated-loading configurations separated below noise floor | Sub-threshold witness at scale: no implemented estimator crosses floor (WP-C1 sub-edge cells) | CONJECTURE; witness-level support only (Phase A) |
| T3 | TW/Jacobi calibration of rank selector + spectral pre-trends test under iid; bootstrap fallback under weak dependence | C3 | Johnstone 2001 (T3), 2008 Jacobi (T4); Onatski 2010 ratios (T5a), 2009 (T5b); universality under weak dependence (Bao-Pan-Zhou line, to anchor) | Serial correlation breaks exact nulls; whitened largest-eigenvalue null for residual projections | WP-C2(iii) size plots; onset-convergence prediction inherited from Witness 1 (onset -> m=1 as n grows) | DORMANT until G3 |
| T4 | SDID deterministic equivalents | C2 zoo comparison | Arkhangelsky et al. 2021 definition (T9); DW machinery (T6) | Coupled weight systems lack closed forms | DE-vs-simulation overlay if attempted | CUT unless reviewers demand (anti-decoration rule) |

New tools registered during Phase B entry reads:

- Mehrotra-Tran-Vu-Zampetakis (S08), Theorem B.2: sharp row-wise ||.||_{2,infinity} truncated-SVD perturbation bound. Candidate auxiliary input for T1's row-targeted error control and possibly T2's converse-side bookkeeping.
- Spiess-Imbens-Venugopal (S07): model-averaging mechanics; positioning only (W-2 text frozen in register), no theory target.
