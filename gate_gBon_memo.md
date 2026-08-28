# Gate G-Bon (Phase E) — Bonander Florida SYG Applied Panel

**Date:** 2026-08-28. **Authority:** `preregistration_bonander_addendum.md` (frozen 2026-08-28, 70001-70004), `preprocessing_frozen.md` §8 leakage, `seeds.yaml` `phase_e.wp_bonander`. **Deterministic driver:** `scripts/wp_bonander_analysis.py` (`MKL 1`).

## 1. Panel

Bonander et al. (2021 AJE) `16×192` monthly `1999-01:2014-12` `syg_data.csv` `sha256 7f8bd93b`. Florida `SYG Oct 1 2005 time82` `T0=81` `n=15` `c=0.185`. Primary `HomicideRates` per 100k monthly (total homicide). Secondaries: `p100khomicide_firearm`, `p100khomicide_exclfirearm`, `Firearm.Suicide.Rates`, `Homicide_count`. Donor pool = 15 states listed in addendum §1; no NA. Provenance E2 CC-BY. BBP scaling fix `lam/sigma2` applied before freeze (addendum D_B1).

## 2. Treated (primary, `HomicideRates`, `code/applications/loaders.py:load_bonander`)

| stat | value |
|---|---|
| `k` | 2 |
| `d` | 23.56 `CI[15.88,126.54]` donor `q05>=1` |
| `p_align` | 0.0040 `r7.44` |
| `label` | **RECOVERABLE** |
| stability | 2/4 FAIL (`gate_lrv` silences) |
| `P3` date alarms | `[27,54]` `p0.896,0.657` `alarm_rate0.0` PASS |

## 3. Placebo battery (15 donors, `B500/200`)

All 15 placebos `FRAGILE-MISALIGNED` (`k1-3 d13.0-40.1 p0.13-0.97`). No `k0`. `d` entirely `>>1`, so channel-1 silent vs channel-2 dead. `placebo_label_counts {"FRAGILE-MISALIGNED":15}`.

Spearman `P1`: `rho=0.529 perm_p=0.978` vs `thr -0.30` **FAIL** (positive, not negative). `P2`: `d_treated 23.56 < Q90 30.33` **FAIL**. Hence `NF` is **only** `N3`: 1 comparable placebo (`New York rmse 0.056` vs treated `0.048` within 20%) classifies opposite (`RECOVERABLE` vs `FRAGILE`) → `N3 true`. `N1/N2` false.

## 4. Secondaries (descriptive, not gated)

| outcome | `k` | `d` | `p` | label |
|---|---|---|---|---|
| `p100khomicide_firearm` | 3 | 35.36 | 0.289 | FRAGILE-MISALIGNED |
| `p100khomicide_exclfirearm` | 1 | 30.56 | 0.489 | FRAGILE-MISALIGNED |
| `Firearm.Suicide.Rates` | 1 | 21.12 | 0.888 | FRAGILE-MISALIGNED |
| `Homicide_count` | 2 | 17.79 | 0.0005 | RECOVERABLE |

Firearm-specific outcomes are misaligned while total homicide (rate and count) is recoverable — outcome-definition heterogeneity; channel-2 flags.

## 5. Gate decision

* **Identification control:** PASS (`Florida` not `FRAGILE-INVISIBLE`).
* **`NF`:** true via `N3` certification asymmetry (first time all placebos misaligned while treated recoverable).
* **`P2/P1`:** both FAIL; `stability` FAIL. Therefore **applied-value `PASS` = FALSE** (requires `NF & P2 & P3 & stability`). `verdict_full` false. This is `certification-with-asymmetry` but not a `FULL PASS`.

Interpretation: Bonander is a **moderate-supercritical certification panel** (`d~23` similar to `Basque 5.8` regime, far from `CA 2856/Germany 4.8e8/Freire 199` saturation) that demonstrates the pipeline's **channel-2 discrimination**: donor factor `k2` exists and Florida loads on it (`p0.004`), donors' placebos do not. It does **not** deliver a `d<1` or `FRAGILE-INVISIBLE` fragile treated case; `N1` remains unfired. No `KILL` trigger; `P3` null model re-validated (`alarm_rate 0`).

## 6. What remains for a fragile panel

The hunt for a `d<1` / `k0` / `p>=0.05`-treated fragile case continues. Bonander shows `misaligned placebos` are easy to find when `n_d=15` monthly homicide; a `misaligned treated` with small `n_d/T0≈1` and noisy micro outcome (school, patent, syphoned health) is still needed to fire `N1`. Next candidate remains the public NBER `Wisconsin 8×161 handgun` proxy (`8 donors, 161 mo, c0.05`) or a staggered RTC early-adopter slice.

## 7. Artifacts

`results_e/bonander/distance_bonander.json`, `placebos_bonander.csv`, `summary_bonander.json`, `figures/fig_bonander_distance_to_frontier.png`.

## 8. Disposition

`PROCEED` to next Phase E candidate with same frozen thresholds; Bonander ships as **certification case study with N3 asymmetry** alongside `CA/Germany/Basque/Freire`. Update `Idea5_Panel_RMT_Research_Plan.md` Phase E: `Bonander primary RECOVERABLE (N3) but P2/P1/stability fail → INCREMENTAL-ONLY`.
