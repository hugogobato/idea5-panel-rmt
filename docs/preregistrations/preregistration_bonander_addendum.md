# Preregistration Bonander Addendum — Florida Stand-Your-Ground (Bonander et al. 2021 AJE) as Phase E applied panel

**Date frozen:** 2026-08-28, BEFORE any decisive Bonander computation beyond the single exploratory Florida HomicideRates probe (2026-08-28) that forced the BBP scaling fix (`lam/sigma2`) and revealed `k=2 d~23` saturation. No placebo `P1/P2`, NF, or gate decision had been computed. **Parent documents:** `docs/preregistrations/preregistration.md` (2026-08-24), `docs/preregistrations/preregistration_c5_addendum.md` (2026-08-25), `docs/preregistrations/preregistration_d2_addendum.md` (2026-08-26, thresholds inherited verbatim), `docs/preregistrations/preregistration_d2b_addendum.md` (2026-08-26), `docs/applied/preprocessing_frozen.md` (leakage rules bind), `docs/gates/gate_g3_memo.md` §7 (restricted claims + G3 GO), `docs/gates/gate_g4_memo.md` + `docs/gates/gate_g4_wave2_memo.md` (Phase D closed, Phase E gate open under restricted claim set). **Authority:** plan §7 Phase E: evidence-earned paper consolidation; this addendum is the first Phase E decisive applied run. Every threshold below is *inherited* from the already-frozen Phase C/D instruments or fixed from theory; none was tuned toward Bonander numbers.

---

## 1. Panel and provenance (E2, open)

**Source:** Bonander, Humphreys & Degli Esposti (2021) *American Journal of Epidemiology* `10.1093/aje/kwab211`; replication archive `rvayc-osfstorage-archive.zip` → `syg_data.csv` `689544 B` `sha256 7f8bd93b6add9e5e1b29d3ae738d30e90c17a2e1c867ebdee3c411f0341095b9` pinned in `data/raw/bonander/syg_data.csv` and `SHA256SUMS`, `CC-BY`. The file was downloaded from the OSF project referenced in the paper's data-availability statement via the `rvayc` OSF storage archive URL (now 500).

**Panel:** `16` US states × `192` months `1999-01:2014-12` (`time 1:192`). The `16` are the paper's analytic sample: `Arkansas (5), Connecticut (9), Delaware (10), Florida (12 - treated), Hawaii (15), Iowa (19), Maine (23), Maryland (24), Massachusetts (25), Nebraska (31), New Jersey (34), New York (36), North Dakota (38), Ohio (39), Rhode Island (44), Wyoming (56)` (FIPS-derived `State.Code`). Verified: `pivot State x time` on `HomicideRates` is `16×192` with `0` NA. `Florida` is `Case==1` (`Case 0:2880 1:192`); `treatdummy 0/1` flips at `Year 2005 Month 10 time 82`.

**Treatment:** Florida Stand-Your-Ground `Oct 1 2005`. `treat_year 2005`, `treat_month 10`, `treat_time 82` (`(2005-1999)*12+10`). `pre = time 1:81` (`1999-01:2005-09`) `T0=81`. `post = 82:192` (`2005-10:2014-12`) `Tp=111`. Donor pool = `15` states above minus Florida, `c=n_d/T0=15/81≈0.185`. No supermajority exclusion; national aggregate not present.

**Outcomes (all monthly `p100k`, pre-centered levels):**
* Primary (gated): `HomicideRates` = total monthly homicide per 100k (the paper's headline `HomicideRates`, monthly `0-1.83` mean `0.44`).
* Registered secondaries (descriptive, never gated): `p100khomicide_firearm` (firearm homicide), `p100khomicide_exclfirearm` (non-firearm homicide), `Firearm.Suicide.Rates`, `Homicide_count` (raw count) — to probe heterogeneity; the gate decision rests *only* on the primary.

## 2. Analysis space and leakage

Primary space: **unit-centered pre-period levels** (each state's pre row minus its own `T0` mean). Rationale identical to D2 §2 (MP bulk, BBP). No common-trend removal. Secondary descriptive space: **first differences** (`T0-1`, each row self-centered) — reported as interpretation, never gated.

Leakage: every spectral quantity uses donor `PRE` + treated `PRE` only. Donor `POST` may enter only through the two shipped instruments `z_shift / pre_trends_post_test`; treated `POST` never enters any statistic except for the incumbent `postgap` description. All panels via `applications.loaders.load_bonander`.

## 3. Channel 1 & 2 (inherited verbatim from D2 §3-4, with 2026-08-28 scaling fix)

* `evals` = eigenvalues of `(1/T0) Yc Yc'` (`Yc` donor centered). `c=n_d/T0`.
* `sigma2_0 = median(evals)/mu_MP(c)`, `mu_MP` = MP median (unit-mean law). `k0 = gate(evals, sigma2_0)` (`0` if `lam1 <=1.05*sigma2_0*(1+sqrt(c))^2`, else largest-gap among top 4).
* Bulk refinement: if `n_d-k0 >= max(10, ceil(n_d/4))` then `sigma2_1 = mean(evals[k0:])` else `sigma2_1=sigma2_0`; floor at `1e-8*lam1`.
* Final `k = gate(evals, sigma2_1)`, `K_MAX=4`.
* **Scaling fix (2026-08-28, before this addendum freeze):** BBP inversion is on `sigma^2=1` scale: `lam_scaled = evals[:k]/sigma2_1`, `s = invert_bbp(lam_scaled,c)` with `lam=1+s+c+c/s`, `m=s/sqrt(c)`. `d = max m` if `k>=1` else `0`. The fix was forced by the single exploratory Florida run where raw `lam~0.35` gave `nan` vs scaled `lam/sigma2~11` gives `m~23`. The fix is theory-mandated and applies to all panels; prior D2/Wave2 panels were supercritical and remain so (magnitudes rescaled but labels unchanged).

Alignment: `e_obs = sum_{j<=k} <y1c, v_j>^2 / ||y1c||^2`, `G=2000` Gaussian rows `~N(0,I_T0)` against *fixed* `Vt`, `p_align=(#{e>=e_obs}+1)/(G+1)`, `r_align=e_obs/median(null)`. Channel-2 confirmed iff `p<0.05`.

## 4. Classification (frozen)

| Label | Rule |
|---|---|
| FRAGILE-INVISIBLE | `k==0` |
| RECOVERABLE | `k>=1 & p<0.05 & donor q05(d)>=1` |
| INCONCLUSIVE-SUBEDGE | `k>=1 & p<0.05 & CI straddles 1` |
| FRAGILE-SUBEDGE | `k>=1 & p<0.05 & q95<1` |
| FRAGILE-MISALIGNED | `k>=1 & p>=0.05` |

## 5. Uncertainty

`B_donor=500` rows-with-replacement donor resamples (full pipeline per draw, `donor_bootstrap_seed + idx*1e6`), `B_time=200` circular moving-block `block 4` time resamples. Percentile `95%` CI for `d`. Both layers for treated and every placebo.

## 6. Placebo battery and separation criteria (frozen, one-sided)

In-space placebos: each donor state `i` as pseudo-treated, donor pool = `all 16 states except {i}` (Florida stays as donor when `i != Florida`; standard contamination guard drops only the placebo's own row). Record per placebo: donor `d_i, p_i, label_i` and a simple pre-fit `rmse_i` (proxy: OLS `w=argmin||y1 - w'Yd||_2` `rmse=sqrt(mean((y1 - w'Yd)^2))` on the centered pre window; the incumbent simplex `SCM rmse` is also reported descriptively). The simplex search, if used, uses the single frozen `V` seed from `config/seeds.yaml` (structure held constant across placebos).

Separation (frozen thresholds, inherited):
* `P1` (fit-ordering): Spearman `rho(d_i, rmse_i)` across `15` placebos `<= -0.30` (smoking-like `n≈15` critical) with permutation `p<=0.05` (`999` reps, seed `perm`); one-sided `H1: rho<0`.
* `P2` (treated separation): `d_treated >= Q90({d_i})` (placebo `d` distribution).
* `P3` (date alarms): `z_shift` (shipped C5c, defaults `B=200`) on donor-pre with pseudo-cutoffs `time {27,54}` (≈ `2001-03, 2003-06` inside the pre) — alarm rate `p<0.05` must be `<=0.20`.
* Stability (gated): `>=3` of `4` combos `selector {primary gate, gate_lrv} x window {full 1:81, trimmed 3:79 (drop 2+2 months)}` agree on label (primary space only). `ungated gap` and `cv_rank` reported.

Descriptive control: `pre_trends_post_test` donor-pre basis vs real donor-post window, simulated iid finite-n null (`G=300`).

## 7. Novel-finding arms and verdict (frozen, inherited from D2 §10)

NF arms (any one suffices for a novel finding):
* `N1`: treated classifies `FRAGILE*` with frontier explaining observed pre-fit structure that informal inspection could not quantify.
* `N2`: poor-looking treated (`rmse_treated > Q75(placebo rmse)`) certified `RECOVERABLE`.
* `N3`: certification asymmetry: placebo units with `comparable` pre-fit (`|rmse - rmse_treated| <=0.20*|rmse_treated|`) classify differently (`treated RECOVERABLE` while comparable placebos `FRAGILE` or vice versa).

Verdict (frozen):
* Identification control: `Florida` may not classify `FRAGILE-INVISIBLE` (positive-control failure → investigate, KILL candidate per plan).
* Applied-value `PASS` requires: `NF` on the primary outcome **and** `P2 PASS` on that panel **and** `P3 PASS` **and** Section 6 stability on the NF panel. If additionally `P1 PASS`, verdict `FULL PASS`; otherwise `PASS-WITH-LIMITATION`.
* If primary classifies comfortably `RECOVERABLE` and no NF fires: `certification-only` → `INCREMENTAL-ONLY` pressure (plan give-up 4).
* If no `P2/P3` separation anywhere: give-up 2 territory.

## 8. Seeds, environment, compute

Seeds registered in `config/seeds.yaml` `phase_e.wp_bonander`: `alignment_null 70001`, `donor_bootstrap 70002`, `time_bootstrap 70003`, `spearman_permutation 70004`, `z_shift 8880001` (C5 convention), `posttest_null 7770001`, `placebo_V 60001` (identical across placebos). Offsets `base + idx*1e6` (Florida `idx 0`, placebos `1..15` ordered alphabetically). Environment inherits `docs/applied/preprocessing_frozen.md` §9; compute local only, minutes.

## 9. Deviation log

`D_B1` (pre-freeze fix, 2026-08-28): BBP scaling `lam/sigma2` mandated before freezing this addendum; exploratory Florida probe that forced it is discounted (not a decisive placebo gate).

`D_B2` (scope, pre-freeze): `RHO_THRESH` for `P1` inherited as `-0.30` (smoking-like `n=16`); the Germany `-0.42` threshold is not imported.

No post-freeze deviation yet.
