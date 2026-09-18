# Preregistration FBI RTC Addendum: Nebraska Monthly Violent Crime Candidate Screen

**Date frozen:** 2026-08-29, after data acquisition and a coverage-only feasibility audit, but before any Nebraska or donor outcome spectrum, alignment statistic, pre-fit RMSE, or frontier label was computed. **Purpose:** determine whether one public, policy-aligned micro panel qualifies for a full trusted-result reproduction. This is a screen, not yet a Gate G4 pass.

**Parent documents:** `docs/preregistrations/preregistration_d2_addendum.md`, `docs/preregistrations/preregistration_bonander_addendum.md`, `docs/gates/gate_g4_memo.md`, `docs/gates/gate_g4_wave2_memo.md`, and `docs/gates/gate_gBon_memo.md`. All frontier thresholds are inherited unchanged.

## 1. Candidate selection without outcome leakage

Donohue, Aneja, and Weber (2019) study the effect of right-to-carry (RTC) laws on violent crime with state-level synthetic controls. Their paper states that Kansas and Nebraska RTC laws took effect at the beginning of 2007 and that the donor design uses states with no RTC law as of 2014, plus sufficiently late adopters when available. The FBI CDE extract supplies monthly violent-crime rates from January 1999 onward.

Before inspecting any outcome values, candidate states and donors were screened solely on the FBI's monthly `Percent of Population Coverage` field over January 1999 through December 2006. The frozen floor is 90% in every pre-month. Kansas fails with a minimum of 17.31% and is excluded. Nebraska passes with a minimum of 91.29% and is selected. The nine no-RTC donor candidates are California, Connecticut, Delaware, Hawaii, Maryland, Massachusetts, New Jersey, New York, and Rhode Island. New York fails with a minimum of 78.21%; the other eight pass and form the frozen donor pool.

This yields `9` units by `228` months (`1999-01:2017-12`), with Nebraska treated at `2007-01`, `T0=96`, `n_d=8`, and `c=8/96=0.083333`. Donor selection may not be revised after the spectrum is observed.

## 2. Data and provenance

**Source:** FBI Crime Data Explorer, summarized state endpoint, `violent-crime`, `type=estimates`. Raw artifacts are `data/raw/fbi_cde/fbi_cde_monthly_1999_2017.csv` (SHA-256 `f6053c3d591a4b39cba3558ec602f9b755c03e8b351de474b9532ebe571f5684`) and its request ledger `fbi_cde_monthly_1999_2017.meta.json` (SHA-256 `04fec2f22a0d0ffdbb9b0623b42764e12eff0fe166309ad21be042debd043b48`). The acquisition is documented in `data/raw/fbi_cde/README.md`.

**Primary gated outcome:** monthly estimated violent-crime rate per 100,000 as returned under the API's `rates` branch. **Registered scale audit:** monthly estimated count under the API's `actuals` branch, analyzed for Nebraska only and never used for the screen decision.

The CDE rate remains subject to agency-participation composition. The 90% floor limits, but does not eliminate, this threat. Coverage trajectories and a 95% floor sensitivity will be reported before any full application is authorized.

## 3. Frozen frontier screen

Use unit-centered pre-period levels and `applications.frontier_d2.analyze_unit` without modification. Inherit `GATE_TOL=1.05`, `K_MAX=4`, `G_NULL=2000`, alignment threshold `p<0.05`, donor bootstrap `B=500`, circular time-block bootstrap `B=200` with block length `4`, and four-way selector/window stability. BBP inversion uses `lambda/sigma2`.

Analyze Nebraska once against the eight frozen donors. Analyze each donor once as pseudo-treated against the other eight units, including Nebraska, using pre-2007 data only. Record `k`, `d`, alignment, donor and time intervals, label, stability, and centered OLS pre-fit RMSE.

P1 is descriptive because eight placebos do not support the inherited `n=16` or `n=38` thresholds. Report Spearman rho and its 999-permutation one-sided p-value without a pass label. P2 retains the frozen rule `d_Nebraska >= Q90(d_placebo)`. P3 uses `z_shift` at pre-window indices `36` and `60` (January 2002 and January 2004 boundaries), with alarm rate at most `0.20`. Stability retains agreement in at least three of four selector/window combinations.

The inherited NF arms are computed on the primary rate. `N1` is a Nebraska `FRAGILE*` label. `N2` requires Nebraska to be `RECOVERABLE` and its pre-fit RMSE to exceed the placebo 75th percentile. `N3` uses the frozen 20% comparable-fit band and opposite certification labels.

## 4. Screen decision and stop rule

`QUALIFIED_FOR_TRUSTED_REPLICATION` requires an NF arm, P2, P3, and stability. This authorizes a separate, preregistered reproduction of the Donohue, Aneja, and Weber Nebraska specification and its application-standard baselines. It does not itself pass Gate G4.

If Nebraska is `FRAGILE` but any control above fails, record `FRAGILE-BUT-SCREEN-FAILED`; the panel cannot anchor the paper without a repair justified independently of its outcome. If Nebraska is `RECOVERABLE` and no NF fires, record `DISCARD-NO-FRAGILE-TREATED` and stop this branch. If coverage composition or missingness invalidates the panel, record `DISCARD-DATA-QUALITY`.

No donor substitution, alternative cutoff, offense switch, coverage-floor relaxation, or window tuning is permitted after the primary result. Any such change is a new exploratory panel and cannot inherit this screen's status.

## 5. Seeds, artifacts, and compute

Seeds in `config/seeds.yaml` under `phase_e.wp_fbi_rtc_nebraska` are alignment null `72001`, donor bootstrap `72002`, time bootstrap `72003`, and Spearman permutation `72004`, with `base + index*10^6` offsets. `z_shift` and post-test pools retain `8880001` and `7770001`. The decisive artifact is `results_e/fbi_rtc_nebraska/summary_fbi_rtc_nebraska.json`. Expected runtime is minutes locally with negligible RAM.

## 6. References

Donohue, John J., Abhay Aneja, and Kyle D. Weber. 2019. “Right-to-Carry Laws and Violent Crime: A Comprehensive Assessment Using Panel Data and a State-Level Synthetic Control Analysis.” *Journal of Empirical Legal Studies* 16(2): 198-247. https://doi.org/10.1111/jels.12219

Federal Bureau of Investigation. *Crime Data Explorer*. https://cde.ucr.cjis.gov/LATEST/
