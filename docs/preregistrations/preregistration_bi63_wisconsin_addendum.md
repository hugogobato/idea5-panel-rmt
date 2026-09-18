# Preregistration bi63 Addendum: Wisconsin Yearly Suicide Proxy (Certification Only)

**Date frozen:** 2026-08-29, before the reproducibility run written by `scripts/wp_bi63_analysis.py`. This addendum is not blinded. The handoff disclosed exploratory point probes on the same data (`Deaths`: approximately `k1`, `d514`, `p0.0005`; age-adjusted rate: approximately `d23`) and disclosed that all eight placebo states looked `RECOVERABLE`. Those probes cannot support a new-finding claim. The only authorized purpose of this run is to reproduce, document, and close the proxy as a certification-only case.

**Parent documents:** `docs/preregistrations/preregistration_d2_addendum.md` (frozen frontier instrument), `docs/preregistrations/preregistration_bonander_addendum.md` (theory-mandated `lambda/sigma2` BBP scaling fix), `docs/gates/gate_g4_memo.md`, `docs/gates/gate_g4_wave2_memo.md`, and `docs/gates/gate_gBon_memo.md`. All numerical thresholds are inherited without tuning.

## 1. Scientific scope and non-claim

The motivating policy is Wisconsin's repeal of its 48-hour handgun purchase waiting period in June 2015, studied by Powell's *Imperfect Synthetic Controls*. The public Powell replication archive identifies the nine-state waiting-period sample as Wisconsin plus California, Hawaii, Indiana, Iowa, Maryland, Mississippi, New Jersey, and Rhode Island. Its monthly handgun-suicide mortality outcome is restricted. The public CDC bi63 table instead contains annual all-suicide deaths and age-adjusted all-suicide death rates.

Consequently, this panel does not reproduce Powell's outcome, time resolution, or estimand. No causal effect of the June 2015 repeal will be estimated or inferred. The panel tests only whether the already-frozen frontier instrument certifies a small-donor annual mortality panel as spectrally recoverable. Regardless of its label, this run cannot rescue Gate G4 and cannot be called a proper flagship applied study.

## 2. Data, provenance, and frozen panel

**Source:** CDC/NCHS, `NCHS - Leading Causes of Death: United States`, Socrata identifier `bi63-dtpu`, public-domain United States government work. Raw file: `data/raw/cdc_bi63_all_1999_2017.csv`, `10,868` data rows plus header, SHA-256 `76458c5abfc7f959fceecf6661fedd36418c4c2d9f2ab297c26910b7ab3b4b42`. The metadata state that bi63 contains age-adjusted death rates for leading causes, based on resident death certificates, with rates per 100,000 standardized to the 2000 United States population.

**Cause filter:** exact `113 Cause Name` equal to `Intentional self-harm (suicide) (*U03,X60-X84,Y87.0)`.

**Units:** the nine Powell waiting-period states listed above. Wisconsin is treated. The other eight states are donors.

**Time:** `1999:2017`. The pre-period is `1999:2014`, so `T0=16`, `n_d=8`, and `c=0.5`. The annual treatment index is `2015`. Because the repeal occurred during June, every 2015 observation is post for the frontier split, but no post-treatment effect is estimated.

**Outcomes:** the primary certification series is `aadr`, the age-adjusted all-suicide death rate per 100,000. The registered scale audit is `deaths`, the annual all-suicide count. Counts are expected to be more strongly dominated by population scale and are not substantively preferred.

## 3. Frozen frontier analysis

For each outcome, use unit-centered pre-period levels and `applications.frontier_d2.analyze_unit` without alteration. The instrument inherits `GATE_TOL=1.05`, `K_MAX=4`, `G_NULL=2000`, alignment significance `p<0.05`, donor bootstrap `B=500`, circular time-block bootstrap `B=200` with block length `4`, and trimmed-window sensitivity dropping two leading and two trailing pre-years. BBP inversion uses `lambda/sigma2` as frozen on 2026-08-28.

The five classification labels and their exact rules are inherited from `docs/preregistrations/preregistration_d2_addendum.md` Section 5. No new threshold is introduced.

## 4. Placebos, outputs, and closure rule

Each of the eight donor states is analyzed once as pseudo-treated against the other eight states, including Wisconsin. Only pre-2015 observations enter all spectral quantities. Record `k`, `d`, donor-bootstrap interval, `p_align`, label, stability, and centered OLS pre-fit RMSE for Wisconsin and every placebo.

The decisive artifact is `results_e/bi63/summary_bi63.json`. Mechanical success requires both outcomes, one Wisconsin result and eight placebo results per outcome, finite diagnostics, and loader tests passing.

The scientific disposition is fixed in advance:

1. If the primary age-adjusted rate and its placebos are comfortably `RECOVERABLE`, record a saturated certification case and close it as `INCREMENTAL-ONLY` for applied value.
2. If Wisconsin is `FRAGILE` or differs from comparable placebos, record the anomaly as exploratory-only because the result was already outcome-inspected and the public proxy does not identify the policy estimand. It may motivate an independently preregistered panel, but it does not fire an NF arm here.
3. If the loader, provenance, or classification is unstable, close the proxy as unusable.

No `P1`, `P2`, `P3`, NF, or Gate G4 pass is computed. This prevents a known, non-blinded proxy from being promoted after the fact.

## 5. Seeds and compute

Seeds are registered in `config/seeds.yaml` under `phase_e.wp_bi63`: alignment null `71001`, donor bootstrap `71002`, time bootstrap `71003`, and the reserved placebo-order/permutation seed `71004`. Standard offsets are `base + index*10^6`, with Wisconsin index `0` and placebo states indexed `1:8` alphabetically. The run is local, single-process, and expected to require minutes and negligible RAM.

## 6. References

CDC/NCHS. *NCHS - Leading Causes of Death: United States*. https://data.cdc.gov/d/bi63-dtpu

Powell, D. *Imperfect Synthetic Controls*. The local replication archive documents the restricted mortality files and public population resource. The paper DOI or stable publisher URL must be added to the paper bibliography when the certification vignette is written.
