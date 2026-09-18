# Gate G-bi63 Memo: Wisconsin Yearly Suicide Proxy Closure

**Date:** 2026-08-29. **Authority:** `docs/preregistrations/preregistration_bi63_wisconsin_addendum.md`, frozen after full disclosure of the exploratory point probes and before the reproducibility driver. **Execution:** `scripts/wp_bi63_analysis.py`, seeds `71001:71004`. **Artifact:** `results_e/bi63/summary_bi63.json`.

## 1. Panel and scope

The CDC/NCHS bi63 extract supplies annual all-suicide deaths and age-adjusted rates for `1999:2017`. Wisconsin is paired with the eight other waiting-period states in Powell's replication code. The split is `1999:2014` pre and `2015:2017` post, so `T0=16`, `n_d=8`, and `c=0.5`.

This is not Powell's restricted monthly handgun-suicide outcome. The annual all-suicide proxy cannot resolve a June 2015 policy change and cannot identify the paper's policy estimand. It was preregistered only to reproduce and close a disclosed certification probe.

## 2. Primary age-adjusted rate

| Statistic | Wisconsin result |
|---|---|
| `k` | `1` |
| `d` | `23.073` |
| donor interval | `[22.556, 192.948]` |
| `p_align` | `0.0005` |
| label | `RECOVERABLE` |
| time-block interval | `[18.937, 41.782]` |
| stability | `2/4`, FAIL |

All eight placebo states are `RECOVERABLE`, with point `d` between `24.256` and `40.009`. The primary selector therefore saturates comfortably above the frontier and supplies no treated-placebo separation. The robust `gate_lrv` selector silences both full and trimmed windows, so Wisconsin and all placebo analyses fail the inherited `3/4` stability criterion.

## 3. Count-scale audit

Wisconsin deaths classify `RECOVERABLE` with `k=1`, `d=514.540`, donor interval `[53.514, 4040.860]`, and `p_align=0.0005`. All eight placebo states are also `RECOVERABLE`, with point `d` between `67.720` and `520.628`. Stability again fails `2/4` because `gate_lrv` selects `k=0`.

Counts magnify the common population-scale factor and move the panel farther into the saturated regime. They do not provide a substantively preferable policy outcome.

## 4. Decision

**Mechanical result:** PASS. The official raw data, pinned checksum, loader, registered seeds, treated run, two outcome scales, and sixteen placebo runs reproduce deterministically; all targeted loader/frontier tests pass.

**Applied result:** `INCREMENTAL-ONLY`, certification closure. The primary rate and count audit are uniformly deep-supercritical under the primary selector, show no treated-placebo asymmetry, fail robust-selector stability, and do not reproduce the restricted policy outcome. No NF, P1, P2, P3, or Gate G4 claim is authorized.

**Consequence:** close bi63 permanently as a paper-side scope illustration at most. Do not spend further effort polishing it into a causal case study.
