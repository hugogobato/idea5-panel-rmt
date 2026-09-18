# Gate G-FBI-RTC Memo: Nebraska Monthly Violent Crime Candidate Screen

**Date:** 2026-08-29. **Authority:** `docs/preregistrations/preregistration_fbi_rtc_nebraska_addendum.md`, frozen after a coverage-only audit and before any outcome spectrum. **Execution:** `scripts/wp_fbi_rtc_nebraska_screen.py`, seeds `72001:72004`. **Artifact:** `results_e/fbi_rtc_nebraska/summary_fbi_rtc_nebraska.json`.

## 1. Leakage-safe candidate construction

The official FBI CDE extract contains 50 states by 228 months (`1999-01:2017-12`) with no missing count, rate, population, or coverage field. Candidate selection used only the pre-2007 reporting-coverage field. Kansas was rejected because its minimum pre-period coverage is `17.31%`. Nebraska passes at `91.29%`. Of the nine no-RTC donor candidates, New York fails at `78.21%`; California, Connecticut, Delaware, Hawaii, Maryland, Massachusetts, New Jersey, and Rhode Island pass the frozen `90%` floor.

The resulting panel is `9×228`, with Nebraska treated in `2007-01`, `T0=96`, `n_d=8`, and `c=0.08333`. This is the first policy-aligned, low-aspect monthly panel in the application search.

## 2. Preregistered primary result

| Statistic | Nebraska violent-crime rate |
|---|---|
| `k` | `1` |
| `d` | `154.891` |
| donor interval | `[73.990, 2636.109]` |
| `p_align` | `0.8141`, `r_align=0.126` |
| label | `FRAGILE-MISALIGNED` |
| stability | `2/4`, FAIL |
| P2 | `154.891 < Q90 363.905`, FAIL |
| P3 | alarm rate `0/2`, PASS |
| NF | `N1=true`; `N2=N3=false` |

This is the first preregistered treated `FRAGILE` result in the project, but it is not a successful applied screen. Seven of eight placebos are also `FRAGILE-MISALIGNED`; Delaware alone is `RECOVERABLE`. The descriptive fit ordering is absent (`rho=0.143`, permutation `p=0.669`), three placebos have comparable pre-fit, and none produces N3 asymmetry. Nebraska is therefore not separated from a generally misaligned panel.

The selector disagreement is systematic: the primary gate returns `k=1` for full and trimmed windows, while `gate_lrv` returns `k=0` for both. The donor-post factor-law control also rejects at `p=0.0033`, indicating that even a pre-period classification would not certify post-2007 transport.

## 3. Registered count-scale audit

The estimated monthly count flips Nebraska to `RECOVERABLE` (`k=1`, `d=172.942`, `p_align=0.0005`) while retaining the same `2/4` stability failure. This scale dependence is scientifically important: population scale supplies alignment that is absent from the rate panel. It also prevents the rate result from being framed as a stable property of the policy design.

## 4. Screen decision

The frozen qualification rule requires an NF arm, P2, P3, and stability. Nebraska has N1 and P3 but fails P2 and stability. The decision is therefore:

**`FRAGILE-BUT-SCREEN-FAILED`.**

No trusted-result reproduction is authorized from this screen, and Gate G4 remains failed. Nebraska may appear in a diagnostic casebook as the first real treated misalignment phenotype, provided the paper states that the phenotype is panel-wide, selector-unstable, scale-dependent, and accompanied by post-window drift. It cannot be the flagship applied study.
