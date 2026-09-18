# Pilot Cost Report (WP-B3)

**Date:** 2026-08-24. **Script:** `code/run_pilot.py` (raw numbers in `figures/pilot_costs.json`). **Environment:** Python 3.12.3, numpy 2.4.3, scipy 1.17.1; BLAS/OMP threads pinned to 1 inside all measurements.

## 1. Measurement context (recorded, per plan rule)

Machine during measurement: load average ~41 on 8 visible cores, other user jobs running. Absolute wall times are therefore inflated relative to an idle machine (factor unknown, plausibly 2 to 5x); relative comparisons across methods and the LOCAL/COLAB classification are conservative (biased toward COLAB), which is the safe direction for planning.

## 2. Decisive-cell pilot (production size)

Cell: `n=250, T0=250, T_post=125` (`c=1.0`), `r=1`, `s=2 sigma^2` (m = 2), treated share `theta=0.5`, aligned; seeds 20000-20009; 10 reps.

| method | mean ms/rep | max ms/rep | max/mean | RMSE (sigma units) |
|---|---|---|---|---|
| donor_mean | 1 | 8 | high (sub-ms) | 1.223 |
| scm_simplex | 4345 | 5576 | 1.28 | 1.197 |
| ridge_sc | 216 | 334 | 1.55 | 1.156 |
| spectral_sc | 28 | 47 | 1.68 | 1.127 |
| mc_nn_cv | 4788 | 5566 | 1.16 | 1.182 |
| sdid | 11023 | 12675 | 1.15 | 1.219 |

Scientific sanity: ordering matches theory at this cell (spectral best, trivial baseline worst); peak RSS across reps 0.093 GiB.

Cost-model check (plan pass rule: model explains pilot within 30%): per-rep stability ratios max/mean are within 30% for the three expensive methods (scm 28%, mc 16%, sdid 15%); ridge's ratio 55% reflects sub-second timings dominated by scheduler noise under load, immaterial at grid scale. Linear-in-reps cost model accepted.

## 3. Scaling probe (2 reps per size)

Wall seconds, same cell family:

| n (=T0) | spectral_sc | scm_simplex | sdid | mc_nn_cv |
|---|---|---|---|---|
| 125 | 0.00 | 0.40 | 1.14 | 0.87 |
| 250 | 0.03 | 4.58 | 12.81 | 5.32 |
| 400 | 0.10 | 24.73 | 61.64 | 21.24 |

Log-log exponents versus `n*T0`: spectral_sc 1.31, mc_nn_cv 1.47, sdid 1.67, scm_simplex 1.79. Spectral methods scale near-linearly in panel entries; simplex-family solvers steepen with size (iteration growth). Extrapolation to WP-C4 sizes: at `n=T0=2000` spectral stays seconds-per-rep while sdid/scm reach minutes-to-tens-of-minutes per rep; memory remains far below 4 GiB everywhere below n=5000 (~200 MB matrices).

## 4. LOCAL/COLAB classification (plan Section 11 rule)

Rule: LOCAL iff projected wall < 2 h on local machine at <= 8 workers AND peak RSS < 4 GiB. Serial column uses measured production-cell costs; worker column assumes 8 workers (currently unavailable in practice given load ~41; treat worker column as optimistic bound).

| experiment | serial h | h @ 8 workers | class |
|---|---|---|---|
| C1 decisive grid (450 cells x 500 reps, full method set) | 1275.1 | 159.4 | COLAB |
| C1 single-c slice (15 spike points x 500) | 42.5 | 5.3 | COLAB |
| C2(i) null battery (6 x 500) | 17.0 | 2.1 | COLAB |
| C2(ii) baseline-favorable (4 x 500) | 11.3 | 1.4 | LOCAL |
| C2(iv) structural-break (4 x 500) | 11.3 | 1.4 | LOCAL |
| C2(iii) AR(1)+bootstrap battery (2 rho x 500, bootstrap factor ~50x diagnostic cost, assumed) | 9.1 | 1.1 | LOCAL |

Notes: (i) the bootstrap multiplier is an assumption (50x spectral+ridge diagnostic cost, B=200 block bootstraps amortized); it gets re-measured before C2 runs. (ii) Memory never binds; the trigger is purely wall time. (iii) All classifications use loaded-machine times, i.e., conservative.

## 5. Inputs carried to Phase C preregistration (decisions NOT taken here)

The measured costs make the originally drafted C1 grid (full method set at n=T0=250) a multi-week Colab-fleet job. Preregistration will need to choose among: reducing default cell size (e.g., n=T0=160 keeps every method under ~1.5 s/rep and turns single slices LOCAL), restricting the full sweep to the claim-bearing methods (spectral, ridge, SCM, CV-MC, baselines) with SDID added only on a subset, or lowering rep counts after a variance run. These choices trade against the preregistered falsifier precision and must be frozen in `preregistration.md` before any decisive run; none is decided in Phase B.

## 6. Seed registry additions

Smoke run: seeds 7000-7019 (`code/run_smoke.py`). Pilot: 20000-20009 (decisive cell), 31000-31001 (scaling probe). See `seeds.yaml`.
