# Preregistration C5 Addendum — repair-and-confirm package

**Date frozen:** 2026-08-25, BEFORE any C.5 replication was generated. **Parent:** preregistration.md (2026-08-24) and gate_g3_memo.md Section 5 (PIVOT). **Purpose:** address the four diagnosed failures with amended instruments, frozen here; outcomes may restore the GO criteria on the amended definitions. Nothing in the parent grids is re-run; all confirmation runs use FRESH seeds (15000-18999 ranges, registered in seeds.yaml).

---

## C5a. Kink-location confirmation (corrected primary estimator)

Defect being repaired: raw second differences on the nonuniform m-grid (spacing 0.2 vs 0.1) are mechanically biased toward wide-spaced points and dominated by descent-tail curvature.

Amended primary estimator (frozen): kink location := breakpoint argmin_b of the two-segment fit (constant plateau on m <= b; OLS line on m > b), b scanned over arange(0.60, 1.60, 0.05), applied to the mean-RMSE curve of GATED spectral_sc (full arm, r = 1). PASS iff |b_hat - 1| <= 0.15 in >= 4/5 c columns. The legacy curvature estimator is reported for continuity only.

Confirmation data (fresh seeds): c in {0.25, 0.5, 1, 2, 4} (equal-entries geometries as parent), m grid linspace(0.6, 1.6, 21), 500 reps, seeds 15000+i, methods {spectral_gated, donor_mean}. Existing merged data are additionally re-analyzed under the amended estimator (labeled exploratory, since estimators post-date those runs).

## C5b. Bite extension (theta = 3 arm)

Rationale: F predicts every method sits at sqrt(1 + sum theta_j) below the edge; the parent grid capped theta at 1.0 (plateau 1.41 sigma), so the >= 2 sigma collapse was unreachable BY CONSTRUCTION. At theta = 3 the predicted sub-edge plateau is exactly 2.0 sigma. This arm tests that prediction and the frontier's practical meaning: no weighting scheme escapes it.

Grid: theta = 3, alignment "first", r = 1, c in {0.5, 1, 2} (equal-entries geometries), m = parent 17-point grid, methods FULL6 (parent sweep set), reps 400, seeds 16000+i, diag light.

Revised bite criterion (PASS iff all three hold):
1. Each incumbent {scm_simplex, ridge_sc, mc_nn_cv} attains mean RMSE >= 1.95 sigma pooled over sub-edge cells (m <= 0.8) of every included c column;
2. Gated selector flag rate (k = 0) on the same cells >= 80%;
3. Flag rate on supercritical theta-3 cells (m >= 1.5) <= 20%.

## C5c. Diagnostics v2 under weak dependence and heteroskedasticity

Defects being repaired: Z_boot duplication inflation (size 0); Z_tw iid-calibration shattering under AR/het; silence gate false-firing ~100% under AR nulls. Prototype triage (40-80 rep pilots, exploratory, 2026-08-25): block-length tuning fails at rho = 0.7; circular-shift null reached 5.0% at rho = 0.7; permutation b = 20 conservative-to-liberal; sieve sketch invalid (implementation error, abandoned).

Frozen v2 instruments:
1. **z_shift (primary)**: p-value from B = 200 uniform random circular offsets of the time index; per offset, resid_statistic on the rotated split (basis = first T0 - Tp_eff columns, pseudo-post = remainder); p = (#\{z* >= z_obs\} + 1)/(B + 1). Validity claim restricted to H0 (position exchangeability); rejection rates under alternatives are reported as detection rates.
2. **Row standardization for het**: before any v2 statistic, divide each donor row by its basis-window robust scale sigma_i = median(|Y_i,t - Y_i,t-1|) / (0.6745 sqrt 2), computed on the basis window only; sigma_eff = median_i(sigma_i) feeds the TW centering. Applied always (approximately identity under homoskedasticity).
3. **gate_lrv**: silence gate v2: robust row-standardization (item 2 scales) followed by the MP edge test with sigma^2 replaced by the median per-row Newey-West long-run-variance estimate (Bartlett kernel, lag 4(T0/100)^(2/9)); k = 0 when lambda1 <= 1.05 (1+sqrt(c))^2 sigma2_lrv, else largest-gap rank among top-4. DESIGN NOTE (amended pre-run): the originally drafted shift-null gate was found degenerate before execution because YY' is column-permutation invariant, making every circular offset produce an identical scatter statistic; the LRV design targets the same failure mode (dependence/het inflating the effective bulk edge) through its parametric cause instead.
4. **z_perm20 (secondary)**: disjoint-block permutation, block = 20, otherwise as prototype.
References carried unchanged: z_tw (simulated iid null), classical trend t-test.

Battery: laws {gaussian, ar1 rho=0.3, ar1 rho=0.7, heteroskedastic ratio 4} x states {null, spiked-no-break m=1.6 theta=1, break delta=1, break delta=2}, geometry (160,160,80), reps 300, seeds 17000+i.

PASS rules: z_shift size in [3%, 8%] for EACH of ar03/ar07/het nulls (gaussian reported against nominal); gate_lrv false-fire <= 6% on every null state; z_shift detection >= 80% at break delta = 2; spiked-no-break rejection documented with 12% liberalism tolerance. z_perm20/gate_mp/z_tw/t-test rows reported for the record.

## C5d. Break-detection formalization (post-window statistic)

Ships the repaired statistic into the package: pre_trends_post_test = resid_statistic(full donor-PRE basis, REAL donor POST window), calibrated by the simulated iid finite-n null (G = 300, seed derivation as parent). Scope note: the parametric null is valid under Gaussian base noise (the C2(iv) DGP family); under non-Gaussian dependence the C5c shift machinery applies to the joint window only under H0.

Cells: structural_break delta in {0, 0.5, 1, 2}, m = 2.0, theta = 1, aligned, gaussian noise, geometry (160,160,80); plus delta = 2 orthogonal-theta-0 control; reps 400, seeds 18000+i.

PASS: power(delta = 2) >= 80% AND size(delta = 0) <= 15% (spiked-benign liberalism documented; projection-leakage inflation is the named mechanism). Full delta dose-response reported.

## Deviation log (amendments superseding parent items)

D1 replaces the parent primary kink estimator (defect diagnosed post hoc, gate_g3_memo Section 3.1).
D2 extends the bite region to theta = 3 (parent cap made the criterion unreachable; F predicted it ex ante).
D3 replaces Z_boot (degenerate) and recalibrates the gate under weak dependence; z_shift/gate_shift are new instruments frozen above.
D4 operationalizes break DETECTION through the donor-post window (parent C2(iv) instrumentation was pre-window-only and vacuous for detection).

All other parent rules (metrics, leakage, merge discipline, no-partial-figures) carry over unchanged.
