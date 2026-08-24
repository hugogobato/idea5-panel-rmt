# Gate G3 Memo (Phase C close-out) — SKELETON, decision pending shard returns

**Date opened:** 2026-08-24. **Preregistration:** frozen same day (`preregistration.md`, deviations logged in its Section 9). **Execution vehicle:** 47 self-contained Colab notebooks (`colab/`, manifest `colab/shard_manifest.yaml`), repo `github.com/hugogobato/idea5-panel-rmt`. **Status:** notebooks built and validated; results pending upload/execution/aggregation. This memo acquires numbers ONLY through `scripts/merge_shards.py` + `scripts/make_figures.py`; nothing here may be filled by hand from partial grids.

## 1. Preregistered expectations (frozen before any run)

1. Kink: gated spectral SC mean-RMSE curvature peak within |m-1| <= 0.15 in >= 4 of 5 c columns (full arm, r = 1).
2. DE overlay: F tracks simulated gated curves inside MC bands across columns.
3. Onset slice: W1-statistic detectability onset at the largest size in [0.90, 1.10].
4. Bite: some incumbent >= 2.0 sigma in a c >= 1 region while gated flag power >= 80%, false alarms <= 20% on m >= 1.5 aligned cells.
5. Null battery: estimator ties; gate size <= 6%; Z_tw/Z_boot size <= 6% iid.
6. Calibration: Z_boot size in [3%, 8%] under AR(1)/het with power >= 80% spiked; classical t-test expected to break size under rho = 0.7.
7. Structural break: Z_boot power >= 80% at delta >= 1; control <= 8%; all estimators degrade vs control pairs.
8. Dense weak factors: no mechanism inversion (spectral must not lose to donor_mean); MC-vs-spectral winner reported honestly either way.

## 2. Results

PENDING merge. Slots: kink table; overlay verdicts per column; onset-by-size table; bite numbers; size/power tables; rank head-to-head accuracy by c; bootstrap multiplier as measured (G2 obligation 5); C2(ii) battery table; runtime/memory audit vs cost model.

## 3. Deviations observed during execution

PENDING (any notebook cell error, shard rerun, or aggregation anomaly gets logged here with cause).

## 4. Strongest baseline's best case / proposed method's failure regions

PENDING analysis (WP-C5 memo requirement).

## 5. Decision (exactly one)

PENDING. Rules per plan Section 7 Phase C give-up rules:
GO to Phase D iff preregistered bite + calibration (+ working fallback) hold; PIVOT if transition real but bite restricted; INCREMENTAL-ONLY if incumbents near-oracle everywhere or diagnostics never win; KILL if no sharp transition, methods dominated everywhere, estimand unidentifiable under maintained assumptions, or miscalibration unrepairable.
