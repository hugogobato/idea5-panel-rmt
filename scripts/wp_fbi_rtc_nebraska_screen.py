"""Preregistered Nebraska RTC monthly violent-crime frontier screen."""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from applications import frontier_d2 as F  # noqa: E402
from applications import loaders as L  # noqa: E402

SEEDS = {"align": 72_001, "donor": 72_002, "time": 72_003, "perm": 72_004}
Z_SHIFT_SEED = 8_880_001
POSTTEST_SEED = 7_770_001
CUTOFF_IDX = [36, 60]
STATES = ("CA", "CT", "DE", "HI", "MD", "MA", "NJ", "NY", "RI", "NE")
OUT = ROOT / "results_e" / "fbi_rtc_nebraska"


def _clean(obj):
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


def _rmse(y1: np.ndarray, Yd: np.ndarray) -> float:
    y1c = F.center_rows(y1[None, :])[0]
    Ydc = F.center_rows(Yd)
    weights, *_ = np.linalg.lstsq(Ydc.T, y1c, rcond=None)
    return float(np.sqrt(np.mean((y1c - weights @ Ydc) ** 2)))


def _compact(result: dict) -> dict:
    keys = (
        "k", "d", "p_align", "r_align", "label", "sigma2_1", "edge",
        "boot_donor", "boot_time", "diff_space", "sensitivity",
        "stability_agree", "stability_pass", "cv_rank_k", "ungated_k",
        "evals_top12",
    )
    return {key: result[key] for key in keys}


def _panel(outcome: str) -> L.Panel:
    return L.load_fbi_crime_proxy(
        treated_state="NE", treat_month=200701, outcome=outcome,
        states=STATES, min_pre_coverage=90.0,
        start_month=199901, end_month=201712)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    panel = _panel("rate")
    T0 = panel.T0
    treated = panel.treated_idx
    donor_idx = [i for i in range(len(panel.units)) if i != treated]
    Yd_pre = panel.Y[donor_idx, :T0]
    y1_pre = panel.Y[treated, :T0]
    Yd_post = panel.Y[donor_idx, T0:]
    primary = _compact(F.analyze_unit(Yd_pre, y1_pre, SEEDS, idx=0))
    primary["unit"] = "NE"
    primary["rmse"] = _rmse(y1_pre, Yd_pre)
    primary["date_alarms"] = F.date_alarm_scan(
        Yd_pre, CUTOFF_IDX, Z_SHIFT_SEED, primary["sigma2_1"])
    primary["alarm_rate"] = float(np.mean(
        [row["alarm"] for row in primary["date_alarms"]]))
    primary["P3_pass"] = bool(primary["alarm_rate"] <= F.ALARM_MAX)
    primary["posttest_control"] = F.posttest_control(
        Yd_pre, Yd_post, primary["sigma2_1"], POSTTEST_SEED)
    print(
        f"Nebraska rate: k={primary['k']} d={primary['d']:.3f} "
        f"p={primary['p_align']:.4f} {primary['label']}", flush=True)

    placebos = []
    for index, unit in enumerate(sorted(u for u in panel.units if u != "NE"), start=1):
        placebo_idx = panel.units.index(unit)
        keep = [i for i in range(len(panel.units)) if i != placebo_idx]
        Yp = panel.Y[keep, :T0]
        yp = panel.Y[placebo_idx, :T0]
        result = _compact(F.analyze_unit(Yp, yp, SEEDS, idx=index))
        result["unit"] = unit
        result["rmse"] = _rmse(yp, Yp)
        placebos.append(result)
        print(
            f"{unit}: k={result['k']} d={result['d']:.3f} "
            f"p={result['p_align']:.4f} {result['label']}", flush=True)

    placebo_d = np.array([row["d"] for row in placebos])
    placebo_rmse = np.array([row["rmse"] for row in placebos])
    labels = [row["label"] for row in placebos]
    q75 = float(np.percentile(placebo_rmse, 75))
    primary["poor_looking"] = bool(primary["rmse"] > q75)
    nf = F.nf_arms(
        primary, primary["rmse"], [row["unit"] for row in placebos],
        placebo_d, placebo_rmse, labels)
    rho, perm_p = F.spearman_with_perm(placebo_d, placebo_rmse, SEEDS["perm"])
    nf["P1_descriptive"] = {"rho": rho, "perm_p": perm_p, "gated": False}
    nf["P3_date_alarms_pass"] = primary["P3_pass"]
    nf["stability_pass"] = primary["stability_pass"]
    qualified = bool(
        nf["NF"] and nf["P2_treated_separation"]
        and primary["P3_pass"] and primary["stability_pass"])
    if qualified:
        disposition = "QUALIFIED-FOR-TRUSTED-REPLICATION"
    elif primary["label"].startswith("FRAGILE"):
        disposition = "FRAGILE-BUT-SCREEN-FAILED"
    else:
        disposition = "DISCARD-NO-FRAGILE-TREATED"

    counts_panel = _panel("actual")
    Yd_counts = counts_panel.Y[donor_idx, :T0]
    y1_counts = counts_panel.Y[treated, :T0]
    counts = _compact(F.analyze_unit(Yd_counts, y1_counts, SEEDS, idx=100))

    summary = {
        "study": "fbi_cde_nebraska_rtc_monthly_violent_crime",
        "preregistration": "docs/preregistrations/preregistration_fbi_rtc_nebraska_addendum.md",
        "panel": {
            "units": panel.units, "T": int(panel.Y.shape[1]), "T0": T0,
            "n_d": panel.n_donors, "c": panel.n_donors / T0,
            "treat_month": 200701, "coverage_floor": 90.0,
            "coverage_excluded_states": panel.meta["coverage_excluded_states"],
        },
        "primary_rate": primary,
        "placebos": placebos,
        "placebo_label_counts": {
            label: sum(row["label"] == label for row in placebos)
            for label in sorted(set(labels))
        },
        "nf": nf,
        "registered_count_scale_audit": counts,
        "qualified_for_trusted_replication": qualified,
        "disposition": disposition,
    }
    path = OUT / "summary_fbi_rtc_nebraska.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(_clean(summary), handle, indent=2)
        handle.write("\n")
    print(f"{disposition}: wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
