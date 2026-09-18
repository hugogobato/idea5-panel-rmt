"""Reproduce and close the disclosed CDC bi63 certification-only proxy."""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from applications import frontier_d2 as F  # noqa: E402
from applications import loaders as L  # noqa: E402

SEEDS = {"align": 71_001, "donor": 71_002, "time": 71_003, "perm": 71_004}
OUT = ROOT / "results_e" / "bi63"


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
        "k", "d", "p_align", "r_align", "label", "sigma2_1",
        "boot_donor", "boot_time", "diff_space", "sensitivity",
        "stability_agree", "stability_pass", "cv_rank_k", "ungated_k",
    )
    return {key: result[key] for key in keys}


def analyze_outcome(outcome: str) -> dict:
    panel = L.load_bi63_wisconsin(outcome=outcome)
    T0 = panel.T0
    treated_name = panel.units[panel.treated_idx]
    donor_indices = [i for i in range(len(panel.units)) if i != panel.treated_idx]
    Yd = panel.Y[donor_indices, :T0]
    y1 = panel.Y[panel.treated_idx, :T0]
    treated = _compact(F.analyze_unit(Yd, y1, SEEDS, idx=0))
    treated["unit"] = treated_name
    treated["rmse"] = _rmse(y1, Yd)
    print(
        f"[{outcome}] Wisconsin: k={treated['k']} d={treated['d']:.3f} "
        f"p={treated['p_align']:.4f} {treated['label']}", flush=True)

    placebos = []
    for idx, placebo_name in enumerate(
            sorted(u for u in panel.units if u != treated_name), start=1):
        placebo_idx = panel.units.index(placebo_name)
        keep = [i for i in range(len(panel.units)) if i != placebo_idx]
        Yp = panel.Y[keep, :T0]
        yp = panel.Y[placebo_idx, :T0]
        result = _compact(F.analyze_unit(Yp, yp, SEEDS, idx=idx))
        result["unit"] = placebo_name
        result["rmse"] = _rmse(yp, Yp)
        placebos.append(result)
        print(
            f"[{outcome}] {placebo_name}: k={result['k']} "
            f"d={result['d']:.3f} p={result['p_align']:.4f} "
            f"{result['label']}", flush=True)
    counts = {
        label: sum(row["label"] == label for row in placebos)
        for label in sorted({row["label"] for row in placebos})
    }
    return {
        "panel": panel.name,
        "outcome": outcome,
        "T0": T0,
        "n_d": panel.n_donors,
        "c": panel.n_donors / T0,
        "treated": treated,
        "placebos": placebos,
        "placebo_label_counts": counts,
        "all_placebos_recoverable": all(
            row["label"] == "RECOVERABLE" for row in placebos),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results = {outcome: analyze_outcome(outcome) for outcome in ("aadr", "deaths")}
    summary = {
        "study": "cdc_bi63_wisconsin_yearly_suicide_proxy",
        "preregistration": "docs/preregistrations/preregistration_bi63_wisconsin_addendum.md",
        "scope": "certification-only; no causal or Gate G4 claim",
        "outcomes": results,
        "closure": {
            "primary": "aadr",
            "certification_only": True,
            "applied_value": "INCREMENTAL-ONLY",
            "reason": (
                "annual all-suicide proxy does not reproduce restricted monthly "
                "handgun-suicide outcome and was inspected before freeze"
            ),
        },
    }
    path = OUT / "summary_bi63.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(_clean(summary), handle, indent=2)
        handle.write("\n")
    print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
