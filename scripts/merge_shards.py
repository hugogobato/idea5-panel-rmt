"""Merge and validate downloaded Colab shard outputs (plan Section 11).

Usage:  python3 scripts/merge_shards.py [--raw results_raw]

Validates EVERY manifest notebook before any figure is produced:
  file present, gzip intact, sha256 matches the self-reported meta JSON,
  row count exact, per-cell seed sets complete, no recorded cell errors.
On success writes merged parquet files per experiment family.
Exit 1 if anything fails; partial grids never reach figures.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from scm_frontier import ROW_COLS  # noqa: E402

FAMILY_OUT = {
    "c1": ROOT / "results_c1" / "risk_curves.parquet",
    "c2ii": ROOT / "results_c2" / "c2ii.parquet",
    "c2iii": ROOT / "results_c2" / "c2iii.parquet",
    "c2iv": ROOT / "results_c2" / "c2iv.parquet",
    "c5a": ROOT / "results_c1" / "c5a_kink_confirm.parquet",
    "c5b": ROOT / "results_c1" / "c5b_bite_extension.parquet",
}


def load_csv_gz(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rb") as fh:
        return pd.read_csv(io.BytesIO(fh.read()))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=str(ROOT / "results_raw"))
    args = ap.parse_args()
    raw = Path(args.raw)
    manifest = yaml.safe_load((ROOT / "colab" / "shard_manifest.yaml").read_text())

    failures = []
    frames: dict[str, list[pd.DataFrame]] = {}
    for entry in manifest["notebooks"]:
        name = entry["file"].replace(".ipynb", "")
        fam = entry["family"]
        if fam in ("onset", "c4", "c5c", "c5d"):
            continue  # special outputs validated separately below
        gz = raw / fam / f"{name}.csv.gz"
        meta_path = raw / fam / f"{name}_meta.json"
        if not gz.exists() or not meta_path.exists():
            failures.append(f"{name}: missing outputs")
            continue
        meta = json.loads(meta_path.read_text())
        # meta records the sha256 of the UNCOMPRESSED csv (notebook finalize)
        try:
            payload = gzip.decompress(gz.read_bytes())
        except Exception as exc:
            failures.append(f"{name}: gzip corrupt ({exc})")
            continue
        sha = hashlib.sha256(payload).hexdigest()
        if sha != meta.get("csv_sha256"):
            failures.append(f"{name}: sha256 mismatch")
            continue
        try:
            df = pd.read_csv(io.BytesIO(payload))
        except Exception as exc:
            failures.append(f"{name}: unreadable ({exc})")
            continue
        if list(df.columns) != list(ROW_COLS):
            failures.append(f"{name}: schema mismatch")
            continue
        expected_rows = sum(c["reps"] * (len(c["methods"])
                                          + (1 if c.get("diag", "none") != "none"
                                             else 0))
                            for c in entry["cells"] if isinstance(c, dict))
        if len(df) != expected_rows or int(meta.get("rows_expected", -1)) != expected_rows:
            failures.append(f"{name}: rows {len(df)} != {expected_rows}")
            continue
        if meta.get("cell_errors"):
            failures.append(f"{name}: cell errors {list(meta['cell_errors'])}")
            continue
        for cspec in entry["cells"]:
            if not isinstance(cspec, dict):
                continue
            sub = df[df["cell_id"] == cspec["cell_id"]]
            seeds = set(sub["rep_seed"].astype(int))
            want = set(range(manifest["seed_base"],
                             manifest["seed_base"] + cspec["reps"]))
            if seeds != want:
                failures.append(f"{name}/{cspec['cell_id']}: seed set incomplete")
                break
        else:
            frames.setdefault(fam, []).append(df)
            print(f"OK   {name}: {len(df)} rows")

    for out in FAMILY_OUT.values():
        out.parent.mkdir(exist_ok=True)
    for fam, dfs in frames.items():
        merged = pd.concat(dfs, ignore_index=True)
        merged.to_parquet(FAMILY_OUT[fam], index=False)
        print(f"merged {fam}: {len(merged)} rows -> {FAMILY_OUT[fam]}")

    # special outputs
    onset = raw / "onset" / "onset_verdict.json"
    if onset.exists():
        v = json.loads(onset.read_text())
        print("onset verdict:", json.dumps(v))
    else:
        print("note: onset verdict not yet downloaded")
    scaling = raw / "c4" / "scaling_probe.csv"
    print("scaling probe:", "present" if scaling.exists() else "not yet downloaded")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(" -", f)
        return 1
    print("\nALL SHARDS VALID. Grid complete; figures may be produced.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
