"""Fetch a reproducible state-month panel from the official FBI CDE endpoint.

The current CDE host exposes the same JSON used by its web application without
an api.data.gov key.  One CSV and one metadata JSON are written.  The API's
``actuals`` key is preserved as the CSV column ``actual``; when
``--series-type estimates`` is used, it is the CDE estimated monthly count.
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
import time
import urllib.parse
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from applications.loaders import FBI_STATE_ABBR_TO_NAME  # noqa: E402

BASE = "https://cde.ucr.cjis.gov/LATEST/summarized/state"
DEFAULT_STATES = tuple(k for k in FBI_STATE_ABBR_TO_NAME if k != "DC")


def _fetch_json(url: str, attempts: int = 3) -> dict:
    error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "Idea5-panel-rmt/1.0"})
            with urllib.request.urlopen(request, timeout=45) as response:
                return json.load(response)
        except Exception as exc:  # pragma: no cover - live network path
            error = exc
            if attempt + 1 < attempts:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"FBI CDE request failed after {attempts} attempts: {url}") from error


def _state_rows(
    state: str, offense: str, from_month: str, to_month: str,
    series_type: str,
) -> tuple[list[dict], dict]:
    query = urllib.parse.urlencode(
        {"from": from_month, "to": to_month, "type": series_type})
    url = f"{BASE}/{state}/{urllib.parse.quote(offense)}?{query}"
    payload = _fetch_json(url)
    state_name = FBI_STATE_ABBR_TO_NAME[state]
    offense_key = f"{state_name} Offenses"
    offenses = payload.get("offenses", {})
    rates = offenses.get("rates", {}).get(offense_key)
    actuals = offenses.get("actuals", {}).get(offense_key)
    populations = payload.get("populations", {}).get("population", {}).get(state_name)
    participated = payload.get("populations", {}).get(
        "participated_population", {}).get(state_name, {})
    coverage = payload.get("tooltips", {}).get(
        "Percent of Population Coverage", {}).get(state_name, {})
    if not isinstance(rates, dict) or not isinstance(actuals, dict) \
            or not isinstance(populations, dict):
        raise ValueError(f"unexpected FBI CDE schema for {state}: {url}")
    if set(rates) != set(actuals) or set(rates) != set(populations):
        raise ValueError(f"FBI CDE month keys disagree for {state}")
    rows = []
    for date in sorted(rates, key=lambda x: (int(x[3:]), int(x[:2]))):
        month, year = (int(x) for x in date.split("-"))
        rows.append({
            "state": state,
            "state_name": state_name,
            "year": year,
            "month": month,
            "date": f"{year:04d}-{month:02d}",
            "offense": offense,
            "series_type": series_type,
            "actual": actuals[date],
            "rate": rates[date],
            "population": populations[date],
            "participated_population": participated.get(date),
            "coverage_pct": coverage.get(date),
        })
    return rows, {"url": url, "cde_properties": payload.get("cde_properties", {})}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--offense", default="violent-crime")
    parser.add_argument("--from-month", default="01-1999", metavar="MM-YYYY")
    parser.add_argument("--to-month", default="12-2017", metavar="MM-YYYY")
    parser.add_argument(
        "--series-type", choices=("estimates", "totals"), default="estimates")
    parser.add_argument(
        "--states", default=",".join(DEFAULT_STATES),
        help="comma-separated state abbreviations")
    parser.add_argument(
        "--output", type=pathlib.Path,
        default=ROOT / "data" / "raw" / "fbi_cde"
        / "fbi_cde_monthly_1999_2017.csv")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    states = [x.strip().upper() for x in args.states.split(",") if x.strip()]
    unknown = sorted(set(states) - set(FBI_STATE_ABBR_TO_NAME))
    if unknown:
        raise ValueError(f"unknown state abbreviations: {unknown}")
    if len(states) != len(set(states)):
        raise ValueError("duplicate states requested")
    rows, requests = [], []
    for i, state in enumerate(states, start=1):
        state_rows, request_meta = _state_rows(
            state, args.offense, args.from_month, args.to_month,
            args.series_type)
        rows.extend(state_rows)
        requests.append({"state": state, **request_meta})
        print(f"[{i:02d}/{len(states):02d}] {state}: {len(state_rows)} months", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "state", "state_name", "year", "month", "date", "offense",
        "series_type", "actual", "rate", "population",
        "participated_population", "coverage_pct",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    metadata_path = args.output.with_suffix(".meta.json")
    metadata = {
        "source": "FBI Crime Data Explorer",
        "base_endpoint": BASE,
        "offense": args.offense,
        "from": args.from_month,
        "to": args.to_month,
        "series_type": args.series_type,
        "states": states,
        "rows": len(rows),
        "requests": requests,
    }
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")
    print(f"Wrote {args.output} ({len(rows)} rows)")
    print(f"Wrote {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
