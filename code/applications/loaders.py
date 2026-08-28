"""Data ingestion for WP-D1 (preprocessing_frozen.md Section 3).

Raw extracts live under data/raw/ and are immutable; every file carries a
sha256 in data/raw/SHA256SUMS. These loaders implement the frozen cleaning
contract EXACTLY; any change requires a new freeze document first.

Panels
------
smoking : California Proposition 99 panel, 39 states x 31 years (1970-2000),
    treated unit California from 1989. Outcome packs per capita; ADH 2010
    covariates (log personal income, beer, retail price, % aged 15-24).
germany : West Germany reunification panel, 17 countries x 44 years
    (1960-2003), treated unit West Germany from 1990. Outcome GDP per capita
    (PPP, current USD) plus the six ADH 2015 predictors.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

RAW = None  # resolved lazily so the module works from any cwd


def _raw_dir() -> str:
    global RAW
    if RAW is None:
        import pathlib

        RAW = pathlib.Path(__file__).resolve().parents[2] / "data" / "raw"
    return str(RAW)


SMOKING_DONOR_ORDER = (
    "Alabama", "Arkansas", "Colorado", "Connecticut", "Delaware",
    "Georgia", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas",
    "Kentucky", "Louisiana", "Maine", "Minnesota", "Mississippi",
    "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire",
    "New Mexico", "North Carolina", "North Dakota", "Ohio", "Oklahoma",
    "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota",
    "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "West Virginia",
    "Wisconsin", "Wyoming", "California",
)
SMOKING_TREATED = "California"
SMOKING_TREAT_YEAR = 1989
GERMANY_TREATED = "West Germany"
GERMANY_TREAT_YEAR = 1990


@dataclass
class Panel:
    """Wide panel bundle. Rows follow donors-then-treated order."""

    name: str
    units: list[str]
    years: np.ndarray
    Y: np.ndarray            # n x T outcomes
    treated_idx: int
    treat_year: int
    meta: dict

    @property
    def T0(self) -> int:
        return int(np.searchsorted(self.years, self.treat_year))

    @property
    def n_donors(self) -> int:
        return self.Y.shape[0] - 1


def load_smoking(outcome_only_csv: str | None = None,
                 covariates_csv: str | None = None) -> tuple[Panel, pd.DataFrame]:
    """Load the Prop 99 panel; validate the two independent raw mirrors
    against each other (synthdid outcome extract vs covariate mirror)."""
    import os

    raw = _raw_dir()
    outcome_only_csv = outcome_only_csv or os.path.join(
        raw, "california_prop99_synthdid.csv")
    covariates_csv = covariates_csv or os.path.join(
        raw, "smoking_adh_covariates_facure.csv")

    sdf = pd.read_csv(outcome_only_csv, sep=";")
    expected_cols = {"State", "Year", "PacksPerCapita", "treated"}
    if not expected_cols.issubset(sdf.columns):
        raise ValueError(f"smoking outcome csv missing columns {expected_cols}")
    cdf = pd.read_csv(covariates_csv)
    expected_cov = {
        "state", "year", "cigsale", "lnincome", "beer", "age15to24",
        "retprice", "california", "after_treatment",
    }
    if not expected_cov.issubset(cdf.columns):
        raise ValueError(f"smoking covariate csv missing columns {expected_cov}")

    # Mirror cross-check: identical outcome values on the shared key.
    cdf["state_name"] = cdf["state"].map(_smoking_state_code_to_name(cdf))
    merged = sdf.merge(
        cdf[["state_name", "year", "cigsale"]],
        left_on=["State", "Year"], right_on=["state_name", "year"],
        how="outer", indicator=True,
    )
    if not (merged["_merge"] == "both").all():
        raise ValueError("smoking mirrors disagree on the state-year grid")
    max_abs_diff = float(
        np.max(np.abs(merged["PacksPerCapita"] - merged["cigsale"]))
    )
    if max_abs_diff > 1e-6:
        raise ValueError(
            f"smoking mirrors disagree on outcome values (max |diff| "
            f"{max_abs_diff})")

    piv = sdf.pivot(index="State", columns="Year",
                    values="PacksPerCapita").loc[list(SMOKING_DONOR_ORDER)]
    years = piv.columns.to_numpy(dtype=int)
    Y = piv.to_numpy(dtype=float)
    if list(piv.index) != list(SMOKING_DONOR_ORDER):
        raise ValueError("smoking unit order mismatch")
    if years[0] != 1970 or years[-1] != 2000 or Y.shape != (39, 31):
        raise ValueError("smoking grid mismatch")

    panel = Panel(
        name="smoking_prop99", units=list(SMOKING_DONOR_ORDER), years=years,
        Y=Y, treated_idx=len(SMOKING_DONOR_ORDER) - 1,
        treat_year=SMOKING_TREAT_YEAR,
        meta={"outcome": "packs per capita",
              "mirrors_max_abs_diff": max_abs_diff},
    )
    cdf_out = cdf.copy()
    cdf_out["State"] = cdf_out["state_name"]
    return panel, cdf_out


def _smoking_state_code_to_name(cdf: pd.DataFrame) -> dict[int, str]:
    codes = sorted(cdf["state"].unique())
    names = SMOKING_DONOR_ORDER
    if len(codes) != len(names):
        raise ValueError("unexpected number of state codes in covariate csv")
    # The mirror orders states alphabetically by FIPS-assigned listing; the
    # mapping below is verified against the treated flag instead of assuming.
    calif_codes = set(cdf.loc[cdf["california"] == True, "state"])  # noqa: E712
    if len(calif_codes) != 1:
        raise ValueError("could not identify California code uniquely")
    calif_code = calif_codes.pop()
    remaining = [c for c in codes if c != calif_code]
    others_alphabetical = sorted(n for n in names if n != SMOKING_TREATED)
    mapping = dict(zip(remaining, others_alphabetical))
    mapping[calif_code] = SMOKING_TREATED
    return mapping


def smoking_adh_predictors(panel: Panel, cdf: pd.DataFrame) -> np.ndarray:
    """ADH 2010 X matrix, rows in paper order (Table 1):

    ln(income pc, 1997$) mean 1980-88; % aged 15-24 mean 1980-88; retail
    price mean 1980-88; beer pc mean 1984-88; cigsale levels 1988, 1980,
    1975. Row m stacks [donors..., treated].
    """
    units = panel.units

    def series(var: str) -> pd.DataFrame:
        sub = cdf.pivot(index="State", columns="year", values=var)
        return sub.reindex(units)

    def win_mean(var: str, y0: int, y1: int) -> np.ndarray:
        s = series(var)
        return s.loc[:, y0:y1].mean(axis=1).to_numpy()

    cols = [
        win_mean("lnincome", 1980, 1988),
        win_mean("age15to24", 1980, 1988) * 100.0,
        win_mean("retprice", 1980, 1988),
        win_mean("beer", 1984, 1988),
        series("cigsale").loc[:, 1988].to_numpy(),
        series("cigsale").loc[:, 1980].to_numpy(),
        series("cigsale").loc[:, 1975].to_numpy(),
    ]
    X = np.column_stack(cols)
    if not np.isfinite(X).all():
        raise ValueError("non-finite ADH predictor entries")
    return X


def load_germany(tab_path: str | None = None) -> Panel:
    """Load the ADH 2015 (updated erratum archive) reunification panel."""
    import os

    raw = _raw_dir()
    tab_path = tab_path or os.path.join(raw, "repgermany_updated.tab")
    df = pd.read_csv(tab_path, sep="\t")
    df["country"] = df["country"].str.strip('"')
    if df.shape != (748, 13):
        raise ValueError(f"germany grid mismatch: {df.shape}")
    order = df.drop_duplicates("index").sort_values("index")[
        "country"].tolist()
    if GERMANY_TREATED not in order:
        raise ValueError("West Germany missing")
    piv = df.pivot(index="country", columns="year", values="gdp").reindex(order)
    Y = piv.to_numpy(float)
    if not np.isfinite(Y).all():
        raise ValueError("non-finite GDP entries")
    years = piv.columns.to_numpy(int)
    return Panel(
        name="german_reunification", units=order, years=years, Y=Y,
        treated_idx=order.index(GERMANY_TREATED), treat_year=GERMANY_TREAT_YEAR,
        meta={"outcome": "gdp per capita PPP current USD"},
    )


def germany_predictors(df_raw: pd.DataFrame, stage: str) -> np.ndarray:
    """ADH 2015 X matrices following rep_updated.r verbatim.

    stage='training': predictors gdp/trade/infrate mean 1971-80; industry
    mean 1971-80; schooling mean{1970,1975}; invest70 at 1980.
    stage='main': gdp/trade/infrate mean 1981-90; industry mean 1981-90;
    schooling mean{1980,1985}; invest80 at 1980.
    """
    d = df_raw.copy()
    d["country"] = d["country"].str.strip('"')
    d = d.set_index(["country", "year"])

    def cell(var: str, c: str, y: int) -> float:
        return float(d.loc[(c, y), var])

    def win_mean(var: str, c: str, y0: int, y1: int) -> float:
        vals = [cell(var, c, y) for y in range(y0, y1 + 1)]
        arr = np.array(vals, float)
        if not np.isfinite(arr).any():
            return np.nan
        return float(np.nanmean(arr))

    def sp_mean(var: str, c: str, ys: tuple[int, ...]) -> float:
        return float(np.nanmean([cell(var, c, y) for y in ys]))

    countries = d.index.get_level_values("country").unique().tolist()
    if stage == "training":
        rows = []
        for c in countries:
            rows.append([
                win_mean("gdp", c, 1971, 1980),
                win_mean("trade", c, 1971, 1980),
                win_mean("infrate", c, 1971, 1980),
                win_mean("industry", c, 1971, 1980),
                sp_mean("schooling", c, (1970, 1975)),
                cell("invest70", c, 1980),
            ])
    elif stage == "main":
        rows = []
        for c in countries:
            rows.append([
                win_mean("gdp", c, 1981, 1990),
                win_mean("trade", c, 1981, 1990),
                win_mean("infrate", c, 1981, 1990),
                win_mean("industry", c, 1981, 1990),
                sp_mean("schooling", c, (1980, 1985)),
                cell("invest80", c, 1980),
            ])
    else:
        raise ValueError(stage)
    return np.array(rows, float)


def germany_raw() -> pd.DataFrame:
    import os

    df = pd.read_csv(os.path.join(_raw_dir(), "repgermany_updated.tab"),
                     sep="\t")
    df["country"] = df["country"].str.strip('"')
    return df


BONANDER_TREATED = "Florida"
BONANDER_TREAT_TIME = 82          # Oct 2005, 1-indexed months since 1999-01
BONANDER_T0 = 81
BONANDER_DONORS = (
    "Arkansas", "Connecticut", "Delaware", "Hawaii", "Iowa", "Maine",
    "Maryland", "Massachusetts", "Nebraska", "New Jersey", "New York",
    "North Dakota", "Ohio", "Rhode Island", "Wyoming",
)
BONANDER_OUTCOMES = (
    "HomicideRates", "p100khomicide_firearm", "p100khomicide_exclfirearm",
    "Firearm.Suicide.Rates", "Homicide_count",
)

BASQUE_TREATED = "Basque Country (Pais Vasco)"
BASQUE_TREAT_YEAR = 1975
BASQUE_NATIONAL_AGGREGATE = "Spain (Espana)"


def load_bonander(csv_path: str | None = None, outcome: str = "HomicideRates") -> Panel:
    """Load the Bonander et al. (2021 AJE) Florida SYG panel (Phase E).

    Source: OSF `rvayc-osfstorage-archive.zip` → `syg_data.csv`
    `sha256 7f8bd93b6add9e5e1b29d3ae738d30e90c17a2e1c867ebdee3c411f0341095b9`
    pinned in `data/raw/bonander/syg_data.csv` (E2, CC-BY). The file is
    16 states × 192 months (1999-01:2014-12, `time 1:192`). The panel
    returned is the *primary* outcome `HomicideRates` unless another
    registered outcome is requested; all registered outcomes share the
    same unit/time grid.
    """
    import os

    if outcome not in BONANDER_OUTCOMES:
        raise ValueError(f"unknown bonander outcome {outcome!r}")
    raw = _raw_dir()
    csv_path = csv_path or os.path.join(raw, "bonander", "syg_data.csv")
    df = pd.read_csv(csv_path)
    if df.shape[0] != 3072:
        raise ValueError(f"bonander row mismatch: {df.shape}")
    for col in ("State", "time", outcome):
        if col not in df.columns:
            raise ValueError(f"bonander missing column {col}")
    # 16-state analytic sample (paper) — order alphabetically, Florida treated
    states = sorted(df["State"].unique().tolist())
    if len(states) != 16 or BONANDER_TREATED not in states:
        raise ValueError(f"bonander state set mismatch: {states}")
    # time 1..192, pre 1:81 (1999-01:2005-09), post 82:192
    piv = df.pivot(index="State", columns="time", values=outcome).reindex(states)
    if piv.shape != (16, 192):
        raise ValueError(f"bonander pivot mismatch: {piv.shape}")
    if piv.isna().any().any():
        raise ValueError("bonander NA in pivot")
    Y = piv.to_numpy(float)
    years = piv.columns.to_numpy(int)  # time, not calendar year
    return Panel(
        name=f"bonander_{outcome}", units=states, years=years, Y=Y,
        treated_idx=states.index(BONANDER_TREATED),
        treat_year=int(BONANDER_TREAT_TIME),  # time index, not year
        meta={"outcome": outcome, "T0": BONANDER_T0, "treat_time": BONANDER_TREAT_TIME,
              "donors": list(BONANDER_DONORS), "provenance": "OSF rvayc E2 CC-BY"},
    )


def load_basque(csv_path: str | None = None) -> Panel:
    """Load the Abadie-Gardeazabal (2003) Basque panel (Wave-2, D2b).

    Source: R `Synth` package dataset (CRAN mirror), sha256-pinned as
    data/raw/basque_synth_ag2003.csv; provenance level E2 (community
    mirror of the authors' released extract). The national aggregate
    row ("Spain (Espana)") is excluded from the donor pool exactly as in
    the published analyses. Outcome: real GDP per capita (thousands of
    1986 USD? units as distributed; used only through centered spectra).
    """
    import os

    raw = _raw_dir()
    csv_path = csv_path or os.path.join(raw, "basque_synth_ag2003.csv")
    df = pd.read_csv(csv_path)
    if df.shape != (774, 17):
        raise ValueError(f"basque grid mismatch: {df.shape}")
    if not np.isfinite(df["gdpcap"]).all():
        raise ValueError("non-finite gdpcap entries")
    order = (df.drop_duplicates("regionno").sort_values("regionno"))
    names = order["regionname"].tolist()
    if BASQUE_TREATED not in names or BASQUE_NATIONAL_AGGREGATE not in names:
        raise ValueError("basque region set mismatch")
    keep = [n for n in names if n != BASQUE_NATIONAL_AGGREGATE]
    piv = df.pivot(index="regionname", columns="year",
                   values="gdpcap").reindex(keep)
    Y = piv.to_numpy(float)
    years = piv.columns.to_numpy(int)
    if Y.shape != (17, 43) or years[0] != 1955 or years[-1] != 1997:
        raise ValueError("basque grid shape/year mismatch")
    return Panel(
        name="basque_ag2003", units=keep, years=years, Y=Y,
        treated_idx=keep.index(BASQUE_TREATED),
        treat_year=BASQUE_TREAT_YEAR,
        meta={"outcome": "real GDP per capita (Synth AG2003 extract)",
              "provenance": "E2 mirror, CRAN Synth data/basque.rda"},
    )
