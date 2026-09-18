"""WP-D2 regression tests (preregistration_d2_addendum.md machinery)."""

import math

import numpy as np
import pytest

from applications import frontier_d2 as F

SEEDS = {"align": 60101, "donor": 60102, "time": 60103}


def _spiked_panel(n, T0, s, seed):
    rng = np.random.default_rng(seed)
    u = rng.normal(size=n)
    u /= np.linalg.norm(u)
    v = rng.normal(size=T0)
    return math.sqrt(s) * np.outer(u, v) \
        + rng.normal(size=(n, T0)), u * math.sqrt(s), v


def test_mp_median_inside_support_and_below_mean():
    for c in (0.25, 0.5, 1.0, 2.0, 4.0):
        sq = math.sqrt(c)
        assert (1 - sq) ** 2 <= F.mp_median(c) <= (1 + sq) ** 2
    assert F.mp_median(1.0) < 1.0


def test_pipeline_recovers_known_supercritical_spike():
    Y, a, v = _spiked_panel(40, 40, 9.0, seed=101)
    y1 = 3.0 * v
    res = F.analyze_unit(Y, y1, SEEDS, idx=0)
    assert res["k"] == 1
    assert 6.0 <= res["d"] <= 14.0
    assert res["p_align"] < F.ALIGN_P
    assert res["label"] == "RECOVERABLE"
    assert res["stability_pass"] is True


def test_visible_spike_without_treated_leverage_flags_misaligned():
    Y, a, v = _spiked_panel(40, 40, 9.0, seed=102)
    rng = np.random.default_rng(7)
    res = F.analyze_unit(Y, rng.normal(size=40), SEEDS, idx=1)
    assert res["k"] >= 1
    assert res["p_align"] >= F.ALIGN_P
    assert res["label"] == "FRAGILE-MISALIGNED"


def test_silent_panel_classifies_invisible():
    rng = np.random.default_rng(103)
    Y = rng.normal(size=(40, 40))
    res = F.analyze_unit(Y, rng.normal(size=40), SEEDS, idx=2)
    assert res["k"] == 0
    assert res["label"] == "FRAGILE-INVISIBLE"
    assert res["d"] == 0.0


def test_alignment_null_is_scale_free():
    rng = np.random.default_rng(104)
    basis = np.linalg.qr(rng.normal(size=(40, 40)))[0].T
    pool = rng.standard_normal((F.G_NULL, 40))
    e_hi, p_hi, r_hi = F.alignment_test(basis, basis[0] * 100.0, 1, pool)
    e_lo, p_lo, r_lo = F.alignment_test(basis, basis[0] * 0.01, 1, pool)
    assert e_hi == pytest.approx(e_lo, rel=1e-12)
    assert p_hi < 0.01 and r_hi > 5
    _, p_noise, _ = F.alignment_test(basis, rng.normal(size=40), 1, pool)
    assert p_noise > 0.05


def test_spearman_permutation_detects_negative_ordering():
    x = np.arange(20, dtype=float)
    y = -2.0 * x + 1.0
    rho, p = F.spearman_with_perm(x, y, seed=105)
    assert rho == pytest.approx(-1.0)
    assert p <= 1.0 / (F.N_PERM + 1.0)


def test_date_alarm_scan_runs_and_bounds_rate():
    rng = np.random.default_rng(106)
    Y = rng.normal(size=(30, 24))
    rows = F.date_alarm_scan(Y, [10, 14], seed=8880001, sigma2=1.0)
    assert len(rows) == 2
    assert all(0.0 <= r["p"] <= 1.0 for r in rows)


def test_nf_arms_logic():
    treated = {"label": "RECOVERABLE", "d": 9.0}
    out = F.nf_arms(treated, 2.0, ["p1"], np.array([1.0, 0.5]),
                    np.array([2.1, 1.9]), ["FRAGILE-INVISIBLE",
                                           "RECOVERABLE"])
    assert out["N3_certification_asymmetry"] is True
    assert out["NF"] is True
    assert out["P2_treated_separation"] is True


def test_load_basque_structure():
    from applications import loaders as L
    panel = L.load_basque()
    assert panel.Y.shape == (17, 43)
    assert L.BASQUE_TREATED in panel.units
    assert "Spain (Espana)" not in panel.units
    assert panel.treat_year == 1975
    assert bool(np.isfinite(panel.Y).all())


def test_load_bi63_wisconsin_structure():
    from applications import loaders as L
    panel = L.load_bi63_wisconsin(outcome="aadr")
    assert panel.Y.shape == (9, 19)
    assert panel.T0 == 16
    assert panel.n_donors == 8
    assert panel.units[-1] == "Wisconsin"
    assert panel.years[[0, -1]].tolist() == [1999, 2017]
    assert bool(np.isfinite(panel.Y).all())


def test_load_fbi_crime_proxy_filters_on_precoverage(tmp_path):
    from applications import loaders as L
    rows = []
    for state, coverage in (("KS", 99.0), ("CA", 98.0), ("NY", 80.0)):
        for year, month in ((2006, 10), (2006, 11), (2006, 12),
                            (2007, 1), (2007, 2), (2007, 3)):
            rows.append({
                "state": state, "year": year, "month": month,
                "actual": year + month, "rate": year / 100 + month,
                "population": 1_000_000, "coverage_pct": coverage,
            })
    path = tmp_path / "fbi.csv"
    import pandas as pd
    pd.DataFrame(rows).to_csv(path, index=False)
    panel = L.load_fbi_crime_proxy(
        str(path), treated_state="KS", treat_month="2007-01",
        states=["KS", "California", "NY"], min_pre_coverage=90.0)
    assert panel.units == ["CA", "KS"]
    assert panel.Y.shape == (2, 6)
    assert panel.T0 == 3
    assert panel.meta["coverage_excluded_states"] == ["NY"]
