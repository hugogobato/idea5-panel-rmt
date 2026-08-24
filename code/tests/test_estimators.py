"""WP-B2 unit tests: reference implementation correctness ladder.

Scientific pass rules (plan Section 7, WP-B2): shapes hold; r=0 sits at the
noise floor; infinite-spike panels are near-oracle; rank selection touches
pre-period data only; fixed seeds reproduce bitwise. Oracle beats every
method in its home turf; the trivial baseline never beats tuned methods by
more than Monte Carlo noise in favorable cells; the simplex solver converges
in >= 99% of favorable-cell runs.
"""

import inspect

import numpy as np
import pytest

import scm_frontier as sf

SIGMA = 1.0
METHODS = ("donor_mean", "scm_simplex", "ridge_sc", "spectral_sc", "mc_nn_cv", "sdid")


def run_cell(reps, **cfg):
    out = {m: [] for m in METHODS}
    infos = {m: [] for m in METHODS}
    oracle = []
    for j in range(reps):
        pan = sf.generate_panel(seed=10_000 + j, **cfg)
        Y, L = pan["Y"], pan["L"]
        T0, T_post = cfg.get("T0", 240), cfg.get("T_post", 100)
        donors, y1_pre = Y[1:], Y[0, :T0]
        y_star = Y[0, T0:]
        oracle.append(float(np.sqrt(np.mean((sf.oracle_predict(L[0, T0:]) - y_star) ** 2))) / SIGMA)
        for m in METHODS:
            info = {}
            pred = getattr(sf, m)(donors, y1_pre, info=info)
            out[m].append(float(np.sqrt(np.mean((pred - y_star) ** 2))) / SIGMA)
            infos[m].append(info)
    return {m: np.array(v) for m, v in out.items()}, infos, np.array(oracle)


def test_shapes():
    pan = sf.generate_panel(n=41, T0=80, T_post=30, seed=1)
    donors, y1_pre = pan["Y"][1:], pan["Y"][0, :80]
    for m in METHODS:
        pred = getattr(sf, m)(donors, y1_pre)
        assert pred.shape == (30,)
        assert np.all(np.isfinite(pred))


def test_determinism():
    p1 = sf.generate_panel(n=41, T0=80, T_post=30, seed=7)
    p2 = sf.generate_panel(n=41, T0=80, T_post=30, seed=7)
    assert np.array_equal(p1["Y"], p2["Y"])
    donors, y1_pre = p1["Y"][1:], p1["Y"][0, :80]
    for m in METHODS:
        a = getattr(sf, m)(donors, y1_pre)
        b = getattr(sf, m)(p2["Y"][1:], p2["Y"][0, :80])
        assert np.array_equal(a, b), m


def test_leakage_interface():
    for m in METHODS:
        params = inspect.signature(getattr(sf, m)).parameters
        assert not any("post" in p and "treated" in p for p in params), m
        assert not any(p in ("y_full", "y_treated") for p in params), m


def test_null_floor():
    res, infos, orc = run_cell(
        reps=15, n=81, T0=160, T_post=60, r=1,
        spike_strengths=(0.2,), treated_share=(0.05,), alignment="none",
    )
    assert abs(orc.mean() - 1.0) < 0.07
    for m in METHODS:
        excess = (res[m] - orc).mean()
        assert -0.02 <= excess <= 0.25, (m, excess)


def test_infinite_spike_near_oracle():
    res, _, orc = run_cell(
        reps=15, n=81, T0=160, T_post=60, r=1,
        spike_strengths=(400.0,), treated_share=(0.5,), alignment="first",
    )
    for m in METHODS:
        excess = (res[m] - orc).mean()
        bound = 0.10 if m != "donor_mean" else 0.45
        assert excess <= bound, (m, excess)


def test_favorable_cell_ordering():
    """Oracle floor holds; trivial baseline cannot beat tuned methods; spectral
    beats donor-mean with paired t >> 3 in the mechanism-favorable cell."""
    reps = 20
    cfg = dict(n=81, T0=160, T_post=60, r=1, spike_strengths=(6.0,),
               treated_share=(1.0,), alignment="first")
    res, _, orc = run_cell(reps=reps, **cfg)
    for m in METHODS:
        excess = (res[m] - orc).mean()
        assert excess > -0.03, (m, "oracle violated beyond MC noise")
        assert excess < 0.45, (m, "unstable in home turf")
    d = res["donor_mean"] - res["spectral_sc"]
    assert d.mean() / (d.std() / np.sqrt(reps)) > 3.0
    for m in ("ridge_sc", "spectral_sc", "mc_nn_cv"):
        assert (res[m] < res["donor_mean"]).mean() >= 0.8, (
            m, "trivial baseline usually beats tuned method"
        )


def test_simplex_solver_stability():
    _, infos, _ = run_cell(
        reps=30, n=81, T0=160, T_post=60, r=1,
        spike_strengths=(6.0,), treated_share=(1.0,), alignment="first",
    )
    ok = [i["scm_success"] for i in infos["scm_simplex"]]
    assert np.mean(ok) >= 0.99, ok


def test_rank_selector():
    null = sf.generate_panel(n=81, T0=160, T_post=10, r=1,
                             spike_strengths=(0.2,), treated_share=(0.0,),
                             alignment="none", seed=3)
    c = 81 / 160
    assert sf.rank_selector(null["Y"][1:, :160], SIGMA, c) == 0
    spik = sf.generate_panel(n=81, T0=160, T_post=10, r=1,
                             spike_strengths=(6.0,), treated_share=(1.0,),
                             alignment="first", seed=4)
    assert sf.rank_selector(spik["Y"][1:, :160], SIGMA, c) >= 1


def test_invert_bbp_roundtrip():
    c = 0.5
    s_vals = np.linspace(1.5 * np.sqrt(c), 20.0, 50)
    lam = 1.0 + s_vals + c + c / s_vals
    s_back = sf.invert_bbp(lam, c)
    assert np.max(np.abs(s_back - s_vals)) < 1e-9
    below = sf.invert_bbp(np.array([(1 + np.sqrt(c)) ** 2 * 0.999]), c)
    assert np.all(np.isnan(below))


def test_tw_statistic_null_scale():
    vals = []
    for j in range(10):
        pan = sf.generate_panel(n=81, T0=160, T_post=5, r=1,
                                spike_strengths=(0.2,), treated_share=(0.0,),
                                alignment="none", seed=200 + j)
        vals.append(sf.tw_statistic(pan["Y"][1:, :160], SIGMA))
    vals = np.array(vals)
    assert np.all(np.isfinite(vals)) and abs(vals.mean()) < 3.0
