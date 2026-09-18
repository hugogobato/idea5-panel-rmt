"""WP-D1 trusted-benchmark reproduction tests (docs/applied/preprocessing_frozen.md).

Each test re-runs the frozen pipeline end-to-end (panels are tiny; total
runtime < 30 s) and asserts the published anchors within the tolerances
frozen in docs/applied/preprocessing_frozen.md Section 5. These are REGRESSION guards:
if any of them moves, the preprocessing or solver changed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from applications import loaders as L  # noqa: E402
from applications import synth_adh as S  # noqa: E402

V_SEED = 60_001  # config/seeds.yaml phase_d.wp_d1.v_multistart_seed


@pytest.fixture(scope="module")
def smoking():
    panel, cdf = L.load_smoking()
    return panel, cdf, L.smoking_adh_predictors(panel, cdf)


@pytest.fixture(scope="module")
def germany():
    panel = L.load_germany()
    raw = L.germany_raw()
    return panel, raw


def test_smoking_mirrors_agree(smoking):
    panel, cdf, _ = smoking
    assert panel.Y.shape == (39, 31)
    assert panel.units[-1] == "California"
    assert panel.T0 == 19
    assert panel.meta["mirrors_max_abs_diff"] <= 1e-6


def test_smoking_table1_treated_row_matches_publication(smoking):
    """The mirror must reproduce the paper's own treated-row means."""
    panel, cdf, _ = smoking
    X = L.smoking_adh_predictors(panel, cdf)
    ca = X[-1]
    assert abs(ca[0] - 10.08) < 0.01      # ln income
    assert abs(ca[2] - 89.42) < 0.01      # retail price
    assert abs(ca[3] - 24.28) < 0.01      # beer 84-88
    assert abs(ca[4] - 90.10) < 0.01      # cigsale 1988
    assert abs(ca[6] - 127.10) < 0.01     # cigsale 1975


def test_adh2010_smoking_weights_and_effects(smoking):
    panel, _, X = smoking
    Y, T0 = panel.Y, panel.T0
    fit = S.synth_adh_smoking(X, Y, T0, seed=V_SEED)
    w = fit["w"]
    donors = [u for u in panel.units if u != "California"]
    named = {u: float(w[i]) for i, u in enumerate(donors) if w[i] > 1e-6}
    target = {"Colorado": 0.164, "Connecticut": 0.069, "Montana": 0.199,
              "Nevada": 0.234, "Utah": 0.334}
    for u, t in target.items():
        assert abs(named.get(u, 0.0) - t) <= 0.01, (u, named.get(u, 0.0), t)
    assert len(named) == 5

    gap = Y[-1] - w @ Y[:-1]
    pre_mspe = float(np.mean(gap[:T0] ** 2))
    assert abs(pre_mspe - 3.0) <= 1.5
    assert abs(gap[T0:].mean() - (-19.6)) <= 1.5
    idx2000 = int(np.where(panel.years == 2000)[0][0])
    assert abs(gap[idx2000] - (-26.0)) <= 3.0
    ratio = float(np.mean(gap[T0:] ** 2) / pre_mspe)
    assert abs(ratio - 130.0) <= 78.0


def test_synthdid_trio_smoking(smoking):
    panel, _, _ = smoking
    Y = panel.Y
    sc = S.synthdid_sc(Y, panel.T0, panel.treated_idx)
    sdid = S.synthdid_att(Y, panel.T0, panel.treated_idx)
    did = S.synthdid_did(Y, panel.T0, panel.treated_idx)
    assert abs(sc["tau"] - (-19.6)) <= 0.4
    assert abs(sdid["tau"] - (-15.6)) <= 0.5
    assert abs(did - (-27.3)) <= 0.5


def test_adh2015_germany_reproduction(germany):
    panel, raw = germany
    Y, T0 = panel.Y, panel.T0
    donor_mask = np.arange(len(panel.units)) != panel.treated_idx
    X_train = L.germany_predictors(raw, "training")
    X_main = L.germany_predictors(raw, "main")
    Z0 = Y[donor_mask][:, T0 - 9:T0 + 1]
    Z1 = Y[panel.treated_idx, T0 - 9:T0 + 1]
    fit = S.synth_adh_germany(X_train, panel.treated_idx, Z0, Z1, X_main,
                              seed=V_SEED)
    w = fit["w"]
    donors = [u for u in panel.units if u != "West Germany"]
    named = {u: float(w[i]) for i, u in enumerate(donors) if w[i] > 1e-6}
    target = {"Austria": 0.42, "USA": 0.22, "Japan": 0.16,
              "Switzerland": 0.11, "Netherlands": 0.09}
    for u, t in target.items():
        assert abs(named.get(u, 0.0) - t) <= 0.02, (u, named.get(u, 0.0), t)
    assert len(named) == 5

    synth = w @ Y[donor_mask]
    tr = Y[panel.treated_idx]
    avg_gap = float((tr[T0:] - synth[T0:]).mean())
    assert abs(avg_gap - (-1600.0)) <= 250.0
    rel_2003 = float((synth[-1] - tr[-1]) / tr[-1])
    assert abs(rel_2003 - 0.12) <= 0.04


def test_germany_published_v_passthrough(germany):
    """Stage-2 at the paper's V must land on the same five donors."""
    panel, _ = germany
    raw = L.germany_raw()
    X_main = L.germany_predictors(raw, "main")
    v_pub = np.array([0.442, 0.134, 0.072, 0.001, 0.107, 0.245])
    st = S.fit_w_given_v(v_pub, *S.split_X(X_main, panel.treated_idx))
    donors = [u for u in panel.units if u != "West Germany"]
    nz = {u for u, x in zip(donors, st["w"]) if x > 1e-6}
    assert nz == {"Austria", "USA", "Japan", "Switzerland", "Netherlands"}
