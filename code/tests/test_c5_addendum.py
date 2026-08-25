"""C5 addendum tests: circular-shift null, permutation null, gate_shift,
post-window break test, kink_breakpoint estimator (frozen C5 designs)."""

import numpy as np

import scm_frontier as sf


def _ar1(rng, n=159, T0=160, rho=0.7):
    e = rng.normal(size=(n, T0))
    e[:, 0] /= np.sqrt(1 - rho**2)
    for t in range(1, T0):
        e[:, t] = rho * e[:, t - 1] + e[:, t]
    return e


def test_kink_breakpoint_recovers_true_corner():
    ms = np.linspace(0.6, 1.6, 21)
    curve = np.where(ms <= 1.0, 1.414, 1.414 - 0.8 * (ms - 1.0))
    assert abs(sf.kink_breakpoint(ms, curve) - 1.0) < 0.051
    # spacing-bias regression test: nonuniform grid must not move the corner
    ms2 = np.sort(np.concatenate([np.linspace(0.2, 3.0, 15), [0.9, 1.1]]))
    curve2 = np.where(ms2 <= 1.0, 1.414, 1.414 - 0.25 * (ms2 - 1.0))
    assert abs(sf.kink_breakpoint(ms2, curve2) - 1.0) <= 0.05


def test_z_shift_size_under_ar07():
    ps = []
    for j in range(30):
        Y = _ar1(np.random.default_rng(500 + j))
        p, _, _ = sf.z_shift_pvalue(Y, 1.0, T_post=80, B=100)
        ps.append(p)
    rate = float(np.mean(np.array(ps) < 0.05))
    assert rate <= 0.15  # loose CI at 30 reps; nominal 0.05


def test_z_perm_no_inflation():
    # with-replacement defect regression: permutation must not pile p at ~1
    ps = []
    for j in range(20):
        Y = np.random.default_rng(700 + j).normal(size=(79, 120))
        p, _, _ = sf.z_perm_pvalue(Y, 1.0, T_post=40, B=100, block=20)
        ps.append(p)
    assert float(np.median(ps)) < 0.95


def test_gate_lrv_under_ar_and_spike():
    # AR(1) null: iid MP edge would false-fire; LRV-adjusted gate must not
    fires = 0
    for j in range(20):
        Y = _ar1(np.random.default_rng(900 + j), n=80, T0=160, rho=0.7)
        fires += sf.gate_lrv(Y) > 0
    assert fires <= 4
    # gaussian spiked panel must fire (spike inflates LRV only mildly)
    rng = np.random.default_rng(13)
    fired = 0
    for j in range(10):
        Y = rng.normal(size=(80, 160))
        a = rng.normal(size=80) * np.sqrt(6.0 / 79)
        Y = Y + np.outer(a, rng.normal(size=160))
        fired += sf.gate_lrv(Y) > 0
    assert fired >= 9


def test_pre_trends_post_test_detects_break():
    rng = np.random.default_rng(3)
    pre, post = rng.normal(size=(159, 160)), rng.normal(size=(159, 80))
    out0 = sf.pre_trends_post_test(pre, post, 1.0, G=100, seed=1)
    a = rng.normal(size=159) * np.sqrt(4.0 / 159)
    f_pre, f_post = rng.normal(size=160), rng.normal(size=80) + 3.0
    post_brk = np.outer(a, f_post) + rng.normal(size=(159, 80))
    pre_spk = np.outer(a, f_pre) + rng.normal(size=(159, 160))
    out1 = sf.pre_trends_post_test(pre_spk, post_brk, 1.0, G=100, seed=1)
    assert out1["z"] > out0["z"] + 5.0
    assert out1["p"] < 0.01


def test_robust_row_scales_het():
    rng = np.random.default_rng(5)
    Y = rng.normal(size=(40, 200)) * np.array(
        [4.0 if i % 2 else 0.25 for i in range(40)])[:, None]
    sc = sf.robust_row_scales(Y)
    ratio = sc[1::2].mean() / sc[::2].mean()
    assert 8.0 < ratio < 24.0  # true 16x; estimator is near-exact here
