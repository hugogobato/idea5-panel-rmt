"""Phase C extension tests: CV-rank comparator, pre-trends statistics,
classical trend test, spectral SC full output (gate + CI), and the shared
cell runner contract (preregistration Sections 2, 5.2, 8.8)."""

import numpy as np

import scm_frontier as sf


def _spiked_panel(seed=11, s=6.0, theta=1.0, n=81, T0=160, T_post=60):
    return sf.generate_panel(n=n, T0=T0, T_post=T_post, r=1,
                             spike_strengths=(s,), treated_share=(theta,),
                             alignment="first", seed=seed)


def test_cv_rank_selector_recovers_rank():
    pan = _spiked_panel()
    k = sf.cv_rank_selector(pan["Y"][1:, :160], pan["Y"][0, :160])
    assert k == 1
    null = sf.generate_panel(n=81, T0=160, T_post=10, r=1,
                             spike_strengths=(0.2,), treated_share=(0.0,),
                             alignment="none", seed=5)
    assert sf.cv_rank_selector(null["Y"][1:, :160], null["Y"][0, :160]) == 0


def test_classical_trend_ttest_null_calibration():
    rng = np.random.default_rng(99)
    ps = []
    for _ in range(60):
        y = rng.normal(0.0, 1.0, size=160)
        ps.append(sf.classical_trend_ttest(y))
    ps = np.array(ps)
    assert np.mean(ps < 0.05) <= 0.15
    assert 0.2 < ps.mean() < 0.8


def test_resid_statistic_detects_break():
    sigma, c = 1.0, 81 / 160
    n_d = 80
    rng = np.random.default_rng(7)
    # H0: independent windows
    z0, _ = sf.resid_statistic(rng.normal(size=(n_d, 100)),
                               rng.normal(size=(n_d, 50)), sigma, c)
    # break: rank-one mean shift in post
    a = rng.normal(size=n_d)
    post = rng.normal(size=(n_d, 50)) + 3.0 * np.outer(a, np.ones(50))
    z1, _ = sf.resid_statistic(rng.normal(size=(n_d, 100)), post, sigma, c)
    assert z1 > z0 + 5.0


def test_z_boot_size_and_power():
    # size on iid Gaussian null panels at nominal 5%
    rejections = 0
    reps = 30
    for j in range(reps):
        Y = np.random.default_rng(300 + j).normal(size=(80, 160))
        p, _, _ = sf.z_boot_pvalue(Y, 1.0, T_post=40, B=100)
        rejections += p < 0.05
    assert rejections <= 6  # <= 20% of 30 reps; nominal ~1.5 expected
    # power under an in-window factor-law break
    rng = np.random.default_rng(4242)
    n_rej = 0
    for j in range(12):
        Yp = rng.normal(size=(80, 120))
        a = rng.normal(size=80)
        Yb = np.hstack([Yp, rng.normal(size=(80, 40))
                        + 4.0 * np.outer(a, np.ones(40))])
        p, _, _ = sf.z_boot_pvalue(Yb, 1.0, T_post=40, B=100)
        n_rej += p < 0.05
    assert n_rej >= 9


def test_spectral_sc_full_gate_and_ci():
    pan_null = sf.generate_panel(n=81, T0=160, T_post=60, r=1,
                                 spike_strengths=(0.2,), treated_share=(0.05,),
                                 alignment="none", seed=21)
    out_n = sf.spectral_sc_full(pan_null["Y"][1:], pan_null["Y"][0, :160],
                                sigma=1.0, c=81 / 160)
    assert out_n["k_gated"] == 0
    assert np.allclose(out_n["pred_gated"], out_n["pred_gated"][0])
    assert (out_n["ci_lo"] < out_n["ci_hi"]).all()

    pan_sp = _spiked_panel(seed=22)
    out_s = sf.spectral_sc_full(pan_sp["Y"][1:], pan_sp["Y"][0, :160],
                                sigma=1.0, c=81 / 160)
    assert out_s["k_gated"] >= 1 and out_s["k_ungated"] == out_s["k_gated"]
    y_star = pan_sp["Y"][0, 160:]
    cover = float(np.mean((y_star >= out_s["ci_lo"]) & (y_star <= out_s["ci_hi"])))
    assert cover >= 0.5  # benign cell; loose single-panel bound
    # gated predictions match the legacy estimator path
    legacy = sf.spectral_sc(pan_sp["Y"][1:], pan_sp["Y"][0, :160],
                            sigma=1.0, c=81 / 160)
    assert np.allclose(out_s["pred_gated"], legacy)


def test_run_rep_schema_and_ordering():
    cell = dict(experiment="TEST", shard="t00", cell_id="c1", c=81 / 160,
                n=81, T0=160, T_post=60, r=1, m=3.75, arm="partial",
                theta=0.25, delta=None, noise="gaussian")
    rows = sf.run_rep(cell, seed=10_000, methods=sf.METHODS[:5],
                      diag_level="light")
    assert len(rows) == 6  # 5 methods + _diag
    by = {r["method"]: r for r in rows}
    assert set(by) == {"donor_mean", "scm_simplex", "ridge_sc",
                       "spectral_gated", "spectral_ungated", "_diag"}
    assert by["_diag"]["cv_rank"] in (0, 1, 2, 3, 4)
    assert 0.0 <= float(by["spectral_gated"]["ci_cover"]) <= 1.0
    # oracle floor respected: every method RMSE within sane band of sigma
    for m in ("donor_mean", "ridge_sc", "spectral_gated"):
        assert 0.9 < float(by[m]["rmse"]) < 2.0


def test_run_rep_full_diag_needs_cached_null():
    cell = dict(experiment="TEST", cell_id="full", c=81 / 160, n=81, T0=160,
                T_post=60, r=1, m=2.0, arm="partial", theta=0.25, delta=None)
    nz = sf.experiment.cell_null_z(cell, G=50)
    rows = sf.run_rep(cell, seed=5, methods=("spectral_gated",),
                      diag_level="full", null_z=nz)
    d = [r for r in rows if r["method"] == "_diag"][0]
    assert 0.0 < float(d["z_boot_p"]) <= 1.0
    assert 0.0 < float(d["z_tw_p"]) <= 1.0


def test_run_rep_orthogonal_arm_matches_alignment_none():
    cell = dict(experiment="TEST", cell_id="orth", c=81 / 160, n=81, T0=160,
                T_post=60, r=1, m=2.0, arm="orthogonal", theta=0.0,
                delta=None)
    rows = sf.run_rep(cell, seed=77, methods=("spectral_gated",),
                      diag_level="none")
    pan = sf.generate_panel(seed=77, **sf.experiment.panel_kwargs(cell))
    assert abs(pan["A"][0]).sum() == 0.0


def test_run_cell_chunking():
    cell = dict(experiment="TEST", cell_id="chunk", c=81 / 160, n=41, T0=80,
                T_post=40, r=1, m=2.0, arm="full", theta=1.0, delta=None)
    chunks = []
    sf.run_cell(cell, methods=("spectral_gated", "donor_mean"), reps=7,
                seed_base=10_000, on_chunk=chunks.append, diag_level="none")
    assert len(chunks) == 1  # reps < chunk size -> single flush
    total = sum(len(c) for c in chunks)
    assert total == 7 * 2


def test_spread_arm_dense_weak():
    cell = dict(experiment="C2II", cell_id="dense", c=81 / 160, n=81, T0=160,
                T_post=60, r=8, m=0.6, arm="spread", theta=1.0, delta=None)
    kw = sf.experiment.panel_kwargs(cell)
    assert kw["alignment"] == "all"
    assert np.allclose(kw["treated_share"], 1.0 / 8)
    rows = sf.run_rep(cell, seed=31, methods=("spectral_gated",),
                      diag_level="light")
    d = [r for r in rows if r["method"] == "_diag"][0]
    # individually subcritical spikes: silence gate should usually stay quiet
    assert int(d["k_selected"]) in (0, 1)
