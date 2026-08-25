# Phase C on Google Colab — operating manual

**Preregistration:** `preregistration.md` (frozen 2026-08-24). Notebooks are
generated artifacts; do not edit their parameters. Source of truth:
`code/scm_frontier/` + `scripts/build_colab_notebooks.py`.

## Fleet (47 notebooks)

| family | notebooks | what it answers |
|---|---|---|
| `nb_c1_shard01..42_of42` | 42 | WP-C1 decisive grid (510 cells) + folded C2(i) null battery |
| `nb_c2ii_dense_weak` | 1 | WP-C2(ii) baseline-favorable dense-weak factors |
| `nb_c2iii_calibration` | 1 | WP-C2(iii) size/power under AR(1)/heteroskedastic; re-measures the bootstrap multiplier first |
| `nb_c2iv_structural_break` | 1 | WP-C2(iv) A4-violation detectability |
| `nb_c4_scaling` | 1 | WP-C4 spectral-family scaling to n = T0 = 5000 |
| `nb_onset_slice` | 1 | inherited falsifier: detectability onset -> m = 1 as n grows |
| `nb_c5a_shard01_of1`, `nb_c5b_shard01..04_of4` | 5 | C.5 confirmation: corrected kink estimator (fresh seeds 15000+i); theta=3 bite extension (seeds 16000+i) |
| `nb_c5c_diagv2_battery`, `nb_c5d_break_formal` | 2 | C.5 diagnostics-v2 calibration battery (z_shift/gate_lrv; seeds 17000+i); break formalization (post-window statistic; seeds 18000+i) |

Each shard targets <= 4.5 h worst case (loaded-unit projection; expect
~2 h on an unloaded Colab VM), checkpoints every 25 reps, and verifies its
own row count before offering downloads.

## Run

1. Upload any `.ipynb` from `colab/` to a Colab account and Run-all.
2. At the end it saves `idea5_out/<name>.csv.gz` + `<name>_meta.json`
   and triggers browser downloads (keep them; also check runtime persistence).
3. Put downloaded files here under `results_raw/<family>/`
   (`c1`, `c2ii`, `c2iii`, `c2iv`, `onset`, `c4`), keeping original filenames.

## Aggregate locally

```bash
python3 scripts/merge_shards.py          # validates EVERY manifest entry
python3 scripts/make_figures.py          # only runs if validation passed
```

`merge_shards.py` checks gzip integrity, sha256 against the notebook's own
meta JSON, schema, exact row counts, per-cell seed completeness, and recorded
cell errors. It writes `results_c1/risk_curves.parquet` and
`results_c2/*.parquet`. Figures are never produced from partial grids.

## Determinism

Panels are bitwise-reproducible (numpy PCG64, registered seeds,
10000+i per replication index). Floating-point output may differ at LAPACK
level across BLAS builds; this is expected and does not affect Monte Carlo
statements. Shard meta JSON records library versions.
