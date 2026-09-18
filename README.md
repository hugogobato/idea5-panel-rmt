# causal-spectral-frontier

Spectral recoverability diagnostics for synthetic control and panel causal inference.

This repository develops, tests, and documents a diagnostic layer for panel causal estimators
(plain synthetic control, ridge and spectral variants, nuclear-norm matrix completion, and
synthetic difference-in-differences). It asks a question that pre-treatment fit alone cannot
answer: is the treated unit's counterfactual statistically recoverable from the donor pool at all?

The answer is built from the random-matrix theory of spiked matrices, in particular the
Baik-Ben Arous-Peche (BBP) phase transition. For a donor panel with aspect ratio `c = n_d / T_0`,
a latent factor is consistently estimable from the noisy data only if its spike strength
`s = ||a||^2 / sigma^2` exceeds `sqrt(c)`. The diagnostic reports three objects:

1. **Distance to the recoverability frontier**: the estimated spike strengths relative to the
   Marchenko-Pastur edge, summarized as `d = s / sqrt(c)` with bootstrap intervals.
2. **Alignment**: whether the treated unit's pre-period trajectory projects onto the recoverable
   factor directions beyond a calibrated pure-noise null.
3. **Transport**: whether the donor factor law is stable across the treatment date, tested with
   simulation-based nulls (circular shift and block permutation) that remain calibrated under
   serial correlation and heteroskedasticity.

**Status.** Methods paper in preparation. The simulation program is complete and preregistered;
the applied diagnostics ship as a casebook (see `docs/applied/applied_study_closure_memo.md`).
The claimed contribution is diagnostic and explanatory. No large estimator improvements are
claimed: in the experiments all incumbent methods sit on the same recoverability frontier, and
the spectral estimator's advantage is at most around 0.2 standard deviations at strong signal.
On canonical macro panels the frontier does not bind; that is itself one of the reported
findings, not a failure of the tool.

## Installation

Tested with Python 3.12 and the versions pinned in `requirements.txt`.

```bash
pip install -r requirements.txt
```

The importable packages live under `code/`, so either export the path or install in editable mode
once packaging is released:

```python
import sys
sys.path.insert(0, "code")
```

Run the test suite with:

```bash
python -m pytest code/tests -q
```

## Quickstart 1: a simulated panel near the frontier

The cell dictionary follows the preregistered experiment schema
(`config/results_schema.yaml`). `m = s / sqrt(c)` is the spike multiplier, so `m = 1` is the
frontier and `m > 1` is recoverable in population.

```python
import sys
sys.path.insert(0, "code")

from scm_frontier import METHODS, run_rep

cell = {
    "experiment": "quickstart", "cell_id": "c0.5_full_m2",
    "c": 121 / 240, "n": 121, "T0": 240, "T_post": 100, "r": 1,
    "m": 2.0, "arm": "full", "theta": 0.5,
}
rows = run_rep(cell, seed=10000, methods=METHODS, diag_level="light")
for row in rows:
    if row["method"] != "_diag":
        print(f'{row["method"]:>16}  rmse={row["rmse"]:.3f}')
```

This generates the panel internally from the seed, fits all seven estimator families
(donor mean, simplex synthetic control, ridge SC, gated and ungated spectral SC, nuclear-norm
matrix completion, SDID), and returns one row per method plus a `_diag` row carrying the
spectral diagnostics, the CV-rank comparator, and the classical pre-trend t-test.

## Quickstart 2: diagnose a real panel

The applied pipeline in `code/applications/frontier_d2.py` implements the frozen rules from the
Phase D preregistration. It expects a donor pre-period matrix and the treated pre-period row.
The full run includes donor and time-block bootstraps, so it takes a few minutes.

```python
import sys
sys.path.insert(0, "code")

import numpy as np
from applications.frontier_d2 import Seeds, analyze_unit

donors_pre = np.loadtxt("donors_pre.csv", delimiter=",")    # shape (n_d, T0)
treated_pre = np.loadtxt("treated_pre.csv", delimiter=",")  # shape (T0,)

seeds = Seeds(align=111, donor=222, time=333, perm=444)
out = analyze_unit(donors_pre, treated_pre, seeds.as_dict(), idx=0)

print("rank k          =", out["k"])
print("distance d      =", round(out["d"], 2))
print("alignment p     =", out["p_align"], "| ratio =", round(out["r_align"], 2))
print("label           =", out["label"])
print("donor CI for d  =", out["boot_donor"]["q05"], out["boot_donor"]["q95"])
```

The label comes from `classify()` in the same module:

| Label | Meaning |
|---|---|
| `RECOVERABLE` | at least one gated factor is retained, the treated unit is aligned, and the donor-bootstrap 5th percentile of `d` is at least 1 |
| `FRAGILE-INVISIBLE` | no donor factor clears the gated Marchenko-Pastur edge: the treated row carries no learnable signal |
| `FRAGILE-MISALIGNED` | factors are visible, but the treated unit does not project onto them (alignment `p >= 0.05`) |
| `FRAGILE-SUBEDGE` / `INCONCLUSIVE-SUBEDGE` | the bootstrap interval for `d` straddles the frontier |

Use `analyze_unit` for the full report (spectrum, `d` intervals, sensitivity windows, transport
control) and `spectrum_pipeline` / `alignment_test` for the individual pieces. All diagnostics use
pre-period donor data only; the treated post-period outcomes never enter model selection.

## Repository structure

```text
code/
  scm_frontier/        core package: DGPs, estimator zoo, spectral diagnostics
  applications/        real-panel ingestion, canonical SCM, distance-to-frontier pipeline
  tests/               pytest regression suite (43 tests)
scripts/               experiment drivers, shard merge and validation, figure builders,
                       Colab notebook generator, applied-study drivers
notebooks/             witness notebooks (Phase A adversarial checks)
colab/                 generated Colab shards + shard_manifest.yaml
data/raw/              immutable raw panels and checksums
results_c1/ results_c2/ results_d1/ results_d2/ results_e/ results_raw/
                       preregistered experiment outputs (parquet/csv/json)
figures/               regenerated figures plus memo input summaries
config/seeds.yaml      global seed registry
config/results_schema.yaml   experiment row schema
docs/                  research documentation, plans, preregistrations, gate memos
paper/powell_summary/  short outreach write-up for the applied collaboration
Applied_Study/         large third-party replication archives (kept local, gitignored)
requirements.txt
```

## Documentation

The full research trail is indexed in [`docs/README.md`](docs/README.md). The most useful entry
points for an external reader are:

1. [`docs/plan/Idea5_Panel_RMT_Research_Plan.md`](docs/plan/Idea5_Panel_RMT_Research_Plan.md):
   the phased program, claim ledger, gates, and give-up rules.
2. [`docs/model/model_card.md`](docs/model/model_card.md) and
   [`docs/model/frontier_ansatz.md`](docs/model/frontier_ansatz.md): the formal model and the
   deterministic-equivalent risk ansatz with its symbolically verified special cases.
3. [`docs/preregistrations/preregistration.md`](docs/preregistrations/preregistration.md):
   the frozen simulation design and pass rules.
4. [`docs/gates/gate_g3_memo.md`](docs/gates/gate_g3_memo.md): what the decisive simulation
   gate actually found, including the honest restrictions on the claims.
5. [`docs/applied/applied_study_closure_memo.md`](docs/applied/applied_study_closure_memo.md):
   the applied casebook scope and its negative results.
6. [`paper/powell_summary/powell_summary.pdf`](paper/powell_summary/powell_summary.pdf):
   a four-page intuitive summary with the main figures.

## Reproducing the results

Shards were run as self-contained Colab notebooks (fleet manual:
[`docs/colab/README_COLAB.md`](docs/colab/README_COLAB.md)). After downloading the raw shard
files into `results_raw/`:

```bash
python3 scripts/merge_shards.py     # validates every manifest entry (hash, schema, row counts)
python3 scripts/make_figures.py     # only runs after validation passes
python3 -m pytest code/tests -q     # regression suite (43 tests)
```

`scripts/merge_shards.py` refuses partial grids: missing or corrupted shards must be re-run
individually before figures are produced. Every table and figure states its seed range; the
global registry is `config/seeds.yaml`.

## Data access

Most panels are public and live in `data/raw/` with checksums. The Wisconsin handgun-suicide
application used by Powell (2025, *Journal of Applied Econometrics*) relies on NCHS-restricted
state-by-month NVSS mortality records that require a DUA and cannot be redistributed.
`Applied_Study/` holds the public replication archives downloaded from their DOIs; the restricted
inputs are not part of this repository, and any diagnostic run on them is designed to export
aggregate summaries only (see `docs/applied/applied_study_closure_memo.md`, Section 5).

## Citation

If you use this code or the diagnostic, please cite the repository until the paper is available:

```bibtex
@misc{gobatosouto2026spectralfrontier,
  author       = {Gobato Souto, Hugo},
  title        = {causal-spectral-frontier: Spectral recoverability diagnostics for panel causal inference},
  year         = {2026},
  howpublished = {\url{https://github.com/hugogobato/causal-spectral-frontier}},
  note         = {Methods paper in preparation}
}
```

## Contact

Hugo Gobato Souto. Website: <https://hugogobato.github.io/Hugo-Professional-Website/>.
Issues and pull requests are welcome; for research collaborations, please use the contact
information on the website.
