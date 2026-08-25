"""Build self-contained Colab notebooks for Phase C (preregistration Section 7).

Outputs (repo root):
  colab/nb_*.ipynb                 45 notebooks (40 main-grid shards + 5 specials)
  colab/shard_manifest.yaml        shard-to-cell map with expected row counts
  results_schema.yaml              row schema description

Run from repo root:  python3 scripts/build_colab_notebooks.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

REPS = 500
SEED_BASE = 10_000
N_MAIN_SHARDS = 42

GEOM = {
    0.25: (80, 320, 160),
    0.5: (113, 226, 113),
    1.0: (160, 160, 80),
    2.0: (226, 113, 56),
    4.0: (320, 80, 40),
}
PROD_NULL = (250, 250, 125)

FULL6 = ["donor_mean", "scm_simplex", "ridge_sc", "spectral_gated",
         "spectral_ungated", "mc_nn_cv"]
FULL7 = FULL6 + ["sdid"]


def np_linspace(a, b, n):
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


# ---------------------------------------------------------------- grid ----

def build_main_cells():
    m_grid = sorted({round(v, 10) for v in np_linspace(0.2, 3.0, 15) + [0.9, 1.1]})
    arms = [("orthogonal", 0.0), ("partial", 0.25), ("full", 1.0)]
    cells = []
    for c, (n, T0, Tp) in GEOM.items():
        for m in m_grid:
            for arm, th in arms:
                for r in (1, 3):
                    methods = list(FULL6)
                    if c in (0.5, 2.0) and arm == "full" and r == 1:
                        methods.append("sdid")
                    cells.append(dict(
                        experiment="c1", cell_id=f"c{c}_m{m:.2f}_{arm}_r{r}",
                        c=c, n=n, T0=T0, T_post=Tp, r=r, m=m,
                        arm=arm, theta=th, delta=None,
                        methods=methods, diag="light", reps=REPS))
    for c, (n, T0, Tp) in GEOM.items():          # C2(i) folded in
        cells.append(dict(
            experiment="c2i", cell_id=f"null_c{c}",
            c=c, n=n, T0=T0, T_post=Tp, r=0, m=0.0,
            arm="orthogonal", theta=0.0, delta=None,
            methods=list(FULL6), diag="full", reps=REPS))
    cells.append(dict(
        experiment="c2i", cell_id="null_prod250",
        c=1.0, n=PROD_NULL[0], T0=PROD_NULL[1], T_post=PROD_NULL[2],
        r=0, m=0.0, arm="orthogonal", theta=0.0, delta=None,
        methods=list(FULL7), diag="full", reps=REPS))
    return cells


def build_c2ii_cells():
    n, T0, Tp = GEOM[1.0]
    return [dict(experiment="c2ii", cell_id=f"dense_r{r}_m{m}",
                 c=1.0, n=n, T0=T0, T_post=Tp, r=r, m=m,
                 arm="spread", theta=1.0, delta=None,
                 methods=list(FULL6), diag="light", reps=REPS)
            for r in (8, 16) for m in (0.6, 0.8)]


def build_c2iv_cells():
    n, T0, Tp = GEOM[1.0]
    spec = [("break_d1_aligned", 1.0, "full", 1.0),
            ("break_d2_aligned", 2.0, "full", 1.0),
            ("break_d2_orthogonal", 2.0, "orthogonal", 0.0),
            ("control_d0_aligned", None, "full", 1.0)]
    return [dict(experiment="c2iv", cell_id=cid, c=1.0, n=n, T0=T0,
                 T_post=Tp, r=1, m=2.0, arm=arm, theta=th, delta=d,
                 methods=list(FULL7), diag="full", reps=REPS)
            for cid, d, arm, th in spec]


def build_c2iii_cells():
    laws = [("gaussian", dict(noise="gaussian")),
            ("ar1_r03", dict(noise="ar1", rho=0.3)),
            ("ar1_r07", dict(noise="ar1", rho=0.7)),
            ("het", dict(noise="heteroskedastic", het_ratio=4.0))]
    n, T0, Tp = GEOM[1.0]
    cells = []
    for lname, lk in laws:
        for kind, m, arm, th in [("null", 0.0, "orthogonal", 0.0),
                                 ("spiked", 1.6, "full", 1.0)]:
            cfg = dict(experiment="c2iii", cell_id=f"{kind}_{lname}",
                       c=1.0, n=n, T0=T0, T_post=Tp, r=1, m=m, arm=arm,
                       theta=th, delta=None,
                       methods=["spectral_gated", "spectral_ungated",
                                "ridge_sc"],
                       diag="full", reps=REPS)
            cfg.update(lk)
            cells.append(cfg)
    return cells


# ----------------------------- C5 addendum families ------------------------

C5A_SEED_BASE = 15_000
C5B_SEED_BASE = 16_000
C5C_SEED_BASE = 17_000
C5D_SEED_BASE = 18_000


def build_c5a_cells():
    """Kink confirmation: corrected estimator, fresh seeds, spectral-only."""
    m_grid = [round(0.6 + i * (1.0 / 20), 10) for i in range(21)]
    cells = []
    for c, (n, T0, Tp) in GEOM.items():
        for m in m_grid:
            cells.append(dict(
                experiment="c5a", cell_id=f"c{c}_m{m:.2f}_full_r1",
                c=c, n=n, T0=T0, T_post=Tp, r=1, m=m,
                arm="full", theta=1.0, delta=None,
                methods=["spectral_gated", "donor_mean"], diag="light",
                reps=REPS))
    return cells


def build_c5b_cells():
    """Bite extension: theta = 3 arm (predicted plateau sqrt(4) = 2 sigma)."""
    m_grid = sorted({round(v, 10) for v in np_linspace(0.2, 3.0, 15) + [0.9, 1.1]})
    cells = []
    for c in (0.5, 1.0, 2.0):
        n, T0, Tp = GEOM[c]
        for m in m_grid:
            cells.append(dict(
                experiment="c5b", cell_id=f"c{c}_m{m:.2f}_full_r1",
                c=c, n=n, T0=T0, T_post=Tp, r=1, m=m,
                arm="full", theta=3.0, delta=None,
                methods=list(FULL6), diag="light", reps=400))
    return cells


# ----------------------------------------------------------- cost model ----

def rep_seconds(cell):
    """Loaded-unit seconds per replication (WP-B3 pilot exponents)."""
    q = (cell["n"] * (cell["T0"] + cell["T_post"])) / 93750.0
    meth = cell.get("methods", [])
    s = 0.001
    if "scm_simplex" in meth:
        s += 4.345 * q ** 1.79
    if "ridge_sc" in meth:
        s += 0.216 * q ** 1.31
    if any(x.startswith("spectral") for x in meth) or cell.get("diag") != "none":
        s += 0.03 * q ** 1.31
    if "mc_nn_cv" in meth:
        s += 4.788 * q ** 1.47
    if "sdid" in meth:
        s += 11.02 * q ** 1.67
    if cell.get("diag") == "light":
        s += 0.02
    elif cell.get("diag") == "full":
        nd = cell["n"] - 1
        tb = cell["T0"] - min(cell["T_post"], cell["T0"] // 2)
        s += 200 * 6.5e-7 * max(nd * tb, 1)
    return s


def cell_seconds(cell):
    return rep_seconds(cell) * cell["reps"] + 20


def lpt_balance(cells, k):
    loads = [0.0] * k
    shards = [[] for _ in range(k)]
    for cell in sorted(cells, key=cell_seconds, reverse=True):
        i = min(range(k), key=lambda j: loads[j])
        shards[i].append(cell)
        loads[i] += cell_seconds(cell)
    order = sorted(range(k), key=lambda j: -loads[j])
    return [shards[i] for i in order], [loads[i] for i in order]


# ------------------------------------------------------------ notebooks ----

def flat_module_source():
    parts = []
    for fname in ("dgps.py", "estimators.py", "diagnostics.py",
                  "experiment.py"):
        src = (ROOT / "code" / "scm_frontier" / fname).read_text()
        out, paren_depth = [], 0
        for line in src.splitlines():
            if line.startswith("from __future__"):
                continue
            if paren_depth:
                paren_depth += line.count("(") - line.count(")")
                continue
            if line.startswith("from ."):
                paren_depth = line.count("(") - line.count(")")
                if paren_depth <= 0:
                    paren_depth = 0
                continue
            out.append(line)
        parts.append(f"# ===== flattened from code/scm_frontier/{fname} "
                     "=====\n" + "\n".join(out).strip() + "\n")
    pre = ('"""scm_frontier_flat: single-file Phase C module '
           '(auto-generated; DO NOT EDIT)."""\n\n'
           "from __future__ import annotations\n\n")
    return pre + "\n\n".join(parts)


def md_cell(text):
    return {"cell_type": "markdown", "metadata": {},
            "source": text.splitlines(keepends=True)}


def code_cell(text):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": text.splitlines(keepends=True)}


def new_nb(name):
    return {"cells": [], "metadata": {"colab": {"name": name},
            "kernelspec": {"name": "python3", "display_name": "Python 3"}},
            "nbformat": 4, "nbformat_minor": 5}


SETUP_CODE = '''import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy

print("python", platform.python_version(), "| numpy", np.__version__,
      "| scipy", scipy.__version__)
'''

RUNNER_LIB = '''

def make_flush(csv_path, cols):
    import csv as _csv
    def flush(rows):
        new = not Path(csv_path).exists()
        with open(csv_path, "a", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=cols)
            if new:
                w.writeheader()
            w.writerows(rows)
    return flush
'''

MAIN_LOOP = '''
CSV = OUT / f"{NB_NAME}.csv"
flush = make_flush(CSV, sf.ROW_COLS)

META = {"notebook": NB_NAME, "family": FAMILY, "shard": SHARD_ID,
        "cells": [{"cell_id": c["cell_id"], "experiment": c["experiment"]}
                  for c in CELLS],
        "reps_per_cell": REPS, "seed_base": SEED_BASE,
        "versions": {"python": platform.python_version(),
                     "numpy": np.__version__},
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "cell_errors": {}, "rows_expected": None}

META["rows_expected"] = sum(
    c["reps"] * (len(c["methods"]) + (1 if c["diag"] != "none" else 0))
    for c in CELLS)
print(f"{NB_NAME}: {len(CELLS)} cells, {META['rows_expected']} expected rows")

for ci, cell in enumerate(CELLS):
    t0 = time.perf_counter()
    try:
        sf.run_cell(cell, tuple(cell["methods"]), cell["reps"], SEED_BASE,
                    on_chunk=flush, diag_level=cell["diag"])
    except Exception as exc:  # record, continue with remaining cells
        META["cell_errors"][cell["cell_id"]] = repr(exc)
        print(f"[{ci + 1}/{len(CELLS)}] {cell['cell_id']} ERROR {exc!r}")
        continue
    print(f"[{ci + 1}/{len(CELLS)}] {cell['cell_id']} done in "
          f"{time.perf_counter() - t0:.1f}s")
print(f"all cells done in {time.perf_counter() - T_NB:.0f}s")
'''

FINALIZE_CODE = '''import gzip
import shutil

rows_written = sum(1 for _ in open(CSV)) - 1
META["rows_written"] = rows_written
META["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
META["wall_s"] = round(time.perf_counter() - T_NB, 1)
if META["cell_errors"]:
    print("WARNING: cells with errors ->", META["cell_errors"])
assert rows_written == META["rows_expected"] or META["cell_errors"], (
    f"incomplete without recorded errors: {rows_written} "
    f"!= {META['rows_expected']}")
META["csv_sha256"] = hashlib.sha256(open(CSV, "rb").read()).hexdigest()
gz = CSV.with_suffix(".csv.gz")
with open(CSV, "rb") as fin, gzip.open(gz, "wb") as fout:
    shutil.copyfileobj(fin, fout)
meta_path = OUT / f"{NB_NAME}_meta.json"
meta_path.write_text(json.dumps(META, indent=1))
print("rows:", rows_written, "| sha256:", META["csv_sha256"][:16], "...")

try:
    from google.colab import files
    files.download(str(gz))
    files.download(str(meta_path))
    print("Downloaded:", gz.name, meta_path.name)
except Exception as e:
    print("(Not on Colab / download skipped):", e)
'''


def module_cell():
    return code_cell("# Auto-generated single-file module (source of truth:\n"
                     "# code/scm_frontier/* via scripts/build_colab_notebooks.py)\n"
                     "%%writefile scm_frontier_flat.py\n" + flat_module_source())


def config_cell(nb_name, family, shard_id, cells):
    return code_cell(
        "import scm_frontier_flat as sf\n"
        f'NB_NAME = "{nb_name}"\nFAMILY = "{family}"\nSHARD_ID = "{shard_id}"\n'
        f'REPS = {REPS}\nSEED_BASE = {SEED_BASE}\n'
        "CELLS = json.loads(r'''\n" + json.dumps(cells, indent=1) + "\n''')\n"
        "OUT = Path('idea5_out'); OUT.mkdir(exist_ok=True)\n"
        "T_NB = time.perf_counter()\n" + RUNNER_LIB)


def build_experiment_nb(nb_name, family, shard_id, cells, title_md):
    nb = new_nb(nb_name)
    nb["cells"].append(md_cell(title_md))
    nb["cells"].append(module_cell())
    nb["cells"].append(code_cell(SETUP_CODE))
    nb["cells"].append(config_cell(nb_name, family, shard_id, cells))
    nb["cells"].append(code_cell(MAIN_LOOP.strip("\n")))
    nb["cells"].append(code_cell(FINALIZE_CODE))
    return nb


C2III_MEASURE_CELL = '''# G2 obligation 5: re-measure the bootstrap-cost multiplier BEFORE any
# calibration output is produced (preregistration Section 4.3).
_geo = CELLS[0]
_rng = np.random.default_rng(12345)
_nd = _geo["n"] - 1
_Tp = min(_geo["T_post"], _geo["T0"] // 2)
_Tb = _geo["T0"] - _Tp
_Yp = _rng.normal(size=(_nd, _geo["T0"]))
_t0 = time.perf_counter()
for _ in range(5):
    sf.resid_statistic(_Yp[:, :_Tb], _Yp[:, _Tb:], 1.0, _nd / _Tb)
_base = (time.perf_counter() - _t0) / 5
_t0 = time.perf_counter()
for _ in range(5):
    sf.z_boot_pvalue(_Yp, 1.0, _Tp, B=200)
_boot = (time.perf_counter() - _t0) / 5
BOOT_MULTIPLIER = _boot / max(_base, 1e-9)
print(f"bootstrap multiplier (B=200 vs obs-only): {BOOT_MULTIPLIER:.1f}x "
      f"(base {_base * 1e3:.1f} ms, boot {_boot:.2f} s)")
'''


C5C_CELLS = [
    dict(law="gaussian", lk=dict(noise="gaussian"), state="null"),
    dict(law="gaussian", lk=dict(noise="gaussian"), state="spiked"),
    dict(law="gaussian", lk=dict(noise="gaussian"), state="break_d1"),
    dict(law="gaussian", lk=dict(noise="gaussian"), state="break_d2"),
    dict(law="ar1_r03", lk=dict(noise="ar1", rho=0.3), state="null"),
    dict(law="ar1_r03", lk=dict(noise="ar1", rho=0.3), state="spiked"),
    dict(law="ar1_r03", lk=dict(noise="ar1", rho=0.3), state="break_d1"),
    dict(law="ar1_r03", lk=dict(noise="ar1", rho=0.3), state="break_d2"),
    dict(law="ar1_r07", lk=dict(noise="ar1", rho=0.7), state="null"),
    dict(law="ar1_r07", lk=dict(noise="ar1", rho=0.7), state="spiked"),
    dict(law="ar1_r07", lk=dict(noise="ar1", rho=0.7), state="break_d1"),
    dict(law="ar1_r07", lk=dict(noise="ar1", rho=0.7), state="break_d2"),
    dict(law="het", lk=dict(noise="heteroskedastic", het_ratio=4.0),
         state="null"),
    dict(law="het", lk=dict(noise="heteroskedastic", het_ratio=4.0),
         state="spiked"),
    dict(law="het", lk=dict(noise="heteroskedastic", het_ratio=4.0),
         state="break_d1"),
    dict(law="het", lk=dict(noise="heteroskedastic", het_ratio=4.0),
         state="break_d2"),
]

C5C_BODY = SETUP_CODE + '''
import scm_frontier_flat as sf

OUT = Path("idea5_out"); OUT.mkdir(exist_ok=True)
CSV = OUT / "c5c_diagv2.csv"

def make_flush(csv_path, cols):
    import csv as _csv
    def flush(rows):
        new = not Path(csv_path).exists()
        with open(csv_path, "a", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=cols)
            if new:
                w.writeheader()
            w.writerows(rows)
    return flush

flush = make_flush(CSV, ["cell_id", "law", "state", "rep_seed",
                         "z_shift_p", "z_perm_p", "z_tw_p", "trend_p",
                         "gate_lrv_k", "gate_mp_k"])
CELLS = json.loads(r\'\'\'
''' + json.dumps(C5C_CELLS, indent=1) + '''
\'\'\')
N, T0, TPOST = 160, 160, 80
REPS_C5C, SEED0 = 300, ''' f"{C5C_SEED_BASE}" '''
T_NB = time.perf_counter()
TB = T0 - TPOST // 2
null_z = sf.simulated_null_z(N - 1, TB, TPOST // 2, 1.0, (N - 1) / TB,
                             G=300, seed=7_770_001)
for cell in CELLS:
    law, lk, state = cell["law"], cell["lk"], cell["state"]
    cid = f"{state}_{law}"
    for j in range(REPS_C5C):
        pan_kw = dict(n=N, T0=T0, T_post=TPOST, sigma=1.0,
                      noise=lk.get("noise", "gaussian"),
                      rho=lk.get("rho", 0.5),
                      het_ratio=lk.get("het_ratio", 4.0))
        if state == "null":
            pan_kw.update(r=0, spike_strengths=(), treated_share=(),
                          alignment="none")
        else:
            s_spk = 1.6 * np.sqrt(N / T0)
            pan_kw.update(r=1, spike_strengths=(s_spk,),
                          treated_share=(1.0,), alignment="first")
        pan_kw["structural_break"] = {"break_d1": 1.0,
                                      "break_d2": 2.0}.get(state)
        pan = sf.generate_panel(seed=SEED0 + j, **pan_kw)
        Ypre = pan["Y"][1:, :T0]
        p_shift, _, _ = sf.z_shift_pvalue(Ypre, 1.0, TPOST, B=200)
        p_perm, _, _ = sf.z_perm_pvalue(Ypre, 1.0, TPOST, B=200, block=20)
        z_obs, _ = sf.resid_statistic(Ypre[:, :TB], Ypre[:, TB:], 1.0,
                                      (N - 1) / TB)
        p_ztw = sf.z_tw_pvalue(z_obs, null_z)
        trend_p = sf.classical_trend_ttest(pan["Y"][0, :T0])
        flush([dict(cell_id=cid, law=law, state=state,
                    rep_seed=SEED0 + j, z_shift_p=round(p_shift, 6),
                    z_perm_p=round(p_perm, 6), z_tw_p=round(p_ztw, 6),
                    trend_p=round(trend_p, 6),
                    gate_lrv_k=sf.gate_lrv(Ypre),
                    gate_mp_k=int(sf.gated_rank(sf.scree(Ypre), 1.0,
                                                N / T0)))])
    print(f"{cid} done ({time.perf_counter() - T_NB:.0f}s)")

(Path("idea5_out") / "c5c_done.json").write_text(json.dumps(
    {"rows_expected": len(CELLS) * REPS_C5C}, indent=1))
try:
    from google.colab import files
    files.download(str(CSV)); files.download(str(OUT / "c5c_done.json"))
except Exception as e:
    print("(Not on Colab / download skipped):", e)
'''


def build_c5c_nb():
    nb = new_nb("nb_c5c_diagv2_battery")
    nb["cells"].append(md_cell(
        "# C5c diagnostics-v2 calibration battery\n\n"
        "Preregistration_c5_addendum C5c. Laws x states; z_shift primary,\n"
        "z_perm20 secondary, z_tw/t-test references, gate_lrv vs gate_mp.\n"))
    nb["cells"].append(module_cell())
    nb["cells"].append(code_cell(C5C_BODY.strip("\n")))
    return nb


C5D_BODY = SETUP_CODE + '''
import scm_frontier_flat as sf

OUT = Path("idea5_out"); OUT.mkdir(exist_ok=True)
CSV = OUT / "c5d_break_formal.csv"

def make_flush(csv_path, cols):
    import csv as _csv
    def flush(rows):
        new = not Path(csv_path).exists()
        with open(csv_path, "a", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=cols)
            if new:
                w.writeheader()
            w.writerows(rows)
    return flush

flush = make_flush(CSV, ["cell_id", "rep_seed", "p_post", "z_post",
                         "rmse_spectral"])
CELLS = json.loads(r\'\'\'
''' + json.dumps([
    dict(cell_id=f"delta{d}", delta=d, arm=("none" if d == "orth" else "first"))
    for d in [0, 0.5, 1.0, 2.0]
] + [dict(cell_id="delta2_orthogonal", delta=2.0, arm="none")], indent=1) + '''
\'\'\')
N, T0, TPOST, REPS_D, SEED0 = 160, 160, 80, 400, ''' f"{C5D_SEED_BASE}" '''
T_NB = time.perf_counter()
for cell in CELLS:
    cid = str(cell["cell_id"])
    pan_null = None
    null_z = sf.simulated_null_z(N - 1, T0, TPOST, 1.0, (N - 1) / T0,
                                 G=300, seed=7_770_001)
    align = cell["arm"]
    th = 0.0 if align == "none" else 1.0
    for j in range(REPS_D):
        pan = sf.generate_panel(seed=SEED0 + j, n=N, T0=T0, T_post=TPOST,
                                r=1, spike_strengths=(2.0,),
                                treated_share=(th,), alignment=align,
                                sigma=1.0, structural_break=cell["delta"])
        Y = pan["Y"]
        out = sf.pre_trends_post_test(Y[1:, :T0], Y[1:, T0:], 1.0,
                                      null_z=null_z)
        spec = sf.spectral_sc_full(Y[1:], Y[0, :T0], sigma=1.0, c=N / T0)
        y_star = Y[0, T0:]
        rmse = float(np.sqrt(np.mean((spec["pred_gated"] - y_star) ** 2)))
        flush([dict(cell_id=cid, rep_seed=SEED0 + j,
                    p_post=round(out["p"], 6), z_post=round(out["z"], 6),
                    rmse_spectral=round(rmse, 6))])
    print(f"{cid} done ({time.perf_counter() - T_NB:.0f}s)")
try:
    from google.colab import files
    files.download(str(CSV))
except Exception as e:
    print("(Not on Colab / download skipped):", e)
'''


def build_c5d_nb():
    nb = new_nb("nb_c5d_break_formal")
    nb["cells"].append(md_cell(
        "# C5d break-detection formalization\n\n"
        "Preregistration_c5_addendum C5d: post-window statistic on real\n"
        "donor post outcomes, simulated iid null; delta dose-response.\n"))
    nb["cells"].append(module_cell())
    nb["cells"].append(code_cell(C5D_BODY.strip("\n")))
    return nb


def build_onset_nb():
    nb_name = "nb_onset_slice"
    sizes = [(81, 81), (121, 121), (161, 161), (241, 241), (361, 361),
             (541, 541)]
    body = SETUP_CODE + '''
import scm_frontier_flat as sf
from scipy import stats

SIGMA = 1.0
KAPPA = 3.0
M_GRID = np.linspace(0.55, 1.65, 23)   # Witness-1 grid (prereg Sec 6.3)
R_GRID, R_NULL = 300, 600
SIZES = ''' + json.dumps(sizes) + '''
SEED_SCAN, SEED_NULL = 51101, 51102
OUT = Path("idea5_out"); OUT.mkdir(exist_ok=True)
CSV = OUT / "onset_slice.csv"

def make_flush(csv_path, cols):
    import csv as _csv
    def flush(rows):
        new = not Path(csv_path).exists()
        with open(csv_path, "a", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=cols)
            if new:
                w.writeheader()
            w.writerows(rows)
    return flush

flush = make_flush(CSV, ["n", "T0", "m", "ks_medp", "gate_rate",
                         "top_mean"])
rows_all = []
rng_scan = np.random.default_rng(SEED_SCAN)

for (n, T0) in SIZES:
    c = n / T0
    rng_null = np.random.default_rng(SEED_NULL + n)
    null_q = np.empty(R_NULL)
    for j in range(R_NULL):
        E = rng_null.normal(0.0, SIGMA, size=(n, T0))
        _, _, Vt = np.linalg.svd(E, full_matrices=False)
        null_q[j] = float(E[0] @ Vt[0]) ** 2
    pools = np.array_split(null_q, 6)
    for gi, m in enumerate(M_GRID):
        s = m * np.sqrt(c)
        sd2 = s * SIGMA ** 2 / (n - 1)
        A = np.concatenate(([KAPPA * np.sqrt(sd2)],
                            rng_scan.normal(0.0, np.sqrt(sd2), size=n - 1)))
        qs, gates, tops = (np.empty(R_GRID) for _ in range(3))
        for j in range(R_GRID):
            f = rng_scan.normal(0.0, 1.0, size=T0)
            E = rng_scan.normal(0.0, SIGMA, size=(n, T0))
            Y = np.outer(A, f) + E
            U, d, Vt = np.linalg.svd(Y, full_matrices=False)
            qs[j] = float(Y[0] @ Vt[0]) ** 2
            evals = np.sort(d ** 2 / T0)[::-1]
            gates[j] = sf.gated_rank(evals, SIGMA, c, k_max=4) > 0
            tops[j] = evals[0]
        medp = float(np.median([stats.ks_2samp(qs, p).pvalue
                                for p in pools]))
        row = dict(n=n, T0=T0, m=float(m), ks_medp=medp,
                   gate_rate=float(np.mean(gates)),
                   top_mean=float(np.mean(tops)))
        flush([row])
        rows_all.append(row)
    print(f"size ({n},{T0}) done; last medp={medp:.4g}")

# verdicts (preregistration Section 6.3)
onsets = {}
for (n, T0) in SIZES:
    rs = [r for r in rows_all if r["n"] == n]
    below = [r["m"] for r in rs if r["ks_medp"] < 0.01]
    onsets[f"{n}x{T0}"] = float(below[0]) if below else None
largest = SIZES[-1]
on_last = onsets[f"{largest[0]}x{largest[1]}"]
ONSET_PASS = bool(on_last is not None and 0.90 <= on_last <= 1.10)
verdict = {"onset_by_size": onsets, "onset_pass": ONSET_PASS,
           "rule": "onset(largest size) in [0.90, 1.10]"}
(Path("idea5_out") / "onset_verdict.json").write_text(
    json.dumps(verdict, indent=1))
print(json.dumps(verdict, indent=1))

try:
    from google.colab import files
    files.download(str(CSV))
    files.download(str(Path("idea5_out") / "onset_verdict.json"))
    print("Downloaded:", CSV.name)
except Exception as e:
    print("(Not on Colab / download skipped):", e)
'''
    nb = new_nb(nb_name)
    nb["cells"].append(md_cell(
        "# Inherited falsifier: detectability-onset convergence to m = 1\n\n"
        "Preregistration Section 6.3. W1 squared-coefficient statistic,\n"
        "Witness-1 m-grid, six sizes at c = 1 with joint (n, T0) growth.\n"
        "PASS iff onset at the largest size lies in [0.90, 1.10].\n"))
    nb["cells"].append(module_cell())
    nb["cells"].append(code_cell(body.strip("\n")))
    return nb


def build_scaling_nb():
    body = SETUP_CODE + '''
import resource
import scm_frontier_flat as sf

OUT = Path("idea5_out"); OUT.mkdir(exist_ok=True)
CSV = OUT / "scaling_probe.csv"
CASES = [((500, 500, 250), 3), ((1000, 1000, 500), 3),
         ((2000, 2000, 1000), 3), ((4000, 4000, 2000), 2),
         ((5000, 5000, 2500), 1)]

def make_flush(csv_path, cols):
    import csv as _csv
    def flush(rows):
        new = not Path(csv_path).exists()
        with open(csv_path, "a", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=cols)
            if new:
                w.writeheader()
            w.writerows(rows)
    return flush

flush = make_flush(CSV, ["n", "T0", "rep", "method", "wall_s", "rss_gib"])
for (n, T0, Tp), reps in CASES:
    c = n / T0
    for j in range(reps):
        pan = sf.generate_panel(n=n, T0=T0, T_post=Tp, r=1,
                                spike_strengths=(2.0,), treated_share=(0.5,),
                                alignment="first", sigma=1.0, seed=32000 + j)
        Y = pan["Y"]
        donors, y1_pre = Y[1:], Y[0, :T0]
        t0 = time.perf_counter()
        out = sf.spectral_sc_full(donors, y1_pre, sigma=1.0, c=c)
        dt = time.perf_counter() - t0
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
        flush([dict(n=n, T0=T0, rep=j, method="spectral_sc_full",
                    wall_s=round(dt, 4), rss_gib=round(rss, 4))])
        print(f"n={n}: spectral {dt:.2f}s rss={rss:.2f} GiB")

import json as _json
summary = []
for (n, T0, Tp), reps in CASES:
    pass  # exponents fitted locally by scripts/make_figures.py
Path("idea5_out/scaling_summary.json").write_text(_json.dumps(
    {"cases": [list(map(list, c)) for c, _ in CASES]}, indent=1))

try:
    from google.colab import files
    files.download(str(CSV))
    print("Downloaded:", CSV.name)
except Exception as e:
    print("(Not on Colab / download skipped):", e)
'''
    nb = new_nb("nb_c4_scaling")
    nb["cells"].append(md_cell(
        "# WP-C4 scaling study (preregistration Section 7)\n\n"
        "Spectral-family timing/memory up to n = T0 = 5000; simplex-family\n"
        "costs are extrapolated analytically from WP-B3 exponents.\n"))
    nb["cells"].append(module_cell())
    nb["cells"].append(code_cell(body.strip("\n")))
    return nb


def write_nb(path, nb):
    path.write_text(json.dumps(nb, indent=1))


def main():
    out_dir = ROOT / "colab"
    out_dir.mkdir(exist_ok=True)
    manifest = {"created": "2026-08-24", "seed_base": SEED_BASE,
                "notebooks": []}

    # ---- main fleet -------------------------------------------------------
    cells = build_main_cells()
    shards, loads = lpt_balance(cells, N_MAIN_SHARDS)
    print(f"main fleet: {len(cells)} cells, shard hours "
          f"(loaded-units): max {loads[0] / 3600:.2f}, "
          f"min {loads[-1] / 3600:.2f}")
    for i, shard in enumerate(shards, start=1):
        name = f"nb_c1_shard{i:02d}_of{N_MAIN_SHARDS}"
        title = (f"# WP-C1 decisive grid, shard {i:02d}/{N_MAIN_SHARDS}\n\n"
                 "Frozen under preregistration.md (2026-08-24). Contains\n"
                 "c1 grid cells plus folded C2(i) null cells; checkpoints\n"
                 "every 25 reps; do not edit parameters.\n\n"
                 f"Cells ({len(shard)}):\n" +
                 "\n".join(f"- {c['cell_id']} ({c['experiment']})"
                           for c in shard))
        write_nb(out_dir / f"{name}.ipynb",
                 build_experiment_nb(name, "c1", f"{i:02d}", shard, title))
        manifest["notebooks"].append(dict(
            family="c1", file=f"{name}.ipynb", shard=i,
            cells=[dict(cell_id=c["cell_id"], experiment=c["experiment"],
                        reps=c["reps"], methods=c["methods"],
                        diag=c["diag"]) for c in shard]))

    # ---- specials ---------------------------------------------------------
    specials = [
        ("nb_c2ii_dense_weak", "c2ii", build_c2ii_cells(),
         "# WP-C2(ii) baseline-favorable battery (dense weak factors)\n\n"
         "Preregistration Section 4.2. Individually subcritical spikes;\n"
         "reports honestly whichever family wins.\n"),
        ("nb_c2iv_structural_break", "c2iv", build_c2iv_cells(),
         "# WP-C2(iv) structural-break detectability (A4 violation)\n\n"
         "Preregistration Section 4.4. Z_boot power >= 80% in break cells;\n"
         "control size <= 8%; every estimator degrades vs its control pair.\n"),
    ]
    for name, fam, cs, title in specials:
        write_nb(out_dir / f"{name}.ipynb",
                 build_experiment_nb(name, fam, "single", cs, title))
        manifest["notebooks"].append(dict(
            family=fam, file=f"{name}.ipynb", shard="single",
            cells=[dict(cell_id=c["cell_id"], experiment=c["experiment"],
                        reps=c["reps"], methods=c["methods"],
                        diag=c["diag"]) for c in cs]))

    # ---- C5a/C5b experiment families (fresh-seed confirmation grids) -------
    c5a_cells = build_c5a_cells()
    c5b_cells = build_c5b_cells()
    for fam_name, cs, sbase in (("c5a", c5a_cells, C5A_SEED_BASE),
                                ("c5b", c5b_cells, C5B_SEED_BASE)):
        shards, loads = lpt_balance(cs, 4 if fam_name == "c5b" else 1)
        n_sh = len(shards)
        for i, shard in enumerate(shards, start=1):
            name = f"nb_{fam_name}_shard{i:02d}_of{n_sh}"
            title = (f"# {fam_name.upper()} confirmation grid, shard "
                     f"{i:02d}/{n_sh}\n\nFrozen under "
                     "preregistration_c5_addendum.md; fresh seed range.\n")
            nb = build_experiment_nb(name, fam_name, f"{i:02d}", shard, title)
            # override the seed base inside the config cell
            for cell in nb["cells"]:
                if cell["cell_type"] == "code":
                    src = "".join(cell["source"])
                    if "SEED_BASE" in src and "CELLS =" in src:
                        src = src.replace(f"SEED_BASE = {SEED_BASE}",
                                          f"SEED_BASE = {sbase}")
                        cell["source"] = src.splitlines(keepends=True)
            write_nb(out_dir / f"{name}.ipynb", nb)
            manifest["notebooks"].append(dict(
                family=fam_name, file=f"{name}.ipynb", shard=i,
                seed_base=sbase,
                cells=[dict(cell_id=c["cell_id"], experiment=c["experiment"],
                            reps=c["reps"], methods=c["methods"],
                            diag=c["diag"]) for c in shard]))

    name, fam = "nb_c2iii_calibration", "c2iii"
    cs = build_c2iii_cells()
    nb = new_nb(name)
    nb["cells"].append(md_cell(
        "# WP-C2(iii) calibration battery (AR(1)/heteroskedastic)\n\n"
        "Preregistration Section 4.3. Sizes at nominal 5% (Z_tw, Z_boot,\n"
        "classical trend t-test); power in spiked cells; bootstrap-cost\n"
        "multiplier re-measured before outputs (G2 obligation 5).\n"))
    nb["cells"].append(module_cell())
    nb["cells"].append(code_cell(SETUP_CODE))
    nb["cells"].append(config_cell(name, fam, "single", cs))
    nb["cells"].append(code_cell(C2III_MEASURE_CELL))
    nb["cells"].append(code_cell(MAIN_LOOP.strip("\n")))
    nb["cells"].append(code_cell(FINALIZE_CODE))
    write_nb(out_dir / f"{name}.ipynb", nb)
    manifest["notebooks"].append(dict(
        family=fam, file=f"{name}.ipynb", shard="single",
        cells=[dict(cell_id=c["cell_id"], experiment=c["experiment"],
                    reps=c["reps"], methods=c["methods"], diag=c["diag"])
               for c in cs], extra="bootstrap multiplier measured first"))

    for builder, fam in ((build_onset_nb, "onset"), (build_scaling_nb, "c4"),
                         (build_c5c_nb, "c5c"), (build_c5d_nb, "c5d")):
        nb = builder()
        fname = nb["metadata"]["colab"]["name"] + ".ipynb"
        write_nb(out_dir / fname, nb)
        manifest["notebooks"].append(dict(family=fam, file=fname,
                                          shard="single", cells="special"))

    for nb_meta in manifest["notebooks"]:
        p = out_dir / nb_meta["file"]
        nb_meta["sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()

    (out_dir / "shard_manifest.yaml").write_text(json.dumps(manifest, indent=1))

    schema = {
        "row_schema": {
            "meta_cols": list(__import__("scm_frontier").ROW_COLS[:16]),
            "value_cols": list(__import__("scm_frontier").ROW_COLS[16:]),
            "units": {"rmse": "sigma", "att_bias": "sigma",
                      "wall_ms": "milliseconds"},
        },
        "files": {"c1": "results_c1/risk_curves.parquet",
                  "c2ii": "results_c2/c2ii.parquet",
                  "c2iii": "results_c2/c2iii.parquet",
                  "c2iv": "results_c2/c2iv.parquet"},
    }
    (ROOT / "results_schema.yaml").write_text(json.dumps(schema, indent=1))

    total_h = sum(cell_seconds(c) for c in cells) / 3600
    print(f"wrote {len(manifest['notebooks'])} notebooks to colab/, "
          f"manifest + schema written; main-fleet serial projection "
          f"{total_h:.1f} h loaded-units")


if __name__ == "__main__":
    main()
