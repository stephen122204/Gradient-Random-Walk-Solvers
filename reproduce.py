#!/usr/bin/env python3
"""
reproduce.py — one command per paper artifact (Paper 1, GRW numerical study).

Wrapper-only entry point: this script calls the existing study and figure code
with the exact configurations used in the manuscript. It does not modify, and
must never modify, any solver or study code.

Targets
    python reproduce.py studies   # Figs 2, 5, 8 data + Tables 1 data (seed 42)
    python reproduce.py figures   # all 8 manuscript figures from archived arrays
    python reproduce.py all       # studies + figures
    python reproduce.py verify    # re-run studies, compare every reported value
                                  # against expected_values.json (PASS/FAIL)
    python reproduce.py verify --deep
                                  # additionally re-run the three representative
                                  # simulations (seed 42) and require bit-identical
                                  # arrays against figure_data/*.npz
    python reproduce.py figures --rerun-arrays
                                  # regenerate the archived arrays from fresh
                                  # seed-42 simulations before plotting
                                  # (overwrites figure_data/representative_figure_arrays.npz)

All paper configurations are pinned inside study_paper_refinement.py and
figure_scripts/regenerate_manuscript_figures.py (seed 42 throughout). The JSON
configs under configs/ and main.py are interactive exploration tools; they are
not the source of the paper's numbers.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STUDY_OUT = ROOT / "output" / "paper_refinement_original_grw"
DATA_DIR = ROOT / "figure_data"
EXPECTED = ROOT / "expected_values.json"

CSV_NAMES = (
    "heat_refinement_summary.csv",
    "fhn_refinement_summary.csv",
    "burgers_domain_sensitivity_summary.csv",
)
# Wall-clock timing is machine-specific and excluded from comparison.
SKIP_FIELDS = {"runtime_s"}


def run_studies() -> float:
    t0 = time.perf_counter()
    subprocess.run([sys.executable, str(ROOT / "study_paper_refinement.py")],
                   cwd=ROOT, check=True)
    return time.perf_counter() - t0


def run_figures(rerun_arrays: bool = False) -> float:
    cmd = [sys.executable, str(ROOT / "figure_scripts" / "regenerate_manuscript_figures.py")]
    if rerun_arrays:
        cmd.append("--rerun")
    t0 = time.perf_counter()
    subprocess.run(cmd, cwd=ROOT, check=True)
    return time.perf_counter() - t0


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def verify(deep: bool = False) -> int:
    expected = json.loads(EXPECTED.read_text())
    failures = 0
    checks = 0

    print("== 1/3  Re-running the three refinement/domain studies (seed 42) ...")
    elapsed = run_studies()
    print(f"        done in {elapsed:.1f}s")

    print("== 2/3  Comparing study outputs against expected_values.json")
    for name in CSV_NAMES:
        exp_rows = expected["study_csvs"][name]
        got_rows = _read_rows(STUDY_OUT / name)
        if len(exp_rows) != len(got_rows):
            print(f"  FAIL  {name}: {len(got_rows)} rows, expected {len(exp_rows)}")
            failures += 1
            continue
        bad = []
        for exp, got in zip(exp_rows, got_rows):
            for key, val in exp.items():
                if key in SKIP_FIELDS:
                    continue
                checks += 1
                if got.get(key) != val:
                    bad.append((exp.get("N", exp.get("L", "?")), key, val, got.get(key)))
        if bad:
            failures += len(bad)
            print(f"  FAIL  {name}: {len(bad)} field mismatches")
            for row_id, key, want, got_v in bad[:10]:
                print(f"          N/L={row_id} {key}: expected {want}, got {got_v}")
        else:
            print(f"  PASS  {name}: every reported field identical "
                  f"({len(exp_rows)} rows)")

    print("== 3/3  Recomputing representative (Table 2) metrics from archived arrays")
    sys.path.insert(0, str(ROOT / "figure_scripts"))
    import regenerate_manuscript_figures as figs  # noqa: E402  (wrapper import, read-only use)

    arrays = figs.generate_or_load_arrays(force=False)
    metrics = figs.compute_metrics(arrays)
    flat = {
        "heat.L2h": metrics["heat"]["L2h"],
        "heat.Linf": metrics["heat"]["Linf"],
        "heat.rel_L2": metrics["heat"]["rel_L2"],
        "fhn_t9.L2h": metrics["fhn"]["final"]["L2h"],
        "fhn_t9.Linf": metrics["fhn"]["final"]["Linf"],
        "fhn_t9.rel_L2": metrics["fhn"]["final"]["rel_L2"],
        "burgers.L2h": metrics["burgers"]["L2h"],
        "burgers.Linf": metrics["burgers"]["Linf"],
        "burgers.rel_L2": metrics["burgers"]["rel_L2"],
        "burgers.E_det": metrics["burgers"]["E_det"],
        "burgers.E_grw": metrics["burgers"]["E_grw"],
        "burgers.E_total": metrics["burgers"]["E_total"],
        "fhn_front.grw_T9": float(arrays["fhn_front_location"][-1]),
        "fhn_front.anchored_ref_T9": float(arrays["fhn_front_reference"][-1]),
    }
    for key, want in expected["representative_metrics"].items():
        checks += 1
        got = flat[key]
        if got == want:
            print(f"  PASS  {key} = {got!r}")
        else:
            failures += 1
            print(f"  FAIL  {key}: expected {want!r}, got {got!r}")

    if deep:
        import numpy as np

        print("== deep  Re-running representative simulations (seed 42), "
              "requiring bit-identical arrays")
        t0 = time.perf_counter()
        fresh: dict = {}
        fresh.update(figs.simulate_heat())
        fresh.update(figs.simulate_fhn())
        fresh.update(figs.simulate_burgers())
        elapsed = time.perf_counter() - t0
        archived = dict(np.load(DATA_DIR / "representative_figure_arrays.npz"))
        for key in sorted(archived):
            checks += 1
            if key in fresh and np.array_equal(fresh[key], archived[key]):
                continue
            failures += 1
            print(f"  FAIL  array {key} differs from archive")
        print(f"  {'PASS  all' if failures == 0 else 'checked'} "
              f"{len(archived)} archived arrays against fresh seed-42 run "
              f"({elapsed:.1f}s)")

    print()
    if failures == 0:
        print(f"VERIFY: PASS — {checks} checks, all identical to the "
              f"values reported in the paper.")
        return 0
    print(f"VERIFY: FAIL — {failures} of {checks} checks differ.")
    return 1


def main() -> int:
    args = sys.argv[1:]
    if not args or args[0] not in {"studies", "figures", "all", "verify"}:
        print(__doc__)
        return 2
    target = args[0]
    if target == "studies":
        elapsed = run_studies()
        print(f"\nstudies target complete in {elapsed:.1f}s "
              f"(outputs in {STUDY_OUT.relative_to(ROOT)}/)")
        return 0
    if target == "figures":
        elapsed = run_figures(rerun_arrays="--rerun-arrays" in args)
        print(f"\nfigures target complete in {elapsed:.1f}s")
        return 0
    if target == "all":
        t_studies = run_studies()
        t_figures = run_figures()
        print(f"\nall target complete (studies {t_studies:.1f}s, "
              f"figures {t_figures:.1f}s)")
        return 0
    return verify(deep="--deep" in args)


if __name__ == "__main__":
    raise SystemExit(main())
