#!/usr/bin/env python3
"""
reproduce.py — one command per artifact of the GRW accuracy paper
(arXiv:2608.22592).

Wrapper-only entry point: this script calls the existing study and figure code
with the exact configurations used in the paper. It does not modify, and
must never modify, any solver or study code.

Targets
    python reproduce.py studies   # single-seed representative study data
                                  # (seed 42; legacy diagnostics kept from the
                                  # original verification study)
    python reproduce.py figures   # representative-figure set from archived
                                  # arrays (legacy set; the combined paper's
                                  # ten figures come from the `paper` target)
    python reproduce.py all       # studies + figures (representative layer)
    python reproduce.py verify    # re-run studies, compare every reported value
                                  # against expected_values.json (PASS/FAIL)
    python reproduce.py verify --deep
                                  # additionally re-run the three representative
                                  # simulations (seed 42) and compare the arrays
                                  # against figure_data/*.npz
    python reproduce.py figures --rerun-arrays
                                  # regenerate the archived arrays from fresh
                                  # seed-42 simulations before plotting
                                  # (overwrites figure_data/representative_figure_arrays.npz)
    python reproduce.py ensembles # multi-seed ensemble, paired-grid, and
                                  # Cole-Hopf control studies (t4 t7 t5 t3 t8)
    python reproduce.py t3|t4|t5|t7|t8
                                  # one ensemble study (see verify_ensembles.py)
    python reproduce.py verify-ensembles
                                  # re-run all five ensemble studies and compare
                                  # every pinned numeric field against
                                  # pinned_ensembles/ (PASS/FAIL)
    python reproduce.py verify-all
                                  # release gate: verify --deep, then re-run and
                                  # compare all five ensemble studies. PASS only
                                  # if both groups pass.
    python reproduce.py paper     # regenerate the ten combined-paper figures
    python reproduce.py paper1-figures
                                  # compatibility alias for ``paper``; writes to
                                  # output/final_prepublication_tests/paper_figures/
                                  # with a SHA-256 provenance manifest

Standalone checks that are not targets of this script live in checks/
(see README, Additional Checks).

All paper configurations are pinned inside study_paper_refinement.py,
figure_scripts/regenerate_paper_figures.py (seed 42 throughout), and the
studies/ scripts (fixed seed lists documented in each study). The JSON
configs under configs/ and main.py are interactive exploration tools; they
are not the source of the paper's numbers.

Numerical tolerances
    Float comparisons use rel/abs tolerance 1e-12. Under the pinned
    environment in requirements.txt the regeneration is bit-identical; under
    other NumPy/BLAS builds, last-digit floating-point noise up to about
    2.4e-14 relative has been observed. The 1e-12 tolerance sits well above
    that platform noise and many orders of magnitude below the precision at
    which any value is reported in the paper, so a genuine change in any
    paper value cannot pass. Integers, row counts, seeds, and identifiers
    are always compared exactly.
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

# Float tolerance: far above cross-platform noise (~2.4e-14 relative observed
# across NumPy builds), far below paper-value precision. See module docstring.
REL_TOL = 1e-12
ABS_TOL = 1e-12


def _float_close(got: float, want: float) -> bool:
    return abs(got - want) <= max(ABS_TOL, REL_TOL * max(abs(got), abs(want)))


def _compare_field(want: str, got: str) -> tuple[bool, str]:
    """Typed comparison of one CSV field. Ints exact, floats within tolerance,
    everything else string-exact. Returns (ok, detail)."""
    try:
        want_i, got_i = int(want), int(got)
        return want_i == got_i, f"int expected {want_i}, got {got_i}"
    except ValueError:
        pass
    try:
        want_f, got_f = float(want), float(got)
        diff = abs(got_f - want_f)
        tol = max(ABS_TOL, REL_TOL * max(abs(got_f), abs(want_f)))
        return diff <= tol, (f"expected {want_f!r}, got {got_f!r}, "
                             f"|diff|={diff:.3e}, tol={tol:.3e}")
    except ValueError:
        return want == got, f"expected {want!r}, got {got!r} (exact string)"


def run_studies() -> float:
    t0 = time.perf_counter()
    subprocess.run([sys.executable, str(ROOT / "study_paper_refinement.py")],
                   cwd=ROOT, check=True)
    return time.perf_counter() - t0


def run_figures(rerun_arrays: bool = False) -> float:
    cmd = [sys.executable, str(ROOT / "figure_scripts" / "regenerate_paper_figures.py")]
    if rerun_arrays:
        cmd.append("--rerun")
    t0 = time.perf_counter()
    subprocess.run(cmd, cwd=ROOT, check=True)
    return time.perf_counter() - t0


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def verify(deep: bool = False) -> int:
    import math

    import numpy as np

    expected = json.loads(EXPECTED.read_text())
    failures = 0
    checks = 0

    print("== 1/4  Re-running the three refinement/domain studies (seed 42) ...")
    elapsed = run_studies()
    print(f"        done in {elapsed:.1f}s")

    print("== 2/4  Comparing study outputs against expected_values.json")
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
                ok, detail = _compare_field(val, got.get(key, ""))
                if not ok:
                    bad.append((exp.get("N", exp.get("L", "?")), key, detail))
        if bad:
            failures += len(bad)
            print(f"  FAIL  {name}: {len(bad)} field mismatches beyond tolerance")
            for row_id, key, detail in bad[:10]:
                print(f"          N/L={row_id} {key}: {detail}")
        else:
            print(f"  PASS  {name}: every reported field matches within "
                  f"tolerance ({len(exp_rows)} rows)")

    print("== 3/4  Recomputing representative metrics from archived arrays")
    sys.path.insert(0, str(ROOT / "figure_scripts"))
    import regenerate_paper_figures as figs  # noqa: E402  (wrapper import, read-only use)

    arrays = figs.generate_or_load_arrays(force=False)
    metrics = figs.compute_metrics(arrays)
    flat = {
        "heat.L2h": metrics["heat"]["L2h"],
        "heat.Linf": metrics["heat"]["Linf"],
        "heat.rel_L2": metrics["heat"]["rel_L2"],
        "heat.L2h_nbins300": metrics["heat"]["L2h_nbins300"],
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
        if _float_close(got, want):
            print(f"  PASS  {key} = {got!r}")
        else:
            failures += 1
            print(f"  FAIL  {key}: expected {want!r}, got {got!r}, "
                  f"|diff|={abs(got - want):.3e}, tol={max(ABS_TOL, REL_TOL * max(abs(got), abs(want))):.3e}")

    print("== 4/4  Checking derived paper values (fitted speed, E_det*sqrt(L), "
          "kernel factor, negative weights, floor)")
    derived = expected["derived_values"]

    def check(name: str, got, want, exact: bool = False) -> None:
        nonlocal checks, failures
        checks += 1
        ok = (got == want) if exact else _float_close(float(got), float(want))
        if ok:
            print(f"  PASS  {name} = {got!r}")
        else:
            failures += 1
            print(f"  FAIL  {name}: expected {want!r}, got {got!r}")

    fitted = float(np.polyfit(arrays["fhn_front_time"], arrays["fhn_front_location"], 1)[0])
    check("fhn fitted front speed (paper: -0.345)", fitted, derived["fhn_fitted_front_speed"])
    exact_speed = -math.sqrt(2.0) * 0.25
    rel_err = abs(fitted - exact_speed) / abs(exact_speed)
    check("fhn fitted-speed relative error <= 2.5% (paper claim)",
          rel_err <= derived["fhn_fitted_speed_rel_err_vs_exact_max"], True, exact=True)

    dom_rows = _read_rows(DATA_DIR / "burgers_domain_sensitivity_summary.csv")
    prods = [float(r["bc_mismatch_RMSE"]) * math.sqrt(float(r["L"])) for r in dom_rows]
    for prod, want in zip(prods, derived["e_det_sqrtL_products"]):
        check("E_det*sqrt(L) product (paper Sec 6.3)", prod, want)

    check("kernel variance-effective bins (paper: ~43)",
          2 * 12 * math.sqrt(math.pi), derived["kernel_variance_effective_bins"])
    check("kernel noise-reduction factor (paper: ~6.5)",
          math.sqrt(2 * 12 * math.sqrt(math.pi)), derived["kernel_noise_reduction_factor"])

    w9 = arrays["fhn_weights_t9"]
    p9 = arrays["fhn_positions_t9"]
    neg = np.where(w9 < 0)[0]
    check("fhn negative-weight count at t=9 (paper: two globs near x~26)",
          int(len(neg)), derived["fhn_negative_weight_count_t9"], exact=True)
    for pos, want in zip(sorted(float(p) for p in p9[neg]),
                         derived["fhn_negative_weight_positions_t9"]):
        check("fhn negative-weight position", pos, want)

    check("burgers representative clipped bins (paper: floor never active)",
          int(arrays["burgers_clipped"].sum()),
          derived["burgers_clipped_bins_representative"], exact=True)

    if deep:
        print("== deep  Re-running representative simulations (seed 42) and "
              "comparing arrays against the archive")
        t0 = time.perf_counter()
        fresh: dict = {}
        fresh.update(figs.simulate_heat())
        fresh.update(figs.simulate_fhn())
        fresh.update(figs.simulate_burgers())
        elapsed = time.perf_counter() - t0
        archived = dict(np.load(DATA_DIR / "representative_figure_arrays.npz"))
        bit_identical = 0
        for key in sorted(archived):
            checks += 1
            if key not in fresh:
                failures += 1
                print(f"  FAIL  array {key} missing from fresh run")
                continue
            if np.array_equal(fresh[key], archived[key]):
                bit_identical += 1
                continue
            a = archived[key].astype(float)
            f = fresh[key].astype(float)
            if a.shape == f.shape and np.allclose(f, a, rtol=REL_TOL, atol=ABS_TOL):
                max_dev = float(np.max(np.abs(a - f)))
                print(f"  PASS  array {key}: within tolerance "
                      f"(max |diff|={max_dev:.3e}, platform floating-point noise)")
                continue
            failures += 1
            max_dev = float(np.max(np.abs(a - f))) if a.shape == f.shape else float("nan")
            print(f"  FAIL  array {key}: differs beyond tolerance "
                  f"(max |diff|={max_dev:.3e}, tol rel/abs {REL_TOL:g}/{ABS_TOL:g})")
        print(f"  checked {len(archived)} archived arrays against a fresh "
              f"seed-42 run ({elapsed:.1f}s): {bit_identical} bit-identical, "
              f"{len(archived) - bit_identical} within tolerance on this "
              f"platform (all {len(archived)} are bit-identical under the "
              f"pinned environment in requirements.txt)")

    print()
    if failures == 0:
        print(f"VERIFY: PASS — {checks} checks, all consistent with the "
              f"values reported in the paper within the documented "
              f"tolerances (rel/abs {REL_TOL:g}).")
        return 0
    print(f"VERIFY: FAIL — {failures} of {checks} checks differ beyond tolerance.")
    return 1


def main() -> int:
    args = sys.argv[1:]
    known = {"studies", "figures", "all", "verify", "ensembles",
             "verify-ensembles", "verify-all", "paper", "paper1-figures",
             "t3", "t4", "t5", "t7", "t8"}
    if not args or args[0] not in known:
        print(__doc__)
        return 2
    target = args[0]
    if target == "verify-all":
        deep_rc = verify(deep=True)
        import verify_ensembles
        ens_rc = verify_ensembles.verify(rerun=True)
        ok = deep_rc == 0 and ens_rc == 0
        print(f"\nVERIFY-ALL: {'PASS' if ok else 'FAIL'} — deep representative "
              f"checks {'passed' if deep_rc == 0 else 'FAILED'}, ensemble "
              f"comparisons {'passed' if ens_rc == 0 else 'FAILED'}.")
        return 0 if ok else 1
    if target in {"t3", "t4", "t5", "t7", "t8"}:
        import verify_ensembles
        verify_ensembles.run_study(target)
        return 0
    if target == "ensembles":
        import verify_ensembles
        for name in verify_ensembles.STUDIES:
            verify_ensembles.run_study(name)
        return 0
    if target == "verify-ensembles":
        import verify_ensembles
        return verify_ensembles.verify(rerun="--no-rerun" not in args)
    if target in {"paper", "paper1-figures"}:
        import hashlib
        import shutil
        import tempfile

        ens_out = ROOT / "output" / "final_prepublication_tests"
        paper_dir = ens_out / "paper_figures"
        if paper_dir.exists():
            shutil.rmtree(paper_dir)
        paper_dir.mkdir(parents=True, exist_ok=True)

        # The representative generator normally writes to a timestamped
        # exploration directory. Give it a private temporary destination here
        # so this target never guesses which prior run is the newest.
        with tempfile.TemporaryDirectory(prefix="paper1_figures_") as tmp:
            representative_dir = Path(tmp)
            subprocess.run(
                [sys.executable,
                 str(ROOT / "figure_scripts" / "regenerate_paper_figures.py"),
                 "--output-dir", str(representative_dir)],
                cwd=ROOT,
                check=True,
            )
            for name in ("heat_comparison.pdf", "fhn_comparison.pdf",
                         "fhn_diagnostics.pdf", "burgers_diagnostics.pdf"):
                shutil.copy2(representative_dir / name, paper_dir / name)

        # The six ensemble/control figures are always redrawn from committed
        # pinned data. A clean checkout therefore needs no prior study output.
        subprocess.run([sys.executable,
                        str(ROOT / "figure_scripts" / "regenerate_ensemble_figures.py")],
                       cwd=ROOT, check=True)

        wanted = {
            "heat_comparison.pdf": "figure_data/representative_figure_arrays.npz (seed 42)",
            "fhn_comparison.pdf": "figure_data/representative_figure_arrays.npz (seed 42)",
            "fhn_diagnostics.pdf": "figure_data/representative_figure_arrays.npz (seed 42)",
            "burgers_diagnostics.pdf": "figure_data/representative_figure_arrays.npz (seed 42)",
            "heat_bias_spread_total_vs_N.pdf": "pinned_ensembles/heat_extended/summary_by_N.csv",
            "fhn_convergence.pdf": "pinned_ensembles/fhn_extended/summary_by_N.csv",
            "heat_grid_paired.pdf": "pinned_ensembles/heat_grid_paired/summary.json",
            "burgers_decoupled.pdf": "pinned_ensembles/burgers_controls/summary.json",
            "burgers_boundary_domain.pdf": "pinned_ensembles/burgers_controls/summary.json",
            "burgers_perturbation_response.pdf": "pinned_ensembles/burgers_controls/summary.json",
        }
        manifest = {}
        for name, source in wanted.items():
            path = paper_dir / name
            if path.exists():
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                manifest[name] = {
                    "source": source,
                    "sha256": digest,
                    "bytes": path.stat().st_size,
                }
            else:
                manifest[name] = {"source": source, "sha256": None,
                                  "error": "not generated"}
        with (paper_dir / "paper1_figure_manifest.json").open("w") as stream:
            json.dump(manifest, stream, indent=2)
        missing = [k for k, v in manifest.items() if v["sha256"] is None]
        print(f"\npaper: {len(manifest) - len(missing)}/{len(manifest)} "
              f"canonical figures in {paper_dir.relative_to(ROOT)}/ "
              f"(manifest: paper1_figure_manifest.json)")
        return 1 if missing else 0
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
