#!/usr/bin/env python3
"""Rerun the five original ensemble/control studies and compare with committed data."""
from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output" / "final_prepublication_tests"
PIN = ROOT / "pinned_ensembles"

STUDIES = {
    "t4": ("studies/study_t4_heat_extended.py", "heat_extended"),
    "t7": ("studies/study_t7_heat_grid_paired.py", "heat_grid_paired"),
    "t5": ("studies/study_t5_fhn_extended.py", "fhn_extended"),
    "t3": ("studies/study_t3_cole_hopf_plateau.py", "cole_hopf_plateau"),
    "t8": ("studies/study_t8_burgers_controls.py", "burgers_controls"),
}

REL_TOL = 1e-9
ABS_TOL = 1e-12
SKIP_TOKENS = ("runtime", "_rt_", "rt_s", "mean_rt", "elapsed", "time_s",
               "identity_residual")


def _skip(name: str) -> bool:
    low = name.lower()
    return any(t in low for t in SKIP_TOKENS)


def _close(a: float, b: float) -> bool:
    if math.isnan(a) and math.isnan(b):
        return True
    return abs(a - b) <= max(ABS_TOL, REL_TOL * max(abs(a), abs(b)))


def run_study(name: str) -> None:
    script, _ = STUDIES[name]
    print(f"== running {name}: {script}")
    subprocess.run([sys.executable, script], cwd=ROOT, check=True)


def _compare_csv(pinned: Path, fresh: Path, stats: dict) -> None:
    rp = list(csv.reader(open(pinned)))
    rf = list(csv.reader(open(fresh)))
    rel = pinned.relative_to(PIN)
    if len(rp) != len(rf):
        stats["fail"] += 1
        print(f"  FAIL  {rel}: row count {len(rp)} vs {len(rf)}")
        return
    hdr = rp[0]
    bad = 0
    for i, (a, b) in enumerate(zip(rp[1:], rf[1:]), 1):
        for j, (x, y) in enumerate(zip(a, b)):
            if j < len(hdr) and _skip(hdr[j]):
                continue
            try:
                if not _close(float(x), float(y)):
                    bad += 1
                    if bad <= 3:
                        print(f"    {rel} row {i} col {hdr[j]}: {x} vs {y}")
            except ValueError:
                if x != y:
                    bad += 1
    stats["checks"] += 1
    if bad:
        stats["fail"] += 1
        print(f"  FAIL  {rel}: {bad} differing fields")
    else:
        print(f"  PASS  {rel} ({len(rp)-1} rows)")


def _walk(a, b, path: str, diffs: list) -> None:
    if isinstance(a, dict):
        for k in a:
            if _skip(str(k)):
                continue
            if not isinstance(b, dict) or k not in b:
                diffs.append(f"{path}/{k}: missing")
            else:
                _walk(a[k], b[k], f"{path}/{k}", diffs)
    elif isinstance(a, list):
        if not isinstance(b, list) or len(a) != len(b):
            diffs.append(f"{path}: list shape differs")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            _walk(x, y, f"{path}[{i}]", diffs)
    elif isinstance(a, float) or isinstance(b, float):
        try:
            if not _close(float(a), float(b)):
                diffs.append(f"{path}: {a} vs {b}")
        except (TypeError, ValueError):
            diffs.append(f"{path}: {a!r} vs {b!r}")
    elif a != b:
        diffs.append(f"{path}: {a!r} vs {b!r}")


def _compare_json(pinned: Path, fresh: Path, stats: dict) -> None:
    rel = pinned.relative_to(PIN)
    diffs: list = []
    _walk(json.load(open(pinned)), json.load(open(fresh)), "", diffs)
    stats["checks"] += 1
    if diffs:
        stats["fail"] += 1
        print(f"  FAIL  {rel}: {len(diffs)} differing fields")
        for d in diffs[:5]:
            print(f"    {d}")
    else:
        print(f"  PASS  {rel}")


def verify(rerun: bool = True) -> int:
    if rerun:
        for name in STUDIES:
            run_study(name)
    print("\n== comparing outputs against pinned_ensembles/")
    stats = {"checks": 0, "fail": 0}
    for name, (_, outdir) in STUDIES.items():
        for pinned in sorted((PIN / outdir).glob("*")):
            fresh = OUT / outdir / pinned.name
            if not fresh.exists():
                stats["checks"] += 1
                stats["fail"] += 1
                print(f"  FAIL  {pinned.relative_to(PIN)}: fresh output missing")
                continue
            if pinned.suffix == ".csv":
                _compare_csv(pinned, fresh, stats)
            elif pinned.suffix == ".json":
                _compare_json(pinned, fresh, stats)
    ok = stats["fail"] == 0
    print(f"\nensemble verification: {stats['checks'] - stats['fail']}/{stats['checks']} "
          f"file comparisons passed -> {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main() -> int:
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        return 2
    if args[0] in STUDIES:
        run_study(args[0])
        return 0
    if args[0] == "verify":
        return verify(rerun="--no-rerun" not in args)
    print(__doc__)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
