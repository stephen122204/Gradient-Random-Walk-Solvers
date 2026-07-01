#!/usr/bin/env python3
"""
Small launcher for running a JSON config without remembering commands.

Usage:
    python run_config.py
    python run_config.py configs/fhn_grw_steady.json --mode verify
    python run_config.py configs/heat_step_dirichlet.json --mode simulate

If no config path is provided, the script tries to open a file picker. If a GUI
is unavailable, it falls back to a terminal prompt.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/grw-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/grw-cache")

import config as config_module
from utils import plot_results
from simulation import (
    simulate_burgers,
    simulate_fitzhugh_nagumo,
    simulate_heat_equation,
)


def choose_config_path() -> str:
    """Choose a JSON config path with a file picker, falling back to input()."""
    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        path = filedialog.askopenfilename(
            title="Choose a GRW JSON config",
            initialdir=str(Path.cwd() / "configs"),
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        root.destroy()
        if path:
            return path
    except Exception:
        pass

    return input("Path to JSON config, e.g. configs/fhn_grw_steady.json: ").strip()


def choose_mode(default: str = "verify") -> str:
    """Ask for run mode when not provided on the command line."""
    answer = input(
        "Run mode [verify/simulate] "
        f"(default: {default}; verify makes comparison/metrics plots): "
    ).strip().lower()
    return answer or default


def infer_equation(cfg) -> str:
    eq = (cfg.equation_type or "").strip().lower()
    if eq in {"fitzhugh-nagumo", "fitzhugh–nagumo", "fitzhugh", "nagumo"}:
        return "fhn"
    if eq in {"burgers", "burger", "burgers'"}:
        return "burgers"
    if eq == "heat":
        return "heat"
    raise ValueError(f"Unknown equation type: {cfg.equation_type!r}")


def run_simulation(cfg) -> None:
    """Run the raw solver path and save the standard output plot."""
    eq = infer_equation(cfg)
    if eq == "heat":
        globs = [{"position": pos, "value": val} for pos, val in cfg.initial_conditions]
        results = simulate_heat_equation(globs, cfg)
    elif eq == "fhn":
        globs = [
            {"position": float(pos), "value": float(val)}
            for pos, val in cfg.initial_conditions
        ]
        results = simulate_fitzhugh_nagumo(globs, cfg)
    else:
        globs = [
            {"position": pos, "value": [float(val)]}
            for pos, val in cfg.initial_conditions
        ]
        results = simulate_burgers(globs, cfg)

    print("[run_config] Simulation complete. Saving standard plot...")
    plot_results(results, cfg.equation_type, cfg)


def run_verification(cfg, output_dir: str | None = None) -> None:
    """Run the verification/comparison path for the config's equation."""
    import verify_solver

    eq = infer_equation(cfg)
    if output_dir is None:
        output_dir = os.path.join("output", "verify", eq)

    if eq == "heat":
        verify_solver.run_heat(cfg, output_dir, do_save_data=False)
    elif eq == "fhn":
        verify_solver.run_fhn(cfg, output_dir, do_save_data=False, ref_factor=5)
    else:
        verify_solver.run_burgers(cfg, output_dir, do_save_data=False, ref_factor=5)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a GRW JSON config and generate plots.",
    )
    parser.add_argument(
        "config",
        nargs="?",
        help="Path to a JSON config. If omitted, a file picker or prompt is used.",
    )
    parser.add_argument(
        "--mode",
        choices=["verify", "simulate"],
        default=None,
        help=(
            "verify creates comparison plots and metrics; simulate creates the "
            "standard solver output plot."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Verification output directory. Only used with --mode verify.",
    )
    args = parser.parse_args()

    config_path = args.config or choose_config_path()
    if not config_path:
        raise SystemExit("No config path selected.")

    cfg = config_module.load_config_from_json(config_path)
    mode = args.mode or choose_mode(default="verify")
    if mode not in {"verify", "simulate"}:
        raise SystemExit("Mode must be 'verify' or 'simulate'.")

    print(f"[run_config] Config: {config_path}")
    print(f"[run_config] Equation: {cfg.equation_type}")
    print(f"[run_config] Mode: {mode}")

    if mode == "verify":
        run_verification(cfg, output_dir=args.output_dir)
    else:
        run_simulation(cfg)


if __name__ == "__main__":
    main()
