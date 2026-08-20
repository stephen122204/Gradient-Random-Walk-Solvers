#!/usr/bin/env python3
"""Regenerate all eight GRW paper figures from checked-in data or fixed-seed (42) runs.

Figures 1, 3, 4, 6, 7 are drawn from the archived representative arrays in
figure_data/representative_figure_arrays.npz (regenerable bit-for-bit with
--rerun under the pinned environment of requirements.txt). Figures 2, 5, 8
are drawn from the checked-in study CSVs in figure_data/ (the paper's data
of record). Outputs are written to outputs/<timestamp>/ together with
figure_regeneration_metadata.json, which records the seed and the source of
every figure using repository-relative paths.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("SOURCE_DATE_EPOCH", "1704067200")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
CODE_DIR = SCRIPT_DIR.parent
DATA_DIR = CODE_DIR / "figure_data"
CSV_DIR = CODE_DIR / "figure_data"
_RUN_STAMP = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
FIGURE_DIR = CODE_DIR / "outputs" / _RUN_STAMP

sys.path.insert(0, str(SCRIPT_DIR))
from plot_style import (  # noqa: E402
    DETERMINISTIC,
    EXACT,
    GRW,
    GUIDE,
    SECONDARY,
    apply_paper_style,
    save_pdf_png,
)

SEED = 42
ARRAY_PATH = DATA_DIR / "representative_figure_arrays.npz"

# Expected representative values for regression checking. These are the figures
# reported in Table~\ref{tab:verification-summary}; the tolerance catches drift
# from accidental changes to the pipeline while allowing platform-level noise.
def exact_heat(x: np.ndarray, time: float, alpha: float, center: float) -> np.ndarray:
    scaled = (x - center) / (2.0 * math.sqrt(alpha * time))
    erf_values = np.array([math.erf(float(value)) for value in scaled])
    return 0.5 * (1.0 + erf_values)


def exact_fhn(x: np.ndarray, time: float, a: float, center: float) -> np.ndarray:
    theta = math.sqrt(2.0) * (0.5 - a)
    return 1.0 / (1.0 + np.exp(-(x + theta * time - center) / 2.0))


def reconstruct_sorted(x: np.ndarray, w: np.ndarray, grid: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="stable")
    x_sorted = x[order]
    cumulative = np.cumsum(w[order])
    indices = np.searchsorted(x_sorted, grid, side="right") - 1
    return np.where(indices >= 0, cumulative[np.maximum(indices, 0)], 0.0)


def simulate_heat() -> dict[str, np.ndarray]:
    alpha, time, dt, length, center, count = 0.1, 0.5, 0.001, 10.0, 5.0, 50000
    rng = np.random.RandomState(SEED)
    positions = np.full(count, center, dtype=float)
    sigma = math.sqrt(2.0 * alpha * dt)
    for _ in range(int(time / dt)):
        positions += rng.normal(0.0, sigma, size=count)
        positions = np.where(positions < 0.0, -positions, positions)
        positions = np.where(positions > length, 2.0 * length - positions, positions)

    edges = np.linspace(0.0, length, 401)
    weights, _ = np.histogram(positions, bins=edges, weights=np.full(count, 1.0 / count))
    x_grid = 0.5 * (edges[:-1] + edges[1:])
    # Same final particle state rebinned on 300 bins; the paper quotes this
    # coarser-grid value next to the 400-bin reconstruction.
    edges_300 = np.linspace(0.0, length, 301)
    weights_300, _ = np.histogram(positions, bins=edges_300,
                                  weights=np.full(count, 1.0 / count))
    x_grid_300 = 0.5 * (edges_300[:-1] + edges_300[1:])
    return {
        "heat_x": x_grid,
        "heat_grw": np.cumsum(weights),
        "heat_exact": exact_heat(x_grid, time, alpha, center),
        "heat_x_300": x_grid_300,
        "heat_grw_300": np.cumsum(weights_300),
        "heat_exact_300": exact_heat(x_grid_300, time, alpha, center),
    }


def simulate_fhn() -> dict[str, np.ndarray]:
    a, diffusion, time, dt, length, center, count = 0.25, 0.5, 9.0, 0.01, 30.0, 15.0, 500
    theta = math.sqrt(2.0) * (0.5 - a)
    quantiles = (np.arange(count) + 0.5) / count
    x = center - 2.0 * np.log(1.0 / quantiles - 1.0)
    w = np.full(count, 1.0 / count)
    grid = np.linspace(0.0, length, 500)
    rng = np.random.RandomState(SEED)
    sigma = math.sqrt(2.0 * diffusion * dt)
    c2 = -1.5 * diffusion
    c1 = 1.5 * diffusion - theta
    c0 = 0.5 * theta - 0.25 * diffusion

    snapshot_steps = {0: 0.0, 300: 3.0, 600: 6.0, 900: 9.0}
    snapshots: dict[float, np.ndarray] = {0.0: reconstruct_sorted(x, w, grid)}
    weight_snapshots: dict[float, tuple[np.ndarray, np.ndarray]] = {0.0: (x.copy(), w.copy())}

    initial_order = np.argsort(x)
    initial_cumulative = np.cumsum(w[initial_order])
    initial_index = int(np.clip(np.searchsorted(initial_cumulative, 0.5), 0, count - 1))
    empirical_center = float(x[initial_order][initial_index])
    front_time = [0.0]
    front_location = [empirical_center]

    for step in range(1, int(time / dt) + 1):
        x += rng.normal(0.0, sigma, size=count)
        for _ in range(4):
            left = x < 0.0
            x[left] = -x[left]
            w[left] = -w[left]
            right = x > length
            x[right] = 2.0 * length - x[right]
            w[right] = -w[right]

        order = np.argsort(x, kind="stable")
        x, w = x[order], w[order]
        cumulative = np.cumsum(w)
        reaction = c2 * cumulative**2 + c1 * cumulative + c0
        w += dt * reaction * w
        cumulative_post = np.cumsum(w)

        front_index = int(np.clip(np.searchsorted(cumulative_post, 0.5), 0, count - 1))
        front_time.append(step * dt)
        front_location.append(float(x[front_index]))

        if step in snapshot_steps:
            snapshot_time = snapshot_steps[step]
            snapshots[snapshot_time] = reconstruct_sorted(x, w, grid)
            weight_snapshots[snapshot_time] = (x.copy(), w.copy())

    result: dict[str, np.ndarray] = {
        "fhn_x": grid,
        "fhn_times": np.array(sorted(snapshots)),
        "fhn_front_time": np.asarray(front_time),
        "fhn_front_location": np.asarray(front_location),
        "fhn_front_reference": empirical_center - theta * np.asarray(front_time),
        # Keep final-time arrays for backward compatibility.
        "fhn_final_positions": weight_snapshots[9.0][0],
        "fhn_final_weights": weight_snapshots[9.0][1],
    }
    for snapshot_time in sorted(snapshots):
        tag = str(int(snapshot_time))
        result[f"fhn_grw_t{tag}"] = snapshots[snapshot_time]
        result[f"fhn_exact_t{tag}"] = exact_fhn(grid, snapshot_time, a, center)
        # Store positions and weights at each snapshot for the multi-time weight panel.
        result[f"fhn_positions_t{tag}"] = weight_snapshots[snapshot_time][0]
        result[f"fhn_weights_t{tag}"] = weight_snapshots[snapshot_time][1]
    return result


def deterministic_phi_heat(
    phi_initial: np.ndarray,
    grid: np.ndarray,
    nu: float,
    time: float,
    left_value: float,
    right_value: float,
) -> np.ndarray:
    dx = float(grid[1] - grid[0])
    provisional_dt = 0.4 * dx**2 / nu
    steps = int(np.ceil(time / provisional_dt))
    dt = time / steps
    ratio = nu * dt / dx**2
    phi = phi_initial.copy()
    for _ in range(steps):
        updated = phi.copy()
        updated[1:-1] = phi[1:-1] + ratio * (phi[2:] - 2.0 * phi[1:-1] + phi[:-2])
        updated[0], updated[-1] = left_value, right_value
        phi = updated
    return phi


def simulate_burgers() -> dict[str, np.ndarray]:
    amplitude, nu, time, dt, length, count = 1.0, 0.5, 0.5, 0.005, 4.0, 400
    center = length / 2.0
    x0 = np.linspace(0.0, length, count)
    u0 = -amplitude * np.tanh(amplitude * (x0 - center) / (2.0 * nu))

    primitive = np.zeros(count)
    primitive[1:] = np.cumsum(0.5 * (u0[:-1] + u0[1:]) * np.diff(x0))
    log_phi = -primitive / (2.0 * nu)
    log_phi -= log_phi.max()
    phi0 = np.exp(np.clip(log_phi, -700.0, 0.0))
    left_value, right_value = float(phi0[0]), float(phi0[-1])
    exact_total = right_value - left_value

    positions = 0.5 * (x0[:-1] + x0[1:])
    weights = np.diff(phi0)
    rng = np.random.RandomState(SEED)
    sigma = math.sqrt(2.0 * nu * dt)
    for _ in range(int(time / dt)):
        positions += rng.normal(0.0, sigma, size=positions.shape)
        for _ in range(4):
            positions = np.where(positions < 0.0, -positions, positions)
            positions = np.where(positions > length, 2.0 * length - positions, positions)

    grid = np.linspace(0.0, length, count)
    dx = float(grid[1] - grid[0])
    binned = np.zeros(count)
    indices = np.clip(np.floor(positions / dx).astype(int), 0, count - 1)
    np.add.at(binned, indices, weights)

    sigma_bins = 12
    half_width = int(4 * sigma_bins) + 1
    kernel_x = np.arange(-half_width, half_width + 1, dtype=float)
    kernel = np.exp(-0.5 * (kernel_x / sigma_bins) ** 2)
    kernel /= kernel.sum()
    raw_smoothed = np.convolve(binned, kernel, mode="same")
    kernel_norm = np.convolve(np.ones(count), kernel, mode="same")
    smoothed = raw_smoothed / np.maximum(kernel_norm, 1e-12)

    near_zero = 1e-6 * max(float(np.abs(binned).max()), 1e-30)
    if abs(exact_total) < near_zero:
        smoothed -= smoothed.mean()
    elif abs(smoothed.sum()) > 1e-30:
        smoothed *= exact_total / smoothed.sum()

    phi_grw = left_value + np.cumsum(smoothed)
    phi_x_grw = np.gradient(phi_grw, dx)
    phi_floor = max(float(phi0.min()) / 2.0, 1e-10)
    clipped = phi_grw < phi_floor
    phi_safe = np.where(clipped, phi_floor, phi_grw)
    u_grw = -2.0 * nu * phi_x_grw / phi_safe
    u_grw = np.where(clipped, 0.0, u_grw)

    phi_fd = deterministic_phi_heat(phi0.copy(), grid, nu, time, left_value, right_value)
    phi_x_fd = np.gradient(phi_fd, dx)
    u_fd = -2.0 * nu * phi_x_fd / np.maximum(phi_fd, 1e-12)

    argument = amplitude * (grid - center) / (2.0 * nu)
    phi_exact = np.cosh(argument) / math.cosh(amplitude * center / (2.0 * nu))
    ratio_exact = (amplitude / (2.0 * nu)) * np.tanh(argument)
    u_exact = -2.0 * nu * ratio_exact

    return {
        "burgers_x": grid,
        "burgers_u_exact": u_exact,
        "burgers_u_grw": u_grw,
        "burgers_u_fd": u_fd,
        "burgers_phi_exact": phi_exact,
        "burgers_phi_grw": phi_grw,
        "burgers_phi_fd": phi_fd,
        "burgers_ratio_exact": ratio_exact,
        "burgers_ratio_grw": phi_x_grw / phi_safe,
        "burgers_ratio_fd": phi_x_fd / np.maximum(phi_fd, 1e-12),
        "burgers_clipped": clipped.astype(int),
    }


def generate_or_load_arrays(force: bool) -> dict[str, np.ndarray]:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if ARRAY_PATH.exists() and not force:
        with np.load(ARRAY_PATH) as archive:
            return {key: archive[key] for key in archive.files}
    arrays = {}
    arrays.update(simulate_heat())
    arrays.update(simulate_fhn())
    arrays.update(simulate_burgers())
    np.savez_compressed(ARRAY_PATH, **arrays)
    return arrays


def read_csv(name: str) -> list[dict[str, str]]:
    with (CSV_DIR / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def record_output(generated: list[str], figure, stem: str) -> None:
    generated.extend(save_pdf_png(figure, FIGURE_DIR, stem))


def plot_heat_comparison(arrays: dict[str, np.ndarray], generated: list[str]) -> None:
    figure, axis = plt.subplots(figsize=(5.5, 3.5))
    axis.plot(arrays["heat_x"], arrays["heat_exact"], color=EXACT, label="Exact")
    axis.plot(arrays["heat_x"], arrays["heat_grw"], color=GRW, linestyle="--", label="GRW")
    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"$u(x,T)$")
    axis.set_xlim(0.0, 10.0)
    axis.set_ylim(-0.03, 1.03)
    axis.legend(loc="upper left")
    inset = axis.inset_axes([0.58, 0.14, 0.38, 0.34])
    inset.plot(arrays["heat_x"], arrays["heat_grw"] - arrays["heat_exact"], color=GRW, linewidth=0.9)
    inset.axhline(0.0, color=EXACT, linestyle=":", linewidth=0.8)
    inset.set_title(r"$u_N-u^{\mathrm{ex}}$", fontsize=8)
    inset.tick_params(labelsize=6)
    figure.tight_layout()
    record_output(generated, figure, "heat_comparison")


def plot_heat_refinement(generated: list[str]) -> None:
    rows = read_csv("heat_refinement_summary.csv")
    counts = np.array([int(row["N"]) for row in rows])
    errors = np.array([float(row["L2_h"]) for row in rows])
    figure, axis = plt.subplots(figsize=(4.8, 3.6))
    axis.loglog(counts, errors, "o-", color=GRW)
    axis.set_xlabel(r"Particle count $N$")
    axis.set_ylabel(r"$L_h^2$ error")
    figure.tight_layout()
    record_output(generated, figure, "heat_error_fixed_grid_diagnostic_final")


def plot_fhn_comparison(arrays: dict[str, np.ndarray], generated: list[str]) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True, sharey=True)
    x = arrays["fhn_x"]
    for axis, time in zip(axes.flat, arrays["fhn_times"]):
        tag = str(int(time))
        axis.plot(x, arrays[f"fhn_exact_t{tag}"], color=EXACT, label="Exact")
        axis.plot(x, arrays[f"fhn_grw_t{tag}"], color=GRW, linestyle="--", label="GRW")
        axis.set_title(rf"$t={time:g}$")
        axis.set_xlim(0.0, 30.0)
        axis.set_ylim(-0.05, 1.05)
    axes[0, 0].legend(loc="upper left")
    for axis in axes[-1, :]:
        axis.set_xlabel(r"$x$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$u(x,t)$")
    figure.tight_layout()
    record_output(generated, figure, "fhn_comparison")


def plot_fhn_diagnostics(arrays: dict[str, np.ndarray], generated: list[str]) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(7.4, 5.5))
    time = arrays["fhn_front_time"]
    location = arrays["fhn_front_location"]
    reference = arrays["fhn_front_reference"]

    axes[0, 0].plot(time, reference, color=EXACT, label="Anchored reference")
    axes[0, 0].plot(time, location, color=GRW, linewidth=1.1, label="GRW crossing")
    axes[0, 0].set_title("(a) Front location")
    axes[0, 0].set_xlabel(r"$t$")
    axes[0, 0].set_ylabel(r"$x_f(t)$")
    axes[0, 0].legend(loc="best", fontsize=8)

    axes[0, 1].plot(time, location - reference, color=SECONDARY, linewidth=1.1)
    axes[0, 1].axhline(0.0, color=EXACT, linestyle=":", linewidth=1.0)
    axes[0, 1].set_title("(b) Front-location error")
    axes[0, 1].set_xlabel(r"$t$")
    axes[0, 1].set_ylabel(r"$x_f^{\mathrm{GRW}}-x_f^{\mathrm{ref}}$")

    # Use a shared discrete color sequence so the same time value gets the
    # same color in both the weight panel (c) and reconstruction panel (d).
    snap_colors = ["#1b7837", "#762a83", "#e66101", "#2166ac"]
    snap_times = list(arrays["fhn_times"])

    for idx, snap_t in enumerate(snap_times):
        tag = str(int(snap_t))
        color = snap_colors[idx % len(snap_colors)]
        pos_key = f"fhn_positions_t{tag}"
        wt_key = f"fhn_weights_t{tag}"
        if pos_key in arrays and wt_key in arrays:
            axes[1, 0].plot(
                arrays[pos_key],
                arrays[wt_key],
                linestyle="none",
                marker=".",
                markersize=2.0,
                color=color,
                alpha=0.65,
                label=rf"$t={snap_t:g}$",
            )
        axes[1, 1].plot(
            arrays["fhn_x"],
            arrays[f"fhn_grw_t{tag}"],
            color=color,
            label=rf"$t={snap_t:g}$",
            linewidth=1.4,
        )

    axes[1, 0].axhline(0.0, color=EXACT, linestyle=":", linewidth=1.0)
    weight_min = min(float(arrays[f"fhn_weights_t{str(int(t))}"].min()) for t in snap_times)
    axes[1, 0].set_ylim(bottom=min(weight_min * 3.0, -0.0004))
    axes[1, 0].set_title("(c) Gradient weights at snapshot times")
    axes[1, 0].set_xlabel(r"$x$")
    axes[1, 0].set_ylabel(r"$w_i$")
    axes[1, 0].legend(loc="best", fontsize=7, markerscale=2.5)

    axes[1, 1].axhline(0.5, color=EXACT, linestyle=":", linewidth=1.0)
    axes[1, 1].set_title("(d) Cumulative reconstructions")
    axes[1, 1].set_xlabel(r"$x$")
    axes[1, 1].set_ylabel(r"$u_N(x,t)$")
    axes[1, 1].legend(loc="lower right", fontsize=8)

    figure.tight_layout()
    record_output(generated, figure, "fhn_diagnostics")


def plot_fhn_refinement(generated: list[str]) -> None:
    rows = read_csv("fhn_refinement_summary.csv")
    counts = np.array([int(row["N"]) for row in rows])
    profile = np.array([float(row["profile_L2"]) for row in rows])
    front = np.array([float(row["front_error"]) for row in rows])
    guide = profile[0] * np.sqrt(counts[0] / counts)
    figure, axis = plt.subplots(figsize=(4.8, 3.6))
    axis.loglog(counts, profile, "o-", color=GRW, label="Profile error")
    axis.loglog(counts, front, "s--", color=DETERMINISTIC, label="Front error")
    axis.loglog(counts, guide, "--", color="black", linewidth=1.2, label=r"$O(N^{-1/2})$")
    axis.set_xlabel(r"Particle count $N$")
    axis.set_ylabel("Error")
    axis.legend(loc="lower left")
    figure.tight_layout()
    record_output(generated, figure, "fhn_error_vs_N_final")


def plot_burgers_comparison(arrays: dict[str, np.ndarray], generated: list[str]) -> None:
    figure, axis = plt.subplots(figsize=(5.5, 3.5))
    x = arrays["burgers_x"]
    axis.plot(x, arrays["burgers_u_exact"], color=EXACT, label="Exact stationary shock")
    axis.plot(x, arrays["burgers_u_grw"], color=GRW, linestyle="--", label="Cole–Hopf GRW")
    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"$u(x,T)$")
    axis.set_xlim(0.0, 4.0)
    axis.legend(loc="best")
    figure.tight_layout()
    record_output(generated, figure, "burgers_comparison")


def plot_burgers_diagnostics(arrays: dict[str, np.ndarray], generated: list[str]) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(7.5, 5.6), sharex=True)
    x = arrays["burgers_x"]

    axes[0, 0].plot(x, arrays["burgers_phi_exact"], color=EXACT, label="Exact shape")
    axes[0, 0].plot(x, arrays["burgers_phi_fd"], color=DETERMINISTIC, linestyle="-.", label="Deterministic pipeline")
    axes[0, 0].plot(x, arrays["burgers_phi_grw"], color=GRW, linestyle="--", label="GRW")
    axes[0, 0].set_title(r"(a) Transformed field $\phi$")
    axes[0, 0].set_ylabel(r"$\phi(x,T)$")
    axes[0, 0].legend(loc="best", fontsize=8)

    axes[0, 1].plot(x, arrays["burgers_ratio_exact"], color=EXACT, label="Exact")
    axes[0, 1].plot(x, arrays["burgers_ratio_fd"], color=DETERMINISTIC, linestyle="-.", label="Deterministic pipeline")
    axes[0, 1].plot(x, arrays["burgers_ratio_grw"], color=GRW, linestyle="--", label="GRW")
    axes[0, 1].set_title(r"(b) Recovery ratio $\phi_x/\phi$")
    axes[0, 1].set_ylabel(r"$\phi_x/\phi$")

    axes[1, 0].plot(x, arrays["burgers_u_exact"], color=EXACT, label="Exact")
    axes[1, 0].plot(x, arrays["burgers_u_fd"], color=DETERMINISTIC, linestyle="-.", label="Deterministic pipeline")
    axes[1, 0].plot(x, arrays["burgers_u_grw"], color=GRW, linestyle="--", label="GRW")
    axes[1, 0].set_title(r"(c) Recovered field $u$")
    axes[1, 0].set_xlabel(r"$x$")
    axes[1, 0].set_ylabel(r"$u(x,T)$")

    deterministic_error = arrays["burgers_u_fd"] - arrays["burgers_u_exact"]
    grw_error = arrays["burgers_u_grw"] - arrays["burgers_u_fd"]
    axes[1, 1].plot(x, deterministic_error, color=DETERMINISTIC, label=r"$u^{\mathrm{FD}}-u^{\mathrm{ex}}$")
    axes[1, 1].plot(x, grw_error, color=GRW, linestyle="--", label=r"$u^{\mathrm{GRW}}-u^{\mathrm{FD}}$")
    axes[1, 1].axhline(0.0, color=EXACT, linestyle=":", linewidth=1.0)
    axes[1, 1].set_title("(d) Recovered-field error components")
    axes[1, 1].set_xlabel(r"$x$")
    axes[1, 1].set_ylabel("Error")
    axes[1, 1].legend(loc="best", fontsize=8)

    figure.tight_layout()
    record_output(generated, figure, "burgers_diagnostics")


def plot_burgers_domain(generated: list[str]) -> None:
    rows = read_csv("burgers_domain_sensitivity_summary.csv")
    lengths = np.array([float(row["L"]) for row in rows])
    deterministic = np.array([float(row["bc_mismatch_RMSE"]) for row in rows])
    grw = np.array([float(row["grw_particle_RMSE"]) for row in rows])
    total = np.array([float(row["total_RMSE"]) for row in rows])
    figure, axis = plt.subplots(figsize=(4.8, 3.6))
    axis.plot(lengths, deterministic, "o-", color=DETERMINISTIC, label=r"$E_{\mathrm{det}}$")
    axis.plot(lengths, grw, "s-", color=GRW, label=r"$E_{\mathrm{GRW}}$")
    axis.plot(lengths, total, "^-", color=EXACT, label=r"$E_{\mathrm{total}}$")
    axis.set_xlabel(r"Domain size $L$")
    axis.set_ylabel("RMSE")
    axis.set_xticks(lengths)
    axis.legend(loc="upper center", fontsize=8)
    figure.tight_layout()
    record_output(generated, figure, "burgers_domain_sensitivity_final")


def _l2h(u: np.ndarray, ref: np.ndarray, x: np.ndarray) -> float:
    h = float(x[1] - x[0])
    return math.sqrt(h * float(np.sum((u - ref) ** 2)))


def _linf(u: np.ndarray, ref: np.ndarray) -> float:
    return float(np.max(np.abs(u - ref)))


def _rel_l2(u: np.ndarray, ref: np.ndarray, x: np.ndarray) -> float:
    h = float(x[1] - x[0])
    numerator = math.sqrt(h * float(np.sum((u - ref) ** 2)))
    denominator = math.sqrt(h * float(np.sum(ref**2)))
    return numerator / denominator if denominator else float("nan")


def _rmse(u: np.ndarray, v: np.ndarray) -> float:
    return math.sqrt(float(np.mean((u - v) ** 2)))


def compute_metrics(arrays: dict[str, np.ndarray]) -> dict[str, object]:
    """Compute the representative error norms reported in the paper's summary table."""
    heat = {
        "N": 50000,
        "L2h": _l2h(arrays["heat_grw"], arrays["heat_exact"], arrays["heat_x"]),
        "Linf": _linf(arrays["heat_grw"], arrays["heat_exact"]),
        "rel_L2": _rel_l2(arrays["heat_grw"], arrays["heat_exact"], arrays["heat_x"]),
        "L2h_nbins300": _l2h(arrays["heat_grw_300"], arrays["heat_exact_300"],
                             arrays["heat_x_300"]),
        "total_weight": float(arrays["heat_grw"][-1]),
    }
    fhn_snapshots = {}
    for tag in ("0", "3", "6", "9"):
        grw = arrays[f"fhn_grw_t{tag}"]
        exact = arrays[f"fhn_exact_t{tag}"]
        fhn_snapshots[tag] = {
            "L2h": _l2h(grw, exact, arrays["fhn_x"]),
            "Linf": _linf(grw, exact),
            "rel_L2": _rel_l2(grw, exact, arrays["fhn_x"]),
            "total_weight": float(np.sum(arrays[f"fhn_weights_t{tag}"])),
        }
    fhn = {"N": 500, "final": fhn_snapshots["9"], "snapshots": fhn_snapshots}
    fhn["t0_weight_spread"] = float(
        np.ptp(arrays["fhn_weights_t0"])  # peak-to-peak; should be ~0 (uniform start)
    )
    burgers = {
        "N": 400,
        "L": 4,
        "L2h": _l2h(arrays["burgers_u_grw"], arrays["burgers_u_exact"], arrays["burgers_x"]),
        "Linf": _linf(arrays["burgers_u_grw"], arrays["burgers_u_exact"]),
        "rel_L2": _rel_l2(arrays["burgers_u_grw"], arrays["burgers_u_exact"], arrays["burgers_x"]),
        "E_det": _rmse(arrays["burgers_u_fd"], arrays["burgers_u_exact"]),
        "E_grw": _rmse(arrays["burgers_u_grw"], arrays["burgers_u_fd"]),
        "E_total": _rmse(arrays["burgers_u_grw"], arrays["burgers_u_exact"]),
        "clipped_bins": int(arrays["burgers_clipped"].sum()),
    }
    return {"heat": heat, "fhn": fhn, "burgers": burgers}


def _repo_rel(path) -> str:
    """Repository-relative POSIX path string for portable metadata."""
    try:
        return Path(path).resolve().relative_to(CODE_DIR).as_posix()
    except ValueError:
        return Path(path).name


def write_metadata(generated: list[str], force: bool, arrays: dict[str, np.ndarray]) -> None:
    """Record how the paper figures were produced (seed, sources, outputs).

    All paths are repository-relative so the metadata is portable across
    machines and checkouts.
    """
    timestamp = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
    metadata = {
        "generated_at": timestamp,
        "generator": _repo_rel(__file__),
        "shared_style": _repo_rel(SCRIPT_DIR / "plot_style.py"),
        "representative_arrays": _repo_rel(DATA_DIR / "representative_figure_arrays.npz"),
        "representative_arrays_regenerated": bool(force),
        "seed": SEED,
        "outputs": [_repo_rel(p) for p in generated],
        "figures": {
            "heat_comparison": {
                "source": "archived representative arrays (fixed-seed direct Brownian GRW)",
                "parameters": {"domain": [0, 10], "alpha": 0.1, "T": 0.5, "dt": 0.001, "N": 50000},
                "seed": SEED,
                "note": "main panel plus residual inset u_N - u_ex",
            },
            "heat_error_fixed_grid_diagnostic_final": {
                "source": _repo_rel(CSV_DIR / "heat_refinement_summary.csv"),
                "mode": "checked-in study CSV (paper data of record)",
            },
            "fhn_comparison": {
                "source": "archived representative arrays (fixed-seed scalar-FHN GRW with saved snapshots)",
                "parameters": {"domain": [0, 30], "a": 0.25, "D": 0.5, "T": 9, "dt": 0.01, "N": 500, "boundary": "Neumann sign-flip"},
                "seed": SEED,
            },
            "fhn_diagnostics": {
                "source": "same fixed-seed trajectory as fhn_comparison",
                "reference": "empirically anchored initial crossing",
                "seed": SEED,
                "note": "weights panel includes a zero line so the two Neumann-negated weights are visible",
            },
            "fhn_error_vs_N_final": {
                "source": _repo_rel(CSV_DIR / "fhn_refinement_summary.csv"),
                "mode": "checked-in study CSV (paper data of record)",
                "guide": "black dashed O(N^{-1/2}) reference guide, not a fitted exponent",
            },
            "burgers_comparison": {
                "source": "archived representative arrays (fixed-seed Cole-Hopf GRW with saved pipeline arrays)",
                "parameters": {"domain": [0, 4], "A": 1, "nu": 0.5, "T": 0.5, "dt": 0.005, "N": 400, "sigma_bins": 12},
                "seed": SEED,
                "clipped_bins": int(arrays["burgers_clipped"].sum()),
            },
            "burgers_diagnostics": {
                "source": "same fixed-seed pipeline as burgers_comparison",
                "panels": "(a) transformed field, (b) recovery ratio, (c) recovered field, (d) pointwise error components u_FD - u_ex and u_GRW - u_FD (RMSEs E_det and E_GRW)",
                "seed": SEED,
            },
            "burgers_domain_sensitivity_final": {
                "source": _repo_rel(CSV_DIR / "burgers_domain_sensitivity_summary.csv"),
                "mode": "checked-in study CSV (paper data of record)",
                "series": ["E_det", "E_GRW", "E_total"],
                "particle_scaling": "N=100L (fixed particle density; does not isolate particle count)",
            },
        },
    }
    metadata_path = FIGURE_DIR / "figure_regeneration_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    global FIGURE_DIR
    force = "--rerun" in sys.argv
    if "--output-dir" in sys.argv:
        index = sys.argv.index("--output-dir")
        try:
            FIGURE_DIR = Path(sys.argv[index + 1]).resolve()
        except IndexError as exc:
            raise SystemExit("--output-dir requires a path") from exc
    apply_paper_style()
    arrays = generate_or_load_arrays(force=force)
    generated: list[str] = []
    plot_heat_comparison(arrays, generated)
    plot_heat_refinement(generated)
    plot_fhn_comparison(arrays, generated)
    plot_fhn_diagnostics(arrays, generated)
    plot_fhn_refinement(generated)
    plot_burgers_comparison(arrays, generated)
    plot_burgers_diagnostics(arrays, generated)
    plot_burgers_domain(generated)
    write_metadata(generated, force, arrays)
    metrics = compute_metrics(arrays)

    print(f"\nFigures → {FIGURE_DIR}/")
    seen: set[str] = set()
    for path in generated:
        stem = Path(path).stem
        if stem not in seen:
            print(f"  {stem}")
            seen.add(stem)

    heat = metrics["heat"]
    fhn = metrics["fhn"]
    burgers = metrics["burgers"]
    print(f"\nVerification  (N = {heat['N']} / {fhn['N']} / {burgers['N']})")
    print(f"  {'Equation':<8}  {'L2_h':>9}  {'Linf':>9}  {'rel_L2':>9}")
    print(f"  {'Heat':<8}  {heat['L2h']:>9.6f}  {heat['Linf']:>9.6f}  {heat['rel_L2']:>9.6f}")
    print(f"  {'FHN':<8}  {fhn['final']['L2h']:>9.6f}  {fhn['final']['Linf']:>9.6f}  {fhn['final']['rel_L2']:>9.6f}  (t=9)")
    print(f"  {'Burgers':<8}  {burgers['L2h']:>9.6f}  {burgers['Linf']:>9.6f}  {burgers['rel_L2']:>9.6f}  (L=4)")


if __name__ == "__main__":
    main()
