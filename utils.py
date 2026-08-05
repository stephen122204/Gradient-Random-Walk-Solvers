# matplotlib.use() must be called before pyplot is imported.
import sys
import os
import math
import datetime
from math import erf as _erf_scalar
import matplotlib

if sys.platform.startswith("darwin"):
    matplotlib.use("MacOSX")   # macOS native backend
else:
    matplotlib.use("TkAgg")    # cross-platform fallback

import matplotlib.pyplot as plt
import numpy as np

try:
    import config as config_module
except Exception:
    config_module = None
    print("[plot] Warning: could not import module-level config; will rely on cfg param.")

def _make_run_dir() -> str:
    stamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    path = os.path.join("outputs", stamp)
    os.makedirs(path, exist_ok=True)
    return path


# --------------------------------------------------------------------------- #
# Publication-quality plot style (matches the paper figures)
# --------------------------------------------------------------------------- #
_EXACT = "#111111"
_GRW = "#2166ac"
_SECONDARY = "#762a83"


def _apply_paper_style():
    """Apply publication rcParams to all subsequent matplotlib figures."""
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 13,
        "lines.linewidth": 1.9,
        "lines.markersize": 5.5,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.7,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "savefig.dpi": 300,
    })


# --------------------------------------------------------------------------- #
# Exact-solution helpers
# --------------------------------------------------------------------------- #

def _erf_vec(arr):
    """Element-wise erf with no scipy dependency."""
    a = np.asarray(arr, dtype=float)
    return np.array([_erf_scalar(float(v)) for v in a.ravel()]).reshape(a.shape)


def _exact_heat_step(x, t, alpha, x0, uL, uR):
    """Exact heat-equation solution for a Heaviside step IC on an infinite domain."""
    if t <= 0.0:
        return np.where(np.asarray(x, dtype=float) >= x0, float(uR), float(uL))
    scaled = (np.asarray(x, dtype=float) - x0) / (2.0 * math.sqrt(alpha * t))
    return float(uL) + float(uR - uL) * 0.5 * (1.0 + _erf_vec(scaled))


def _exact_fhn_wave(x, t, a_, x_center):
    """Exact FHN traveling-wave: u(x,t) = 1 / (1 + exp(-(x + theta*t - center)/2))."""
    theta = math.sqrt(2.0) * (0.5 - float(a_))
    xi = np.asarray(x, dtype=float) + theta * float(t) - float(x_center)
    return 1.0 / (1.0 + np.exp(-xi / 2.0))


def _exact_burgers_stationary_shock(x, nu, x_center, amplitude):
    """Exact stationary-shock solution: u = -A * tanh(A*(x - xc) / (2*nu))."""
    arg = float(amplitude) * (np.asarray(x, dtype=float) - float(x_center)) / (2.0 * float(nu))
    return -float(amplitude) * np.tanh(arg)


# --------------------------------------------------------------------------- #
# Config / BC helpers (unchanged)
# --------------------------------------------------------------------------- #

def _safe_str(x):
    try:
        return (x or "").strip().lower()
    except Exception:
        return ""


def _get_attr(obj, keys, default=None):
    """Return first existing attribute among possible names on obj."""
    if obj is None:
        return default
    for k in keys:
        if hasattr(obj, k):
            return getattr(obj, k)
    return default


def _bc_info(cfg=None):
    """
    Read BC info from cfg.boundary_conditions (preferred).
    Fall back to legacy attribute names if needed.
    """
    source = cfg if cfg is not None else config_module

    L = float(_get_attr(source, ["domain_size", "L", "length"], 1.0))

    bc = _get_attr(source, ["boundary_conditions"], None)
    if isinstance(bc, dict) and "LEFT" in bc and "RIGHT" in bc:
        ltype = _safe_str(bc["LEFT"].get("type", ""))
        rtype = _safe_str(bc["RIGHT"].get("type", ""))
        lval = float(bc["LEFT"].get("value", 0.0))
        rval = float(bc["RIGHT"].get("value", 0.0))
        print(f"[plot] BCs detected -> left: {ltype} ({lval}), right: {rtype} ({rval}), L={L}")
        return L, ltype, lval, rtype, rval

    # Legacy attribute-style configs (kept for compatibility)
    ltype = _safe_str(_get_attr(
        source,
        ["left_boundary_condition_type", "left_boundary_type", "left_bc_type", "left_bc"],
        ""
    ))
    rtype = _safe_str(_get_attr(
        source,
        ["right_boundary_condition_type", "right_boundary_type", "right_bc_type", "right_bc"],
        ""
    ))
    lval = float(_get_attr(source, ["left_boundary_value", "left_bc_value", "left_value"], 0.0))
    rval = float(_get_attr(source, ["right_boundary_value", "right_bc_value", "right_value"], 0.0))

    print(f"[plot] BCs detected -> left: {ltype} ({lval}), right: {rtype} ({rval}), L={L}")
    return L, ltype, lval, rtype, rval


def _both_dirichlet(cfg=None):
    _, ltype, _, rtype, _ = _bc_info(cfg)
    return ltype.startswith("dirichlet") and rtype.startswith("dirichlet")


def _both_neumann0(cfg=None):
    _, ltype, lval, rtype, rval = _bc_info(cfg)
    return (
        ltype.startswith("neumann")
        and rtype.startswith("neumann")
        and abs(float(lval)) == 0.0
        and abs(float(rval)) == 0.0
    )


def _extract_positions(results):
    n = len(results)
    print(f"[plot] number of globs: {n}")
    if n == 0:
        print("[plot] No results to plot.")
        return None

    try:
        positions = np.array([g["position"] for g in results], dtype=float)
    except Exception as e:
        print(f"[plot] Failed to read positions from results: {e}")
        return None

    if positions.size == 0 or not np.isfinite(positions).any():
        print("[plot] Positions are empty or non-finite.")
        return None

    return positions


# --------------------------------------------------------------------------- #
# Heat: gradient density (Neumann / mixed BCs)
# --------------------------------------------------------------------------- #
def _plot_heat_density(results, title_extra="", out_dir="outputs"):
    """
    Plot the weighted glob density, which approximates u_x(x, t).

    Each glob is weighted by its signed value when building the histogram, so the plot
    shows the empirical gradient density rather than a raw position count. This is the
    direct output of the GRW method before reconstruction.
    """
    positions = _extract_positions(results)
    if positions is None:
        return False

    try:
        weights = np.array([g["value"] for g in results], dtype=float)
    except Exception:
        weights = np.ones(positions.size, dtype=float)

    xmin, xmax = float(np.nanmin(positions)), float(np.nanmax(positions))
    if not np.isfinite(xmin) or not np.isfinite(xmax):
        print("[plot] x-range not finite; aborting plot.")
        return False
    if xmin == xmax:
        xmin -= 1e-6
        xmax += 1e-6

    n = positions.size
    nbins = int(max(25, min(200, n // 20)))
    counts, edges = np.histogram(positions, bins=nbins, range=(xmin, xmax), weights=weights)
    dx = edges[1] - edges[0]
    density = counts / (np.sum(np.abs(counts)) * dx) if np.sum(np.abs(counts)) > 0 else counts
    centers = 0.5 * (edges[:-1] + edges[1:])

    figure, axis = plt.subplots(figsize=(5.5, 3.5))
    axis.plot(centers, density, color=_GRW, linewidth=1.9)
    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"Glob density ($\approx u_x$)")
    title = "Heat GRW: gradient density"
    if title_extra:
        title += f"  {title_extra}"
    axis.set_title(title)

    ymax = (np.nanmax(np.abs(density)) if np.isfinite(density).any() else 1.0) * 1.1
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(-ymax, ymax)
    figure.tight_layout()

    out = os.path.join(out_dir, "heat_density.png")
    figure.savefig(out, dpi=300)
    print(f"[plot] Saved -> {out}")
    plt.show(block=True)
    return True


# --------------------------------------------------------------------------- #
# Heat: field reconstruction (Dirichlet BCs)
# --------------------------------------------------------------------------- #
def _plot_heat_dirichlet_as_field(results, cfg=None, out_dir="outputs"):
    """
    Reconstruct u(x, T) via cumulative sum of sorted glob values; overlay exact
    error-function solution when the IC is a step function.
    """
    positions = _extract_positions(results)
    if positions is None:
        return False

    try:
        values = np.array([g["value"] for g in results], dtype=float)
    except Exception as e:
        print(f"[plot] Failed to read glob values: {e}")
        return False

    L, _, uL, _, uR = _bc_info(cfg)
    L = float(L)
    uL = float(uL)
    uR = float(uR)

    n = positions.size
    nbins = int(max(100, min(400, n // 10)))
    edges = np.linspace(0.0, L, nbins + 1)
    bin_weights, _ = np.histogram(positions, bins=edges, weights=values)
    u_grw = uL + np.cumsum(bin_weights)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Overlay exact erf solution for step IC
    has_exact = (
        getattr(cfg, 'heat_ic_type', '') == 'step'
        and getattr(cfg, 'heat_jump_position', None) is not None
        and float(getattr(cfg, 'total_time', 0.0)) > 0.0
        and float(getattr(cfg, 'diff_constant', 0.0)) > 0.0
    )

    figure, axis = plt.subplots(figsize=(5.5, 3.5))

    if has_exact:
        x0 = float(cfg.heat_jump_position)
        alpha = float(cfg.diff_constant)
        T = float(cfg.total_time)
        u_exact = _exact_heat_step(centers, T, alpha, x0, uL, uR)
        axis.plot(centers, u_exact, color=_EXACT, label="Exact")
        axis.plot(centers, u_grw, color=_GRW, linestyle="--", label="GRW")
        axis.set_title("Heat GRW profile verification")
    else:
        axis.plot(centers, u_grw, color=_GRW, linewidth=1.9, label="GRW")
        axis.set_title(r"Heat GRW: reconstructed $u(x,T)$")

    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"$u(x,T)$")
    axis.set_xlim(0.0, L)
    u_lo, u_hi = min(uL, uR), max(uL, uR)
    axis.set_ylim(u_lo - 0.05, u_hi + 0.05)
    axis.legend(loc="upper left")
    figure.tight_layout()

    out = os.path.join(out_dir, "heat_field_dirichlet.png")
    figure.savefig(out, dpi=300)
    print(f"[plot] Saved -> {out}")
    plt.show(block=True)
    return True


# --------------------------------------------------------------------------- #
# Burgers: u(x,T)
# --------------------------------------------------------------------------- #
def _plot_burgers_field(results, cfg=None, out_dir="outputs"):
    """Plot u(x,T); overlay exact solution for stationary_shock and traveling_wave ICs."""
    positions = _extract_positions(results)
    if positions is None:
        return False

    try:
        u = np.array([g["value"][0] for g in results], dtype=float)
    except Exception as e:
        print(f"[plot] Failed to read Burgers u from results: {e}")
        return False

    order = np.argsort(positions)
    xs = positions[order]
    us = u[order]

    ic_type = (getattr(cfg, 'burgers_ic_type', '') or '').lower()
    nu = float(getattr(cfg, 'diff_constant', 0.5))
    L = float(getattr(cfg, 'domain_size', float(xs[-1]) if len(xs) > 0 else 4.0))

    figure, axis = plt.subplots(figsize=(5.5, 3.5))

    if ic_type == 'stationary_shock':
        amplitude = float(getattr(cfg, 'burgers_ic_amplitude', 1.0) or 1.0)
        x_center = float(getattr(cfg, 'burgers_ic_x_center', None) or L / 2.0)
        u_exact = _exact_burgers_stationary_shock(xs, nu, x_center, amplitude)
        axis.plot(xs, u_exact, color=_EXACT, label="Exact stationary shock")
        axis.plot(xs, us, color=_GRW, linestyle="--", label="Cole–Hopf GRW")
        axis.set_title("Cole–Hopf Burgers shock reconstruction")

    elif ic_type == 'traveling_wave':
        T = float(getattr(cfg, 'total_time', 0.5))
        x_center = float(getattr(cfg, 'burgers_ic_x_center', None) or L / 2.0)
        u_exact = 1.0 - 2.0 * np.sqrt(nu) * np.tanh((xs - x_center - T) / np.sqrt(nu))
        axis.plot(xs, u_exact, color=_EXACT, label="Exact traveling wave")
        axis.plot(xs, us, color=_GRW, linestyle="--", label="Cole–Hopf GRW")
        axis.set_title("Cole–Hopf Burgers traveling wave")

    else:
        axis.plot(xs, us, color=_GRW, linewidth=1.9, label="GRW")
        axis.set_title(r"Burgers: $u(x,T)$")

    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"$u(x,T)$")
    axis.set_xlim(0.0, L)
    axis.legend(loc="best")
    figure.tight_layout()

    out = os.path.join(out_dir, "burgers_u.png")
    figure.savefig(out, dpi=300)
    print(f"[plot] Saved -> {out}")
    plt.show(block=True)
    return True


# --------------------------------------------------------------------------- #
# FHN helpers: reconstruction and multi-/single-panel figures
# --------------------------------------------------------------------------- #
def _reconstruct_u_grid(xs_sorted, ws_sorted, x_grid):
    """Interpolate the sorted cumulative-weight sum onto a uniform output grid."""
    u_cum = np.cumsum(ws_sorted)
    return np.interp(x_grid, xs_sorted, u_cum,
                     left=0.0, right=float(u_cum[-1]) if len(u_cum) > 0 else 0.0)


def _plot_fhn_multi_panel(snapshots, a_, x_center, L, out_dir="outputs"):
    """2x2 panel: GRW vs exact traveling wave at 4 time snapshots."""
    snap_times = sorted(snapshots.keys())[:4]
    x_grid = np.linspace(0.0, L, 500)

    figure, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True, sharey=True)
    axes_list = list(axes.flat)

    for i, t in enumerate(snap_times):
        ax = axes_list[i]
        xs_t, ws_t = snapshots[t]
        u_grw_grid = _reconstruct_u_grid(xs_t, ws_t, x_grid)
        u_exact_grid = _exact_fhn_wave(x_grid, t, a_, x_center)
        ax.plot(x_grid, u_exact_grid, color=_EXACT, label="Exact")
        ax.plot(x_grid, u_grw_grid, color=_GRW, linestyle="--", label="GRW")
        ax.set_title(rf"$t={t:g}$")
        ax.set_xlim(0.0, L)
        ax.set_ylim(-0.05, 1.05)

    for i in range(len(snap_times), len(axes_list)):
        axes_list[i].set_visible(False)

    axes[0, 0].legend(loc="upper left")
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$x$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$u(x,t)$")

    figure.suptitle("Scalar FHN traveling-front verification")
    figure.tight_layout(rect=(0, 0, 1, 0.96))

    out = os.path.join(out_dir, "fhn_uv.png")
    figure.savefig(out, dpi=300)
    print(f"[plot] Saved -> {out}")
    plt.show(block=True)
    return True


def _plot_fhn_single(xs_sorted, ws_sorted, a_, x_center, L, cfg=None, out_dir="outputs"):
    """Single-panel final state: GRW reconstruction with exact overlay."""
    T = float(getattr(cfg, 'total_time', 0.0))
    x_grid = np.linspace(0.0, L, 500)
    u_grw_grid = _reconstruct_u_grid(xs_sorted, ws_sorted, x_grid)

    figure, axis = plt.subplots(figsize=(5.5, 3.5))
    if T > 0.0:
        u_exact_grid = _exact_fhn_wave(x_grid, T, a_, x_center)
        axis.plot(x_grid, u_exact_grid, color=_EXACT, label="Exact")
        axis.plot(x_grid, u_grw_grid, color=_GRW, linestyle="--", label="GRW")
        axis.set_title(rf"FHN GRW profile ($t={T:g}$)")
    else:
        axis.plot(x_grid, u_grw_grid, color=_GRW, linewidth=1.9, label="GRW")
        axis.set_title("FitzHugh-Nagumo: GRW reconstruction")

    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"$u(x,t)$")
    axis.set_xlim(0.0, L)
    axis.set_ylim(-0.05, 1.05)
    axis.legend(loc="upper left")
    figure.tight_layout()

    out = os.path.join(out_dir, "fhn_uv.png")
    figure.savefig(out, dpi=300)
    print(f"[plot] Saved -> {out}")
    plt.show(block=True)
    return True


# --------------------------------------------------------------------------- #
# FHN: main dispatcher
# --------------------------------------------------------------------------- #
def _plot_fhn_fields(results, cfg=None, out_dir="outputs"):
    """
    Plot the FHN solution u(x,t).

    Scalar GRW: uses multi-time 2x2 panel when snapshot data is available
    (set by simulate_fitzhugh_nagumo_grw in simulation.py), otherwise shows
    the final state. Overlays the exact traveling-wave solution in both cases.

    Legacy two-component: plots u and v on a single panel.
    """
    positions = _extract_positions(results)
    if positions is None:
        return False

    order = np.argsort(positions)
    xs = positions[order]
    raw_values = [g["value"] for g in results]
    is_scalar = isinstance(raw_values[0], (float, int, np.floating))

    if not is_scalar:
        # Legacy two-component path
        try:
            uv = np.array(raw_values, dtype=float)[order]
            u = uv[:, 0]
            v = uv[:, 1]
        except Exception as e:
            print(f"[plot] Failed to read FHN (u,v) from results: {e}")
            return False
        figure, axis = plt.subplots(figsize=(6.0, 4.0))
        axis.plot(xs, u, color=_GRW, linewidth=1.9, label="u")
        axis.plot(xs, v, color=_SECONDARY, linewidth=1.9, label="v")
        axis.set_xlabel(r"$x$")
        axis.set_ylabel("state")
        axis.set_title(r"FitzHugh-Nagumo: $u(x,t)$, $v(x,t)$")
        axis.legend()
        figure.tight_layout()
        out = os.path.join(out_dir, "fhn_uv.png")
        figure.savefig(out, dpi=300)
        print(f"[plot] Saved -> {out}")
        plt.show(block=True)
        return True

    # Scalar GRW path
    w = np.array(raw_values, dtype=float)[order]
    a_ = float(getattr(cfg, '_fhn_a', None) or getattr(cfg, 'a', 0.25) or 0.25)
    L = float(getattr(cfg, 'domain_size', float(xs[-1]) if len(xs) > 0 else 30.0))

    x_center = getattr(cfg, '_fhn_x_center', None)
    if x_center is None:
        u_cum = np.cumsum(w)
        idx05 = int(np.clip(np.searchsorted(u_cum, 0.5), 0, len(u_cum) - 1))
        x_center = float(xs[idx05]) if len(xs) > 0 else L / 2.0

    snapshots = getattr(cfg, '_fhn_snapshots', None)
    if snapshots is not None and len(snapshots) >= 2:
        print("[plot] Plotting FHN multi-time verification (2×2 panel)...")
        return _plot_fhn_multi_panel(snapshots, a_, x_center, L, out_dir=out_dir)
    else:
        print("[plot] Plotting FHN single-panel final state...")
        return _plot_fhn_single(xs, w, a_, x_center, L, cfg, out_dir=out_dir)


# --------------------------------------------------------------------------- #
# Main entry point
# --------------------------------------------------------------------------- #
def plot_results(results, equation_type: str, cfg=None):
    """
    Plotting entry point used by main.py.

    Heat routing:
    - Dirichlet–Dirichlet BCs -> reconstructed u(x,T); overlay exact erf for step IC.
    - Neumann or mixed BCs    -> gradient density plot (approx. u_x).

    Burgers:
    - u(x,T) with exact overlay for stationary_shock and traveling_wave ICs.

    FitzHugh–Nagumo (scalar GRW):
    - 2x2 panel at 4 time snapshots when snapshot data is available.
    - Single-panel final state with exact overlay as fallback.
    """
    _apply_paper_style()
    out_dir = _make_run_dir()
    eq = (equation_type or "").strip().lower()

    if eq == "heat":
        if _both_dirichlet(cfg):
            print("[plot] Dirichlet–Dirichlet BCs → reconstructing field u(x,T).")
            ok = _plot_heat_dirichlet_as_field(results, cfg, out_dir=out_dir)
        else:
            extra = "(Neumann-0 reflecting)" if _both_neumann0(cfg) else ""
            print("[plot] Non-Dirichlet BCs → plotting gradient density (≈ u_x).")
            ok = _plot_heat_density(results, title_extra=extra, out_dir=out_dir)
        if not ok:
            print("[plot] Heat plot failed or had no data.")
        else:
            print(f"[plot] Output → {out_dir}/")
        return

    if eq in {"fitzhugh-nagumo", "fitzhugh–nagumo", "fitzhugh", "nagumo"}:
        print("[plot] Plotting FitzHugh–Nagumo fields...")
        ok = _plot_fhn_fields(results, cfg, out_dir=out_dir)
        if not ok:
            print("[plot] FHN plot failed or had no data.")
        else:
            print(f"[plot] Output → {out_dir}/")
        return

    if eq in {"burgers", "burger", "burgers'"}:
        print("[plot] Plotting Burgers field...")
        ok = _plot_burgers_field(results, cfg, out_dir=out_dir)
        if not ok:
            print("[plot] Burgers plot failed or had no data.")
        else:
            print(f"[plot] Output → {out_dir}/")
        return

    print(f"[plot] Unknown equation type for plotting: {equation_type!r}")
