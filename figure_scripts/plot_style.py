"""Shared matplotlib style for the GRW paper figures."""

from pathlib import Path

import matplotlib.pyplot as plt


EXACT = "#111111"
GRW = "#2166ac"
DETERMINISTIC = "#d6604d"
SECONDARY = "#762a83"
GUIDE = "#777777"


def apply_paper_style() -> None:
    """Apply one publication style without requiring a LaTeX installation."""
    plt.rcParams.update(
        {
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
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.bbox": "tight",
            "savefig.dpi": 300,
        }
    )


def save_pdf_png(figure, output_directory: Path, stem: str) -> list[str]:
    """Save matching vector and high-resolution raster copies."""
    output_directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix in ("pdf", "png"):
        path = output_directory / f"{stem}.{suffix}"
        figure.savefig(path, format=suffix, dpi=300)
        paths.append(str(path))
    plt.close(figure)
    return paths
