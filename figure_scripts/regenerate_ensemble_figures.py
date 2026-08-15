#!/usr/bin/env python3
"""Regenerate all six ensemble/control paper figures from pinned data.

  heat_bias_spread_total_vs_N   from pinned_ensembles/heat_extended/summary_by_N.csv
  fhn_convergence               from pinned_ensembles/fhn_extended/summary_by_N.csv
  heat_grid_paired              from pinned_ensembles/heat_grid_paired/summary.json
  burgers_decoupled             from pinned_ensembles/burgers_controls/summary.json
  burgers_boundary_domain       from pinned_ensembles/burgers_controls/summary.json
  burgers_perturbation_response from pinned_ensembles/burgers_controls/summary.json

Titleless, legends only, publication style matching regenerate_paper_figures.
The study scripts also write diagnostic copies, but this script is the
canonical paper-figure path. Output:
output/final_prepublication_tests/paper_figures/
"""
import csv
import json
import os

os.environ.setdefault('SOURCE_DATE_EPOCH', '1704067200')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIN = os.path.join(ROOT, 'pinned_ensembles')
OUT = os.path.join(ROOT, 'output', 'final_prepublication_tests', 'paper_figures')

EXACT = '#111111'
GRW = '#2166ac'
DET = '#d6604d'
SECONDARY = '#762a83'
GUIDE = '#777777'


def _apply_style():
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Computer Modern Roman', 'Times New Roman', 'DejaVu Serif'],
        'mathtext.fontset': 'cm',
        'axes.labelsize': 11, 'axes.titlesize': 12,
        'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
        'lines.linewidth': 1.9, 'lines.markersize': 5.5,
        'axes.grid': True, 'grid.alpha': 0.25, 'grid.linewidth': 0.7,
        'axes.spines.top': True, 'axes.spines.right': True,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.bbox': 'tight',
    })


def _save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + '.pdf'))
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=150)
    plt.close(fig)


def _csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def _json(path):
    with open(path) as f:
        return json.load(f)


def heat_bias_spread_total_vs_N():
    rows = _csv(os.path.join(PIN, 'heat_extended', 'summary_by_N.csv'))
    N = np.array([float(r['N']) for r in rows])
    bias = np.array([float(r['E_bias']) for r in rows])
    spr = np.array([float(r['E_spread']) for r in rows])
    tot = np.array([float(r['E_total']) for r in rows])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(N, tot, 's-', color=EXACT, label=r'$E_\mathrm{total}$')
    ax.loglog(N, spr, 'o-', color=GRW, label=r'$E_\mathrm{spread}$')
    ax.loglog(N, bias, '^-', color=DET, label=r'$E_\mathrm{bias}$')
    xg = np.array([N.min(), N.max()])
    i_mid = len(N) // 2
    ax.loglog(xg, spr[i_mid] * (xg / N[i_mid]) ** (-0.5), ':', color=GUIDE,
              label=r'$N^{-1/2}$')
    ax.set_xlabel('N')
    ax.set_ylabel(r'$L^2$ error')
    ax.legend()
    _save(fig, 'heat_bias_spread_total_vs_N')


def fhn_convergence():
    rows = _csv(os.path.join(PIN, 'fhn_extended', 'summary_by_N.csv'))
    N = np.array([float(r['N']) for r in rows])
    panels = [
        (np.array([float(r['l2_mean']) for r in rows]),
         r'profile $L^2$ error', GRW, 'o'),
        (np.array([float(r['ce_mean']) for r in rows]),
         r'front-center error', DET, 's'),
        (np.array([float(r['ap_mean']) for r in rows]),
         r'aligned-profile error', SECONDARY, 'd'),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    for i, (y, lab, c, m) in enumerate(panels):
        ax = axes[i]
        ax.loglog(N, y, marker=m, color=c, ls='-', label=lab)
        ax.loglog(N, y[0] * (N / N[0]) ** -0.5, '--', color=EXACT, lw=1.2,
                  label=r'$N^{-1/2}$ guide')
        ax.set_xlabel(r'$N$')
        if i == 0:
            ax.set_ylabel('error')
        ax.legend(fontsize=8)
        ax.text(0.03, 0.03, f'({chr(97 + i)})', transform=ax.transAxes,
                fontsize=11)
    fig.tight_layout()
    _save(fig, 'fhn_convergence')


def heat_grid_paired():
    data = _json(os.path.join(PIN, 'heat_grid_paired', 'summary.json'))
    results = data['per_N']
    control = data['operator_control']
    N = np.array([float(r['N']) for r in results['coupled']])
    styles = {
        'coupled':   dict(color='#1b7837', marker='^', ls='-'),
        'fixed300':  dict(color=GRW, marker='s', ls='-'),
        'fixed300e': dict(color='#4393c3', marker='o', ls='--'),
        'fixed400':  dict(color=SECONDARY, marker='d', ls='-'),
        'fixed400e': dict(color='#c51b7d', marker='v', ls='--'),
    }
    labels = {
        'coupled': 'coupled grid ($M=N$)',
        'fixed300': 'fixed $M=300$, center-compare',
        'fixed300e': 'fixed $M=300$, aligned',
        'fixed400': 'fixed $M=400$, center-compare',
        'fixed400e': 'fixed $M=400$, aligned',
    }
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2))
    for arm in ('coupled', 'fixed300', 'fixed300e'):
        y = np.array([float(r['E_total']) for r in results[arm]])
        axL.loglog(N, y, lw=1.6, ms=5, label=labels[arm], **styles[arm])
    guide = np.array([N.min(), N.max()])
    y0 = float(results['coupled'][0]['E_total'])
    axL.loglog(guide, y0 * (guide / N[0]) ** -0.5, ':', color=GUIDE,
               lw=1.2, label=r'$N^{-1/2}$ guide')
    axL.set_xlabel(r'$N$')
    axL.set_ylabel(r'$E_{\mathrm{total}}$')
    axL.legend(fontsize=8)
    axL.text(0.02, 0.02, '(a)', transform=axL.transAxes, fontsize=11)

    for arm in ('fixed300', 'fixed400', 'fixed300e', 'fixed400e', 'coupled'):
        y = np.array([float(r['E_bias']) for r in results[arm]])
        axR.loglog(N, y, lw=1.4, ms=4, label=labels[arm], **styles[arm])
    axR.axhline(float(control['M300']['center_compare']), color=GRW,
                ls=':', lw=1.1, label='operator control, center')
    axR.axhline(float(control['M400']['center_compare']), color=SECONDARY,
                ls=':', lw=1.1)
    axR.axhline(float(control['M300']['edge_compare']), color=EXACT,
                ls=':', lw=1.1, label='operator control, aligned')
    axR.set_xlabel(r'$N$')
    axR.set_ylabel(r'$E_{\mathrm{bias}}$')
    axR.legend(fontsize=7)
    axR.text(0.02, 0.02, '(b)', transform=axR.transAxes, fontsize=11)
    fig.tight_layout()
    _save(fig, 'heat_grid_paired')


def burgers_decoupled():
    data = _json(os.path.join(PIN, 'burgers_controls', 'summary.json'))
    rows = data['decoupled']
    groups = [[r for r in rows if r['part'] == key] for key in ('A', 'B', 'C')]
    x_keys = ('P', 'M', 'sigma_x')
    x_labels = (r'initialization-grid points $P$',
                r'output-grid points $M$',
                r'physical bandwidth $\sigma_x$')
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    for i, (ax, group, x_key, x_label) in enumerate(
            zip(axes, groups, x_keys, x_labels)):
        x = np.array([float(r[x_key]) for r in group])
        ax.errorbar(x, [float(r['rmse_total_mean']) for r in group],
                    yerr=[float(r['rmse_total_std']) for r in group],
                    fmt='o-', color=EXACT, ms=4, lw=1.5, capsize=2,
                    label='vs exact shock')
        ax.errorbar(x, [float(r['rmse_particle_mean']) for r in group],
                    yerr=[float(r['rmse_particle_std']) for r in group],
                    fmt='s--', color=GRW, ms=4, lw=1.5, capsize=2,
                    label='vs deterministic reference')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(x_label)
        if i == 0:
            ax.set_ylabel('RMSE')
        ax.legend(fontsize=8)
        ax.text(0.03, 0.03, f'({chr(97 + i)})', transform=ax.transAxes,
                fontsize=11)
    fig.tight_layout()
    _save(fig, 'burgers_decoupled')


def burgers_boundary_domain():
    data = _json(os.path.join(PIN, 'burgers_controls', 'summary.json'))
    boundary = data['boundary_controls']
    domain = data['domain_multiseed']
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.0))
    L = [float(r['L']) for r in boundary]
    axL.semilogy(L, [float(r['E_det_pinned_rmse']) for r in boundary],
                 's-', color=DET, ms=5, lw=1.6, label='pinned endpoints')
    axL.semilogy(L, [float(r['E_det_exactBC_rmse']) for r in boundary],
                 'o-', color='#1b7837', ms=5, lw=1.6,
                 label='exact transformed boundary data')
    axL.semilogy(L, [float(r['exact_phi_recovery_rmse']) for r in boundary],
                 'd--', color=GUIDE, ms=5, lw=1.2,
                 label='exact transformed field')
    axL.set_xlabel(r'domain size $L$')
    axL.set_ylabel('deterministic RMSE')
    axL.legend(fontsize=8)
    axL.text(0.03, 0.03, '(a)', transform=axL.transAxes, fontsize=11)

    L = [float(r['L']) for r in domain]
    axR.errorbar(L, [float(r['E_total_mean']) for r in domain],
                 yerr=[float(r['E_total_std']) for r in domain],
                 fmt='^-', color=EXACT, ms=5, lw=1.6, capsize=3,
                 label=r'$E_{\mathrm{total}}$ (mean $\pm$ std)')
    axR.errorbar(L, [float(r['E_grw_mean']) for r in domain],
                 yerr=[float(r['E_grw_std']) for r in domain],
                 fmt='o-', color=GRW, ms=5, lw=1.6, capsize=3,
                 label=r'$E_{\mathrm{GRW}}$ (mean $\pm$ std)')
    axR.plot(L, [float(r['E_det_rmse']) for r in domain],
             's--', color=DET, ms=5, lw=1.4,
             label=r'$E_{\mathrm{det}}$ (deterministic)')
    axR.set_xlabel(r'domain size $L$')
    axR.set_ylabel('RMSE')
    axR.legend(fontsize=8)
    axR.text(0.03, 0.03, '(b)', transform=axR.transAxes, fontsize=11)
    fig.tight_layout()
    _save(fig, 'burgers_boundary_domain')


def burgers_perturbation_response():
    data = _json(os.path.join(PIN, 'burgers_controls', 'summary.json'))
    response = data['perturbation_response']
    rows = response['rows']
    amp = [float(r['amplitude']) for r in rows]
    fig, ax = plt.subplots(figsize=(6.6, 4.5))
    ax.errorbar(amp, [float(r['white_l2_mean']) for r in rows],
                yerr=[float(r['white_l2_std']) for r in rows],
                fmt='o-', color=GRW, ms=5, lw=1.6, capsize=3,
                label='white perturbation')
    ax.errorbar(amp, [float(r['smooth_l2_mean']) for r in rows],
                yerr=[float(r['smooth_l2_std']) for r in rows],
                fmt='s-', color='#1b7837', ms=5, lw=1.6, capsize=3,
                label='kernel-smoothed perturbation')
    ax.plot([float(response['measured_phi_rmse'])],
            [float(response['measured_u_l2'])], marker='*', color=DET,
            ms=14, ls='none', label='measured particle pipeline')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'perturbation amplitude on $\phi$ (RMS)')
    ax.set_ylabel(r'$L^2$ error of recovered $u$')
    ax.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, 'burgers_perturbation_response')


def main():
    _apply_style()
    heat_bias_spread_total_vs_N()
    fhn_convergence()
    heat_grid_paired()
    burgers_decoupled()
    burgers_boundary_domain()
    burgers_perturbation_response()
    print(f'ensemble-derived figures written to {OUT}')


if __name__ == '__main__':
    main()
