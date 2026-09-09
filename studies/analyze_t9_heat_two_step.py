"""Prediction-versus-observation analysis for the two-step heat study.

Reads output/heat_two_step/{predictions,ensembles,validation}.json, compares the
measured ensemble statistics with the predictions of the weighted identity, and
writes Figure_11.pdf. Each prediction carries a finite-ensemble fluctuation scale
computed from the exact covariance of the reconstruction, so that agreement is
assessed against the scatter expected at S realizations rather than against zero.

Run from the repository root after study_t9_heat_two_step.py run (and validate).
"""
import json, os, sys
from math import sqrt
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'studies'))
sys.path.insert(0, os.path.join(ROOT, 'figure_scripts'))
from plot_style import apply_paper_style           # noqa: E402
from study_t9_heat_two_step import (grid, covariance, m_finite, u_inf,  # noqa: E402
                                    M_LIST, N_LIST)
apply_paper_style()

OUT = os.path.join(ROOT, 'output', 'heat_two_step')
pred = json.load(open(f'{OUT}/predictions.json')); ens = json.load(open(f'{OUT}/ensembles.json'))
pM = {p['M']: p for p in pred['per_M']}
rows = ens['rows']; S = rows[0]['S']
REF = {'finite': m_finite, 'infinite': u_inf}

def predicted(M, ref, conv, N, stat):
    p = pM[M]; B2 = p[f'{ref}_{conv}']['B2']; V = p['V_h']
    return {'E_total': sqrt(B2 + V/N), 'E_bias': sqrt(B2 + V/(N*S)),
            'E_spread': sqrt((S-1)/S * V/N)}[stat]

def sd_total(M, ref, conv, N):
    """Standard deviation of the measured E_total from the exact covariance.

    Under a Gaussian approximation for the reconstruction errors, the squared
    statistic averaged over S realizations has variance (4 q'Sig q + 2 tr Sig^2)/S,
    and the delta method converts that to a scale for the norm itself.
    """
    edges, be, ce, h = grid(M)
    Sig = covariance(M, N)
    q = sqrt(h) * (m_finite(be) - REF[ref](be if conv == 'edge' else ce))
    var_sq = (4.0 * float(q @ Sig @ q) + 2.0 * float(np.sum(Sig * Sig))) / S
    Et = predicted(M, ref, conv, N, 'E_total')
    return sqrt(var_sq) / (2.0 * Et)

print("=== measured E_total against prediction, in units of the predicted fluctuation sd")
worst = 0.0; devs = []
for M in M_LIST:
    for ref in ('finite', 'infinite'):
        for conv in ('edge', 'center'):
            sel = sorted([r for r in rows if r['M'] == M and r['reference'] == ref and r['convention'] == conv],
                         key=lambda r: r['N'])
            z = [(r['E_total'] - predicted(M, ref, conv, r['N'], 'E_total')) / sd_total(M, ref, conv, r['N']) for r in sel]
            devs += z; worst = max(worst, max(abs(np.array(z))))
            print(f"M={M:3d} {ref:8s} {conv:6s}: " + " ".join(f"{v:+.2f}" for v in z))
print(f"\nall {len(devs)} comparisons: max |deviation| {worst:.2f} sd, rms {np.sqrt(np.mean(np.square(devs))):.2f} sd")

print("\n=== predicted floor B, transition count N*, and measured E_total (finite reference)")
for M in M_LIST:
    p = pM[M]
    for conv in ('center', 'edge'):
        B = p[f'finite_{conv}']['B']; Ns = p[f'finite_{conv}']['N_star']
        sel = sorted([r for r in rows if r['M'] == M and r['reference'] == 'finite' and r['convention'] == conv],
                     key=lambda r: r['N'])
        near = min(sel, key=lambda r: abs(np.log(r['N']/Ns))) if Ns else None
        s = f"M={M:3d} {conv:6s} B={B:.3e}"
        s += f" N*={Ns:8.0f} nearest N={near['N']:6d}: measured {near['E_total']:.5f} predicted {predicted(M,'finite',conv,near['N'],'E_total'):.5f}" if Ns else "  no floor (B=0), sampling only"
        print(s)

print("\n=== ratio of measured to predicted, aggregated")
for stat in ('E_total', 'E_spread', 'E_bias'):
    r_ = [r[stat]/predicted(r['M'], r['reference'], r['convention'], r['N'], stat) for r in rows]
    print(f"{stat:9s} min {min(r_):.3f}  median {np.median(r_):.3f}  max {max(r_):.3f}")

# ---- figure ----
fig, axes = plt.subplots(1, 2, figsize=(10, 4.0))
colors = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd']
Nline = np.logspace(np.log10(min(N_LIST)), np.log10(max(N_LIST)), 200)
for ax, ref, title in ((axes[0], 'finite', '(a) finite-interval reference'),
                       (axes[1], 'infinite', '(b) infinite-line reference')):
    for ci, M in enumerate(M_LIST):
        for conv, mk, ls in (('center', 'o', '-'), ('edge', 's', '--')):
            sel = sorted([r for r in rows if r['M'] == M and r['reference'] == ref and r['convention'] == conv],
                         key=lambda r: r['N'])
            ax.plot([r['N'] for r in sel], [r['E_total'] for r in sel], mk, color=colors[ci], ms=4, ls='none')
            ax.plot(Nline, [predicted(M, ref, conv, n, 'E_total') for n in Nline], ls, color=colors[ci], lw=1.0)
        if ref == 'finite':
            ax.plot([], [], 'o-', color=colors[ci], ms=4, lw=1.0, label=f'$M={M}$')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('$N$')
    ax.set_title(title, fontsize=10); ax.grid(True, which='both', alpha=0.3)
axes[0].set_ylabel(r'$E_{\mathrm{total}}$')
h1, l1 = axes[0].get_legend_handles_labels()
h1 += [plt.Line2D([], [], color='0.35', ls='-', marker='o', ms=4, lw=1.0),
       plt.Line2D([], [], color='0.35', ls='--', marker='s', ms=4, lw=1.0)]
l1 += ['bin centers', 'bin edges']
axes[0].legend(h1, l1, fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig(f'{OUT}/Figure_11.pdf'); print("\nFigure_11.pdf written")

if os.path.exists(f'{OUT}/validation.json'):
    val = json.load(open(f'{OUT}/validation.json')); spec = val['spec']; M = spec['M']; N = spec['N']
    print(f"\n=== held-out accuracy target: RMS error {spec['target_rms_error']} against the infinite-line reference")
    print(f"model rules out the bin-center convention before running: floor {spec['center_convention']['B']:.5f} exceeds the target")
    for r in val['rows']:
        if r['M'] == M and r['reference'] == 'infinite':
            pr = predicted(M, 'infinite', r['convention'], N, 'E_total')
            print(f"  {r['convention']:6s} N={N}: measured {r['E_total']:.5f}  predicted {pr:.5f}  "
                  f"ratio {r['E_total']/pr:.3f}  ({(r['E_total']-pr)/sd_total(M,'infinite',r['convention'],N):+.2f} sd)")
