#!/usr/bin/env python3
"""
One-at-a-time sweeps: hold oracle DE optimum fixed, vary each BOUNDS variable on a uniform grid,
compare global quadratic RSM (lstsq on BBD) vs low_fidelity ODE surrogate; boxplot observed Po.
Objective and symbols match VARIABLE_SELECTION.md / wetlab_sim.BOUNDS order:
  (T_prod, IPTG, X_ind, F_s, t_f).
"""
import argparse
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from wetlab_sim import BOUNDS, run_batch
from optimization import _bb_Z01, _normalize, _quad_features, low_fidelity, method_de_sim
warnings.filterwarnings('ignore')
os.makedirs('plots', exist_ok=True)
os.makedirs('runs', exist_ok=True)
VAR_LABELS = ['T_prod', 'IPTG', 'X_ind', 'F_s', 't_f']  # 1:1 with MODEL.md / optimization decision vector


def fit_global_quadratic_bbd(y_train_kind, seed_po):
    """
    Full-factor quadratic in coded [0,1]^5 on global BOUNDS; training Y from LF (default) or wet-lab Po
    mirrors method_doe_rsm's design row order but avoids the Stage-2 shrink box (study is global sensitivity).
    """
    z01 = _bb_Z01()  # canonical 46-run Box-Behnken on [0,1]^5 — same lattice as OPTIMIZATION_PROTOCOL §5.2 Stage 1
    lo, hi = BOUNDS[:, 0], BOUNDS[:, 1]
    xbbd = lo + z01 * (hi - lo)  # physical corners/face-centres of global box — training sites for lstsq surface
    if y_train_kind == 'lf':
        y = np.array([low_fidelity(xbbd[i]) for i in range(len(xbbd))], dtype=float)  # aligns RSM vs LF comparison on slices
    else:
        y = np.array([run_batch(xbbd[i], seed=int(seed_po) + i, scatter=False, noisy=True)[1] for i in range(len(xbbd))], dtype=float)  # Method-2 noisy target — stochastic training
    feats = _quad_features(z01)  # same basis as optimization.method_doe_rsm one_stage (coded coords, not _normalize rows)
    coefs, *_ = np.linalg.lstsq(feats, y, rcond=None)
    return coefs


def rsm_predict(coefs, x):
    """Evaluate quadratic at physical x ∈ R^5 (global BOUNDS encoding)."""
    z = _normalize(np.asarray(x, dtype=float).reshape(1, -1))  # coded unit box — predictions match coefs fitted on z01 lattice
    return float(max(0.0, _quad_features(z) @ coefs))  # titre cannot go negative — clip for plotting


def sweep_one_at_a_time(x_opt, coefs, n_grid, seed_base):
    """
    Equal spacing in each BOUNDS dimension; freeze other coords at oracle optimum.
    Seeds ladder keeps measure() noise repeatable across runs (wet-lab stand-in).
    """
    po_slices, lf_slices, rsm_slices = [], [], []
    k = 0
    for j in range(5):
        vals = np.linspace(BOUNDS[j, 0], BOUNDS[j, 1], int(n_grid))  # inclusive endpoints — scans full protocol range per variable
        po_j, lf_j, rsm_j = [], [], []
        for v in vals:
            x = np.asarray(x_opt, dtype=float).copy()
            x[j] = float(v)
            lf_j.append(low_fidelity(x))  # biased ODE surrogate — Kennedy–O'Hagan LF track in MFBO (OPTIMIZATION.md)
            rsm_j.append(rsm_predict(coefs, x))  # least-square quadratic — DoE+RSM analogue on global bounds
            po_j.append(run_batch(x, seed=int(seed_base) + k, scatter=False, noisy=True)[1])  # observed soluble titre only — excludes lot scatter noise for cleaner OAT attribution
            k += 1
        po_slices.append(np.asarray(po_j))
        lf_slices.append(np.asarray(lf_j))
        rsm_slices.append(np.asarray(rsm_j))
    return po_slices, lf_slices, rsm_slices


def write_csv(path, x_opt, po_slices, lf_slices, rsm_slices, n_grid):
    lines = []
    for j in range(5):
        vals = np.linspace(BOUNDS[j, 0], BOUNDS[j, 1], int(n_grid))
        for ii, v in enumerate(vals):
            lines.append(','.join([str(j), VAR_LABELS[j], f'{v:.8g}', f'{po_slices[j][ii]:.8g}', f'{lf_slices[j][ii]:.8g}', f'{rsm_slices[j][ii]:.8g}']))
    with open(path, 'w') as fp:
        fp.write('dim,variable,x_axis_value,Po_observed,P_low_fidelity,P_rsm\n')
        fp.write('\n'.join(lines) + '\n')


def figures(po_slices, lf_slices, rsm_slices, n_grid, out_box, out_lines):
    fig1, ax1 = plt.subplots(figsize=(9, 4.5))
    ax1.boxplot(po_slices, labels=VAR_LABELS, notch=False, showmeans=True)  # one box per variable — spread is grid-driven + measurement noise only
    ax1.set_ylabel('Soluble proinsulin titre Po (observed, g/L)')
    ax1.set_title('OAT sensitivity: distribution of noisy Po along each dimension (oracle-opt fixed elsewhere)')
    fig1.tight_layout()
    fig1.savefig(out_box, dpi=140, bbox_inches='tight')
    plt.close(fig1)

    fig2, axes = plt.subplots(2, 3, figsize=(12, 7))
    axes = axes.ravel()
    for j in range(5):
        ax = axes[j]
        vals = np.linspace(BOUNDS[j, 0], BOUNDS[j, 1], int(n_grid))
        ax.plot(vals, lf_slices[j], color='C0', lw=1.8, label='ODE surrogate (biased LF)')  # mechanisticcheap curve
        ax.plot(vals, rsm_slices[j], '--', color='C1', lw=1.5, label='RSM quadratic (lstsq)')  # polynomial cross-check vs LF surface
        ax.scatter(vals, po_slices[j], s=14, alpha=0.55, color='k', zorder=0, label='Po observed')  # measurement scatter on simulated assay
        ax.set_title(VAR_LABELS[j])
        ax.set_xlabel(VAR_LABELS[j])
        ax.set_ylabel('titre g/L')
        if j == 0:
            ax.legend(fontsize=7, loc='best')
    axes[5].axis('off')
    fig2.suptitle('Surrogate verification on OAT grids (LF vs RSM vs observed Po)')
    fig2.tight_layout()
    fig2.savefig(out_lines, dpi=140, bbox_inches='tight')
    plt.close(fig2)


def print_verification(po_slices, lf_slices, rsm_slices):
    print('\nPearson r (LF vs RSM) along each OAT sweep:')
    for j in range(5):
        r_pv, _ = pearsonr(lf_slices[j], rsm_slices[j])
        rmse = float(np.sqrt(np.mean((lf_slices[j] - rsm_slices[j]) ** 2)))
        print(f"  {VAR_LABELS[j]:8s}  r={r_pv:+.4f}  RMSE(LF-RSM)={rmse:.4f}")
    print('\nPearson r (LF vs Po) - observed noise dilutes correlation:')
    for j in range(5):
        r_pv, _ = pearsonr(lf_slices[j], po_slices[j])
        print(f"  {VAR_LABELS[j]:8s}  r={r_pv:+.4f}")


def main():
    pa = argparse.ArgumentParser()
    pa.add_argument('--oracle-seed', type=int, default=0, help='RNG seed passed to differential_evolution in oracle DE — matches method_de_sim')
    pa.add_argument('--oat-seed', type=int, default=70_000, help='base seed ladder for noisy Po across OAT evaluations')
    pa.add_argument('--grid', type=int, default=35, help='uniform points per axis (bounds endpoints included)')
    pa.add_argument('--rsm-target', choices=('lf', 'wetlab'), default='lf', help='BBD fitting target: LF (default) vs noisy Po (Method 2 style)')
    pa.add_argument('--rsm-po-seed', type=int, default=9000, help='starting seed index for wetlab BBd training Po when --rsm-target wetlab')
    args = pa.parse_args()

    r_or = method_de_sim(seed=args.oracle_seed, scatter=False)  # scatter off — oracle trajectory uses nominal parameters only per protocol §5.3
    x_opt = np.asarray(r_or['best_x'], dtype=float)

    coefs = fit_global_quadratic_bbd(args.rsm_target, args.rsm_po_seed)
    po_sl, lf_sl, rsm_sl = sweep_one_at_a_time(x_opt, coefs, args.grid, args.oat_seed)

    csv_path = 'runs/sensitivity_surrogate_slices.csv'
    write_csv(csv_path, x_opt, po_sl, lf_sl, rsm_sl, args.grid)
    print(f'oracle x_opt = {np.round(x_opt, 6)}')
    print(f'wrote {csv_path}')

    figures(po_sl, lf_sl, rsm_sl, args.grid, 'plots/sensitivity_boxplot_observed.png', 'plots/sensitivity_surrogate_verify.png')
    print('wrote plots/sensitivity_boxplot_observed.png')
    print('wrote plots/sensitivity_surrogate_verify.png')
    print_verification(po_sl, lf_sl, rsm_sl)


if __name__ == '__main__':
    main()
