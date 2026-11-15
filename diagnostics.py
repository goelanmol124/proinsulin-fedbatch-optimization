import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
from wetlab_sim import run_batch, simulate, BOUNDS
warnings.filterwarnings('ignore')
os.makedirs('plots', exist_ok=True)

def plot_trajectories():
    x_naive = (29.0, 0.1, 10.0, 0.1, 30.0)
    x_opt = (25.13, 0.3, 14.0, 0.043, 36.0)
    fig, axes = plt.subplots(3, 3, figsize=(11, 8), sharex=True)
    labels = ['V (L)', 'Xv (g/L)', 'Xd (g/L)', 'S (g/L)', 'A (g/L)', 'P_sol (g/L)', 'P_ib (g/L)', 'DO (mg/L)', 'T (°C)']
    for x, color, label in [(x_naive, 'C0', 'naive (T=29, IPTG=0.1, X_ind=10)'), (x_opt, 'C3', 'DE-sim opt (T=25, IPTG=0.30, X_ind=14)')]:
        tr = simulate(x, seed=0, scatter=False)
        for k in range(9):
            ax = axes[k // 3, k % 3]
            ax.plot(tr['t'], tr['y'][k], color=color, lw=1.5, label=label if k == 0 else None)
            ax.set_ylabel(labels[k], fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel('t (h)')
    axes[0, 0].legend(loc='upper left', fontsize=8, frameon=False, bbox_to_anchor=(0, 1.4))
    fig.suptitle('9-state ODE: naive baseline vs DE-sim optimum', y=1.02)
    fig.tight_layout()
    fig.savefig('plots/trajectory.png', dpi=120, bbox_inches='tight')
    plt.close(fig)
    print('wrote plots/trajectory.png')

def plot_convergence(data):
    from optimization import oracle_score
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = dict(trial_error='C0', doe_rsm='C1', gp_bo='C3', mfbo_ude='C4')
    for m in ['trial_error', 'doe_rsm', 'gp_bo', 'mfbo_ude']:
        hists = data[f'hist_{m}']
        curves = []
        max_len = max((len(h) for h in hists))
        for h in hists:
            xs = [pt[0] for pt in h]
            pos = [pt[1] for pt in h]
            curve = np.zeros(len(h))
            best_idx = 0
            best_po = -np.inf
            for k, po in enumerate(pos):
                if po > best_po:
                    best_po = po
                    best_idx = k
                curve[k] = oracle_score(xs[best_idx])
            curve = np.concatenate([curve, np.full(max_len - len(curve), curve[-1])])
            curves.append(curve)
        curves = np.array(curves)
        x = np.arange(1, max_len + 1)
        med = np.median(curves, axis=0)
        q1, q3 = np.percentile(curves, [25, 75], axis=0)
        ax.plot(x, med, color=colors[m], lw=2, label=m)
        ax.fill_between(x, q1, q3, color=colors[m], alpha=0.18)
    ax.set_xlabel('wet-lab batches consumed')
    ax.set_ylabel('oracle P_sol at running argmax-Po (g/L)')
    ax.set_title('Convergence — oracle-scored running recommendation')
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('plots/convergence.png', dpi=120)
    plt.close(fig)
    print('wrote plots/convergence.png')

def plot_comparison(arr):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    methods = ['trial_error', 'doe_rsm', 'gp_bo', 'mfbo_ude']
    data = [arr['P_true'][arr['method'] == m] for m in methods]
    bp = ax.boxplot(data, labels=methods, patch_artist=True, widths=0.5)
    for patch, c in zip(bp['boxes'], ['C0', 'C1', 'C3', 'C4']):
        patch.set_facecolor(c)
        patch.set_alpha(0.4)
    for i, d in enumerate(data, start=1):
        jitter = np.random.default_rng(i).normal(0, 0.04, size=len(d))
        ax.scatter(i + jitter, d, color='k', s=18, zorder=3)
    ax.set_ylabel('final P_sol_true (g/L)')
    ax.set_title('Final titre by method (one point per seed)')
    ax.grid(alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig('plots/comparison.png', dpi=120)
    plt.close(fig)
    print('wrote plots/comparison.png')

def plot_sensitivity():
    x_star = np.array([25.13, 0.3, 14.0, 0.043, 36.0])
    names = ['T_prod', 'IPTG', 'X_ind', 'F_s', 't_f']
    fig, axes = plt.subplots(1, 5, figsize=(15, 3.2), sharey=True)
    for k, ax in enumerate(axes):
        lo, hi = BOUNDS[k]
        xs = np.linspace(lo, hi, 25)
        ys = []
        for v in xs:
            x_q = x_star.copy()
            x_q[k] = v
            Pt, _, _ = run_batch(x_q, seed=0, scatter=False, noisy=False)
            ys.append(Pt)
        ax.plot(xs, ys, 'C2-')
        ax.axvline(x_star[k], color='k', linestyle='--', alpha=0.4)
        ax.set_xlabel(names[k])
        ax.grid(alpha=0.3)
    axes[0].set_ylabel('P_sol_true (g/L)')
    fig.suptitle('1-D sensitivity around DE-sim optimum (other vars held)', y=1.03)
    fig.tight_layout()
    fig.savefig('plots/sensitivity.png', dpi=120, bbox_inches='tight')
    plt.close(fig)
    print('wrote plots/sensitivity.png')

def plot_best_points(arr):
    names = ['T_prod', 'IPTG', 'X_ind', 'F_s', 't_f']
    methods = ['trial_error', 'doe_rsm', 'gp_bo', 'mfbo_ude']
    colors = ['C0', 'C1', 'C3', 'C4']
    fig, axes = plt.subplots(1, 5, figsize=(16, 3.0))
    for k, ax in enumerate(axes):
        for m, c in zip(methods, colors):
            sub = arr[arr['method'] == m]
            ax.scatter(sub[names[k]], sub['P_true'], color=c, label=m if k == 0 else None, s=40, alpha=0.7)
        ax.set_xlabel(names[k])
        ax.set_xlim(*BOUNDS[k])
        ax.grid(alpha=0.3)
    axes[0].set_ylabel('P_sol_true (g/L)')
    axes[0].legend(loc='upper left', fontsize=8, frameon=False, bbox_to_anchor=(0, 1.35), ncol=4)
    fig.suptitle('Best operating point per seed, by method — where do methods land?', y=1.05)
    fig.tight_layout()
    fig.savefig('plots/best_points.png', dpi=120, bbox_inches='tight')
    plt.close(fig)
    print('wrote plots/best_points.png')

def reproducibility(arr):
    print('\n=== Reproducibility ===')
    for m in ['trial_error', 'doe_rsm', 'gp_bo', 'mfbo_ude']:
        sub = arr[arr['method'] == m]
        ref = float(np.median(sub['P_true']))
        within = np.mean(np.abs(sub['P_true'] - ref) / max(ref, 1e-06) < 0.1)
        print(f'  {m:12s}: median={ref:.3f} g/L, fraction within ±10% = {within:.2f}')

def wilcoxon_test(arr):
    from scipy.stats import wilcoxon
    seeds = np.sort(np.unique(arr['seed']))
    champ = np.array([float(arr[(arr['method'] == 'mfbo_ude') & (arr['seed'] == s)]['P_true']) for s in seeds])
    print('\n=== Wilcoxon signed-rank (mfbo_ude > other) ===')
    for m in ['trial_error', 'doe_rsm', 'gp_bo']:
        other = np.array([float(arr[(arr['method'] == m) & (arr['seed'] == s)]['P_true']) for s in seeds])
        try:
            stat, p = wilcoxon(champ, other, alternative='greater')
            print(f'  mfbo_ude vs {m:12s}: stat={stat:.2f}  p={p:.4f}  (n={len(seeds)})')
        except ValueError as e:
            print(f'  mfbo_ude vs {m:12s}: {e}')

def plot_pareto(arr):
    fig, ax = plt.subplots(figsize=(7, 5))
    methods = ['trial_error', 'doe_rsm', 'gp_bo', 'mfbo_ude']
    colors = ['C0', 'C1', 'C3', 'C4']
    for m, c in zip(methods, colors):
        sub = arr[arr['method'] == m]
        ax.scatter(sub['n_wetlab'], sub['P_true'], color=c, label=m, s=70, alpha=0.7, edgecolor='k', linewidth=0.5)
        ax.scatter(np.median(sub['n_wetlab']), np.median(sub['P_true']), color=c, marker='*', s=300, edgecolor='k', linewidth=1.0)
    ax.set_xlabel('wet-lab batches consumed')
    ax.set_ylabel('final P_sol_true (g/L)')
    ax.set_title('Pareto frontier: titre vs wet-lab cost (stars = medians)')
    ax.set_xscale('log')
    ax.grid(alpha=0.3, which='both')
    ax.legend(loc='lower left')
    fig.tight_layout()
    fig.savefig('plots/pareto.png', dpi=120)
    plt.close(fig)
    print('wrote plots/pareto.png')

def main():
    data = np.load('results.npz', allow_pickle=True)
    arr = data['table']
    plot_trajectories()
    plot_convergence(data)
    plot_comparison(arr)
    plot_sensitivity()
    plot_best_points(arr)
    plot_pareto(arr)
    reproducibility(arr)
    wilcoxon_test(arr)
if __name__ == '__main__':
    main()
