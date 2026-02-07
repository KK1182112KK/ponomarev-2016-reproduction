"""
Generate all figures for Ponomarev (2016) reproduction.

Figures:
    Fig 1: Paper's Fig. 1 -- u(t) and x1(t) for Example C, 3 initial conditions
    Fig 2: Example A -- scalar state and control trajectories
    Fig 3: Example B -- state, cascade coordinates, and Lyapunov function
    Fig 4: Example C -- compensated vs uncompensated comparison
    Fig 5: Example C -- delay sweep h in [0.1, pi/2]
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .example_a import run_example_a
from .example_b import run_example_b
from .example_c import run_example_c
from .uncompensated import run_uncompensated


def _get_results_dir():
    """Return the results directory, creating it if needed."""
    results_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'results')
    results_dir = os.path.abspath(results_dir)
    os.makedirs(results_dir, exist_ok=True)
    return results_dir


def make_fig1(results_dir=None, dpi=300):
    """Fig 1: Paper's Figure 1 (Example C, 3 initial conditions)."""
    if results_dir is None:
        results_dir = _get_results_dir()
    print('  Generating Fig 1 (Paper Fig. 1 -- Example C, 3 ICs)...')

    x1_inits = [np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0]
    colors = ['b', 'r', 'k']
    labels = [r'$x_1(0) = \pi/2$', r'$x_1(0) = \pi$', r'$x_1(0) = 3\pi/2$']

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))

    # Left: control u(t)
    ax = axes[0]
    for i in range(3):
        p = {'x0': [x1_inits[i], np.pi / 2.0], 'u0': 1.0}
        res = run_example_c(p)
        ax.plot(res['t'], res['u_hist'], colors[i], linewidth=1.2, label=labels[i])
    ax.set_xlabel('t')
    ax.set_ylabel('u(t)')
    ax.set_title('Control Input')
    ax.legend(loc='best')
    ax.grid(True)

    # Right: state x1(t)
    ax = axes[1]
    for i in range(3):
        p = {'x0': [x1_inits[i], np.pi / 2.0], 'u0': 1.0}
        res = run_example_c(p)
        ax.plot(res['t'], res['x_hist'][:, 0], colors[i], linewidth=1.2, label=labels[i])
    ax.set_xlabel('t')
    ax.set_ylabel(r'$x_1(t)$')
    ax.set_title('Angular Position')
    ax.legend(loc='best')
    ax.grid(True)

    fig.suptitle('Figure 1: Example C -- Predictor Feedback (Ponomarev 2016)')
    fig.tight_layout()

    filepath = os.path.join(results_dir, 'fig1_paper_figure1.png')
    fig.savefig(filepath, dpi=dpi)
    plt.close(fig)
    print(f'    Saved fig1_paper_figure1.png')


def make_fig2(results_dir=None, dpi=300):
    """Fig 2: Example A -- scalar state and control."""
    if results_dir is None:
        results_dir = _get_results_dir()
    print('  Generating Fig 2 (Example A -- scalar case)...')

    res = run_example_a()

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))

    axes[0].plot(res['t'], res['x_hist'], 'b', linewidth=1.5)
    axes[0].axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    axes[0].set_xlabel('t')
    axes[0].set_ylabel('x(t)')
    axes[0].set_title('State')
    axes[0].grid(True)

    axes[1].plot(res['t'], res['u_hist'], 'r', linewidth=1.5)
    axes[1].axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    axes[1].set_xlabel('t')
    axes[1].set_ylabel('u(t)')
    axes[1].set_title('Control Input')
    axes[1].grid(True)

    fig.suptitle('Figure 2: Example A -- Scalar Predictor Feedback')
    fig.tight_layout()

    filepath = os.path.join(results_dir, 'fig2_example_a.png')
    fig.savefig(filepath, dpi=dpi)
    plt.close(fig)
    print(f'    Saved fig2_example_a.png')


def make_fig3(results_dir=None, dpi=300):
    """Fig 3: Example B -- state trajectories and Lyapunov function."""
    if results_dir is None:
        results_dir = _get_results_dir()
    print('  Generating Fig 3 (Example B -- explicit prediction)...')

    res = run_example_b()

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    # States x1, x2
    axes[0].plot(res['t'], res['x_hist'][:, 0], 'b', linewidth=1.2, label=r'$x_1$')
    axes[0].plot(res['t'], res['x_hist'][:, 1], 'r', linewidth=1.2, label=r'$x_2$')
    axes[0].set_xlabel('t')
    axes[0].set_ylabel('x(t)')
    axes[0].set_title('State Trajectories')
    axes[0].legend(loc='best')
    axes[0].grid(True)

    # Cascade coordinates z1, z2
    axes[1].plot(res['t'], res['z_hist'][:, 0], 'b', linewidth=1.2, label=r'$z_1$')
    axes[1].plot(res['t'], res['z_hist'][:, 1], 'r', linewidth=1.2, label=r'$z_2$')
    axes[1].set_xlabel('t')
    axes[1].set_ylabel('z(t)')
    axes[1].set_title('Cascade Coordinates')
    axes[1].legend(loc='best')
    axes[1].grid(True)

    # Lyapunov function
    axes[2].semilogy(res['t'], np.maximum(res['V_hist'], 1e-16), 'k', linewidth=1.5)
    axes[2].set_xlabel('t')
    axes[2].set_ylabel('V(z)')
    axes[2].set_title('Lyapunov Function')
    axes[2].grid(True)

    fig.suptitle('Figure 3: Example B -- Explicit Prediction')
    fig.tight_layout()

    filepath = os.path.join(results_dir, 'fig3_example_b.png')
    fig.savefig(filepath, dpi=dpi)
    plt.close(fig)
    print(f'    Saved fig3_example_b.png')


def make_fig4(results_dir=None, dpi=300):
    """Fig 4: Example C -- compensated vs uncompensated comparison."""
    if results_dir is None:
        results_dir = _get_results_dir()
    print('  Generating Fig 4 (Example C -- compensated vs uncompensated)...')

    # Compensated (predictor feedback)
    p_comp = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0}
    res_c = run_example_c(p_comp)

    # Uncompensated (open-loop, u=0)
    p_unc = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0, 'control_mode': 'zero'}
    res_u = run_uncompensated(p_unc)

    # Uncompensated (naive proportional)
    p_naive = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0, 'control_mode': 'proportional'}
    res_n = run_uncompensated(p_naive)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    # x1(t)
    axes[0].plot(res_c['t'], res_c['x_hist'][:, 0], 'b', linewidth=1.5, label='Predictor')
    axes[0].plot(res_u['t'], res_u['x_hist'][:, 0], 'r--', linewidth=1.2, label='Open-loop')
    axes[0].plot(res_n['t'], res_n['x_hist'][:, 0], 'g-.', linewidth=1.2, label='Proportional')
    axes[0].set_xlabel('t')
    axes[0].set_ylabel(r'$x_1(t)$')
    axes[0].set_title('Angular Position')
    axes[0].legend(loc='best')
    axes[0].grid(True)

    # x2(t)
    axes[1].plot(res_c['t'], res_c['x_hist'][:, 1], 'b', linewidth=1.5, label='Predictor')
    axes[1].plot(res_u['t'], res_u['x_hist'][:, 1], 'r--', linewidth=1.2, label='Open-loop')
    axes[1].plot(res_n['t'], res_n['x_hist'][:, 1], 'g-.', linewidth=1.2, label='Proportional')
    axes[1].set_xlabel('t')
    axes[1].set_ylabel(r'$x_2(t)$')
    axes[1].set_title('Angular Velocity')
    axes[1].legend(loc='best')
    axes[1].grid(True)

    # u(t)
    axes[2].plot(res_c['t'], res_c['u_hist'], 'b', linewidth=1.5, label='Predictor')
    axes[2].plot(res_n['t'], res_n['u_hist'], 'g-.', linewidth=1.2, label='Proportional')
    axes[2].set_xlabel('t')
    axes[2].set_ylabel('u(t)')
    axes[2].set_title('Control Input')
    axes[2].legend(loc='best')
    axes[2].grid(True)

    fig.suptitle('Figure 4: Example C -- Compensated vs Uncompensated')
    fig.tight_layout()

    filepath = os.path.join(results_dir, 'fig4_compensated_comparison.png')
    fig.savefig(filepath, dpi=dpi)
    plt.close(fig)
    print(f'    Saved fig4_compensated_comparison.png')


def make_fig5(results_dir=None, dpi=300):
    """Fig 5: Example C -- delay sweep h in [0.1, pi/2]."""
    if results_dir is None:
        results_dir = _get_results_dir()
    print('  Generating Fig 5 (Example C -- delay sweep)...')

    h_vals = np.linspace(0.1, np.pi / 2.0, 8)
    final_norms = np.zeros(len(h_vals))

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))

    cmap = plt.cm.get_cmap('tab10', len(h_vals))

    ax = axes[0]
    for i, h_val in enumerate(h_vals):
        p = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0, 'h': h_val}
        res = run_example_c(p)
        final_norms[i] = np.linalg.norm(res['x_hist'][-1])
        ax.plot(res['t'], res['x_hist'][:, 0], color=cmap(i), linewidth=1.0,
                label=f'h={h_val:.2f}')
    ax.set_xlabel('t')
    ax.set_ylabel(r'$x_1(t)$')
    ax.set_title(r'State $x_1$ for various $h$')
    ax.legend(loc='best', fontsize=7)
    ax.grid(True)

    ax = axes[1]
    ax.bar(np.arange(len(h_vals)), np.log10(np.maximum(final_norms, 1e-16)),
           tick_label=[f'{h_val:.2f}' for h_val in h_vals])
    ax.set_xlabel('Delay h')
    ax.set_ylabel(r'$\log_{10}(|x(t_{end})|)$')
    ax.set_title('Final State Norm vs Delay')
    ax.grid(True)

    fig.suptitle('Figure 5: Example C -- Delay Sweep')
    fig.tight_layout()

    filepath = os.path.join(results_dir, 'fig5_delay_sweep.png')
    fig.savefig(filepath, dpi=dpi)
    plt.close(fig)
    print(f'    Saved fig5_delay_sweep.png')


def generate_figures(fig_set='all', dpi=300):
    """Generate figures.

    Parameters
    ----------
    fig_set : str
        'all' (default), or specific figure: 'fig1','fig2','fig3','fig4','fig5'
    dpi : int
        Resolution for saving (default 300)
    """
    results_dir = _get_results_dir()

    dispatch = {
        'fig1': make_fig1,
        'fig2': make_fig2,
        'fig3': make_fig3,
        'fig4': make_fig4,
        'fig5': make_fig5,
    }

    fig_set = fig_set.lower()
    if fig_set == 'all':
        for func in dispatch.values():
            func(results_dir, dpi)
    elif fig_set in dispatch:
        dispatch[fig_set](results_dir, dpi)
    else:
        raise ValueError(f'Unknown figure set: {fig_set}')
