#!/usr/bin/env python3
"""
One-click reproduction of Ponomarev (2016) results.

Usage:
    python run_all.py              # default 'sim' mode
    python run_all.py --mode sim   # run simulations only
    python run_all.py --mode fig   # generate figures only
    python run_all.py --mode test  # run validation tests only
    python run_all.py --mode all   # everything: sim + fig + test

Reference:
    Ponomarev (2016), "Nonlinear Predictor Feedback for Input-Affine
    Systems with Distributed Input Delays", IEEE TAC.
"""

import argparse
import sys
import os
import numpy as np

# Ensure the package root is on the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.example_a import run_example_a
from src.example_b import run_example_b
from src.example_c import run_example_c
from src.uncompensated import run_uncompensated
from src.figures import generate_figures


def print_banner():
    """Print identification banner."""
    print('=============================================')
    print('  Ponomarev (2016) -- Reproduction')
    print('  Nonlinear Predictor Feedback for')
    print('  Input-Affine Systems with Distributed')
    print('  Input Delays')
    print('=============================================')
    print()


def run_simulation():
    """Run all simulation examples."""
    print('--- Simulation ---\n')

    # Example A: Scalar case
    print('[Example A] Scalar predictor feedback (Eq. 73-74)')
    res_a = run_example_a()
    print(f'  Final: |x({res_a["t"][-1]:.0f})| = {abs(res_a["x_hist"][-1]):.4e}\n')

    # Example B: Explicit prediction
    print('[Example B] Explicit prediction (Eq. 75-88)')
    res_b = run_example_b()
    print(f'  Final: |x({res_b["t"][-1]:.0f})| = '
          f'{np.linalg.norm(res_b["x_hist"][-1]):.4e}, '
          f'V = {res_b["V_hist"][-1]:.4e}\n')

    # Example C: Numerical prediction (inverted pendulum)
    print('[Example C] Numerical prediction -- Inverted pendulum (Eq. 89-104)')
    x1_inits = [np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0]
    for x1_init in x1_inits:
        p = {'x0': [x1_init, np.pi / 2.0], 'u0': 1.0}
        res_c = run_example_c(p)
        print(f'  x1(0) = {x1_init:.4f}: |x({res_c["t"][-1]:.0f})| = '
              f'{np.linalg.norm(res_c["x_hist"][-1]):.4e}, '
              f'max|u| = {np.max(np.abs(res_c["u_hist"])):.2f}')

    # Cross-validation: compensated vs uncompensated
    print(f'\n[Comparison] Compensated vs Uncompensated (x0 = [pi, pi/2])')
    p_comp = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0}
    res_comp = run_example_c(p_comp)

    p_unc = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0, 'control_mode': 'zero'}
    res_unc = run_uncompensated(p_unc)

    p_naive = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0, 'control_mode': 'proportional'}
    res_naive = run_uncompensated(p_naive)

    print(f'  Predictor:     |x(t_end)| = '
          f'{np.linalg.norm(res_comp["x_hist"][-1]):.4e}')
    print(f'  Open-loop:     |x(t_end)| = '
          f'{np.linalg.norm(res_unc["x_hist"][-1]):.4e}')
    print(f'  Proportional:  |x(t_end)| = '
          f'{np.linalg.norm(res_naive["x_hist"][-1]):.4e}')


def run_figures():
    """Generate all figures."""
    print('--- Generating Figures ---\n')
    generate_figures(fig_set='all', dpi=300)


def run_tests():
    """Run validation tests via pytest."""
    print('--- Running Tests ---\n')
    import pytest
    test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tests')
    exit_code = pytest.main([test_dir, '-v', '--tb=short'])
    if exit_code != 0:
        print(f'\n{exit_code} test(s) failed.')
        sys.exit(exit_code)


def main():
    """Entry point."""
    parser = argparse.ArgumentParser(
        description='Ponomarev (2016) reproduction -- Python implementation')
    parser.add_argument('--mode', default='sim',
                        choices=['sim', 'fig', 'test', 'all'],
                        help='Execution mode (default: sim)')
    args = parser.parse_args()

    print_banner()

    if args.mode == 'sim':
        run_simulation()
    elif args.mode == 'fig':
        run_figures()
    elif args.mode == 'test':
        run_tests()
    elif args.mode == 'all':
        run_simulation()
        run_figures()
        run_tests()

    print('\n=== Done ===')


if __name__ == '__main__':
    main()
