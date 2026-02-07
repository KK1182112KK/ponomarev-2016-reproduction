"""
Validation tests for Example A (scalar predictor feedback).
Ponomarev (2016), Section VI-A, Eq. 73-74.
"""

import sys
import os
import numpy as np
import pytest

# Ensure src is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.example_a import run_example_a


class TestExampleAConvergence:
    """test_convergence -- Run with default params, verify convergence."""

    def test_state_converges(self):
        """State should converge to zero."""
        res = run_example_a()
        final_x = abs(res['x_hist'][-1])
        print(f'    |x(t_end)| = {final_x:.4e} (threshold: 1e-2)')
        assert final_x < 1e-2, 'State x did not converge to zero within tolerance.'

    def test_control_settles(self):
        """Control should also settle near zero."""
        res = run_example_a()
        final_u = abs(res['u_hist'][-1])
        print(f'    |u(t_end)| = {final_u:.4e}')
        assert final_u < 1.0, 'Control u did not settle near zero.'

    def test_dimensions(self):
        """Basic dimensional checks."""
        res = run_example_a()
        assert len(res['t']) == len(res['x_hist']), \
            'Dimension mismatch: t and x_hist.'
        assert len(res['t']) == len(res['u_hist']), \
            'Dimension mismatch: t and u_hist.'


class TestExampleALinearPredictor:
    """test_linear_predictor -- f(x)=0, system is purely linear."""

    def test_tight_convergence(self):
        """With f(x)=0, predictor should give exact prediction."""
        params = {
            'f': lambda x: 0.0,
            'b0': 1.0, 'b1': 0.5, 'bint': 1.0,
            'h': 0.5, 'x0': 1.0, 'u0': 0.0,
            't_end': 10.0, 'dt': 0.001,
        }
        res = run_example_a(params)
        final_x = abs(res['x_hist'][-1])
        print(f'    |x(t_end)| = {final_x:.4e} (threshold: 1e-4)')
        assert final_x < 1e-4, 'Linear case should converge more tightly.'

    def test_monotone_decay_tail(self):
        """After transient, state should decay roughly monotonically."""
        params = {
            'f': lambda x: 0.0,
            'b0': 1.0, 'b1': 0.5, 'bint': 1.0,
            'h': 0.5, 'x0': 1.0, 'u0': 0.0,
            't_end': 10.0, 'dt': 0.001,
        }
        res = run_example_a(params)
        t = res['t']
        x = res['x_hist']

        idx_start = np.searchsorted(t, 2.0 * params['h'])
        x_tail = np.abs(x[idx_start:])
        diffs = np.diff(x_tail)
        frac_decreasing = np.sum(diffs <= 1e-10) / len(diffs)
        print(f'    Fraction decreasing (t > 2h): {frac_decreasing*100:.2f}%')
        assert frac_decreasing > 0.90, \
            'Linear case tail should be mostly monotonically decreasing.'


class TestExampleAExponentialDecay:
    """test_exponential_decay -- f(x)=x (linear f), qualitative convergence."""

    def test_converges_with_unstable_f(self):
        """Should converge despite unstable open-loop dynamics."""
        params = {
            'f': lambda x: x,
            'b0': 1.0, 'b1': 0.5, 'bint': 1.0,
            'h': 0.5, 'x0': 1.0, 'u0': 0.0,
            't_end': 15.0, 'dt': 0.001,
        }
        res = run_example_a(params)
        final_x = abs(res['x_hist'][-1])
        print(f'    |x(t_end)| = {final_x:.4e} (threshold: 1e-2)')
        assert final_x < 1e-2, \
            'f(x)=x case should converge with predictor feedback.'

    def test_bounded(self):
        """State should not blow up."""
        params = {
            'f': lambda x: x,
            'b0': 1.0, 'b1': 0.5, 'bint': 1.0,
            'h': 0.5, 'x0': 1.0, 'u0': 0.0,
            't_end': 15.0, 'dt': 0.001,
        }
        res = run_example_a(params)
        max_x = np.max(np.abs(res['x_hist']))
        print(f'    max|x| = {max_x:.4e}')
        assert max_x < 100, 'State should remain bounded.'


class TestExampleADtRefinement:
    """test_dt_refinement -- O(dt) Euler convergence."""

    def test_convergence_order(self):
        """Convergence order should be approximately 1.0 (Euler)."""
        params_base = {'t_end': 5.0}

        # Coarse run: dt = 0.002
        params_coarse = {**params_base, 'dt': 0.002}
        res_c = run_example_a(params_coarse)

        # Fine run: dt = 0.001
        params_fine = {**params_base, 'dt': 0.001}
        res_f = run_example_a(params_fine)

        # Reference run: dt = 0.0005
        params_ref = {**params_base, 'dt': 0.0005}
        res_r = run_example_a(params_ref)

        # Interpolate all onto the coarse time grid for comparison
        x_f_interp = np.interp(res_c['t'], res_f['t'], res_f['x_hist'])
        x_r_interp = np.interp(res_c['t'], res_r['t'], res_r['x_hist'])

        # Errors relative to reference
        err_coarse = np.max(np.abs(res_c['x_hist'] - x_r_interp))
        err_fine   = np.max(np.abs(x_f_interp - x_r_interp))

        print(f'    err(dt=0.002) = {err_coarse:.4e}')
        print(f'    err(dt=0.001) = {err_fine:.4e}')

        if err_fine > 1e-14:
            order = np.log2(err_coarse / err_fine)
            print(f'    Estimated order = {order:.2f} (expected ~1.0)')
            assert order > 0.7, \
                'Convergence order should be at least 0.7 (close to O(dt)).'
        else:
            print('    Error too small to estimate order; both are near reference.')
