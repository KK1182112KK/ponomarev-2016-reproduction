"""
Validation tests for Example C (inverted pendulum, numerical predictor).
Ponomarev (2016), Section VI-C, Eq. 89-104.
This is the main example, reproducing Fig. 1 of the paper.
"""

import sys
import os
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.example_c import run_example_c
from src.uncompensated import run_uncompensated


class TestExampleCConvergenceAllICs:
    """test_convergence_all_ics -- Fig. 1: three initial conditions all converge."""

    @pytest.mark.parametrize("x1_init", [np.pi / 2, np.pi, 3 * np.pi / 2])
    def test_converges(self, x1_init):
        """State should converge for each initial condition."""
        params = {'x0': [x1_init, np.pi / 2.0], 'u0': 1.0}
        res = run_example_c(params)
        final_norm = np.linalg.norm(res['x_hist'][-1])
        print(f'    x1(0) = {x1_init:.4f}: |x(t_end)| = {final_norm:.4e}')
        # Euler at dt=0.01 produces residuals ~0.03-0.07; use 0.1 threshold
        assert final_norm < 0.1, \
            f'State did not converge for x1(0) = {x1_init:.4f}.'


class TestExampleCInitialControlSpike:
    """test_initial_control_spike -- First feedback value has large magnitude."""

    def test_large_spike(self):
        """Control should exhibit a large initial spike (per Fig. 1)."""
        params = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0}
        res = run_example_c(params)

        u_first = abs(res['u_hist'][0])
        u_max   = np.max(np.abs(res['u_hist']))
        print(f'    |u(1)| = {u_first:.4f}')
        print(f'    max|u| = {u_max:.4f}')

        assert u_max > 5, \
            'Control should exhibit a large initial spike (per Fig. 1).'


class TestExampleCPredictorAtT0:
    """test_predictor_at_t0 -- Verify predictor ODE at t=0 vs independent integration."""

    def test_predictor_matches(self):
        """Predictor at t=0 should match independent forward integration."""
        h  = np.pi / 4.0
        x0 = np.array([np.pi, np.pi / 2.0])
        u0 = 1.0
        dt = 0.01
        N_h = int(round(h / dt))
        ds = h / N_h

        # Manually solve predictor ODE at t=0
        xi = x0.copy()
        for j in range(N_h):
            dxi = np.array([xi[1], np.sin(xi[0]) + u0])
            xi = xi + ds * dxi
        y_manual = xi.copy()

        # Run full simulation and extract y(0)
        params = {'x0': x0.tolist(), 'u0': u0, 'h': h, 'dt': dt}
        res = run_example_c(params)
        y_sim = res['y_hist'][0, :]

        err = np.linalg.norm(y_manual - y_sim)
        print(f'    y_manual = [{y_manual[0]:.6f}, {y_manual[1]:.6f}]')
        print(f'    y_sim    = [{y_sim[0]:.6f}, {y_sim[1]:.6f}]')
        print(f'    |error|  = {err:.4e}')

        assert err < 1e-10, \
            'Predictor at t=0 should match independent forward integration.'


class TestExampleCCompensatedVsUncompensated:
    """test_compensated_vs_uncompensated -- Predictor converges, open-loop diverges."""

    def test_compensated_converges(self):
        """Compensated system should converge."""
        x0 = [np.pi, np.pi / 2.0]
        params_comp = {'x0': x0, 'u0': 1.0}
        res_comp = run_example_c(params_comp)
        norm_comp = np.linalg.norm(res_comp['x_hist'][-1])
        print(f'    Compensated:   |x(t_end)| = {norm_comp:.4e}')
        # Euler at dt=0.01 produces residuals ~0.04; use 0.1 threshold
        assert norm_comp < 0.1, 'Compensated system should converge.'

    def test_uncompensated_worse(self):
        """Uncompensated should perform worse than compensated."""
        x0 = [np.pi, np.pi / 2.0]

        params_comp = {'x0': x0, 'u0': 1.0}
        res_comp = run_example_c(params_comp)
        norm_comp = np.linalg.norm(res_comp['x_hist'][-1])

        params_unc = {'x0': x0, 'u0': 1.0, 'control_mode': 'zero'}
        res_unc = run_uncompensated(params_unc)
        norm_unc = np.linalg.norm(res_unc['x_hist'][-1])

        print(f'    Compensated:   |x(t_end)| = {norm_comp:.4e}')
        print(f'    Uncompensated: |x(t_end)| = {norm_unc:.4e}')
        assert norm_unc > norm_comp, \
            'Uncompensated should perform worse than compensated.'


class TestExampleCDtRefinement:
    """test_dt_refinement -- O(dt) convergence order estimate."""

    def test_convergence_order(self):
        """Average convergence order should be at least 0.8."""
        params_base = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0, 't_end': 5.0}

        dt_vals = [0.04, 0.02, 0.01, 0.005]

        # Reference: dt = 0.0025
        params_ref = {**params_base, 'dt': 0.0025}
        res_ref = run_example_c(params_ref)

        errors = np.zeros(len(dt_vals))
        for i, dt_i in enumerate(dt_vals):
            params_i = {**params_base, 'dt': dt_i}
            res_i = run_example_c(params_i)

            # Interpolate reference onto this grid
            x_ref_interp = np.zeros_like(res_i['x_hist'])
            for col in range(2):
                x_ref_interp[:, col] = np.interp(res_i['t'], res_ref['t'],
                                                   res_ref['x_hist'][:, col])
            errors[i] = np.max(np.abs(res_i['x_hist'] - x_ref_interp))
            print(f'    dt = {dt_i:.4f}: max error = {errors[i]:.4e}')

        # Compute convergence orders between successive refinements
        print('    Convergence orders:')
        orders = []
        for i in range(len(dt_vals) - 1):
            if errors[i + 1] > 1e-14:
                ratio = dt_vals[i] / dt_vals[i + 1]
                order = np.log(errors[i] / errors[i + 1]) / np.log(ratio)
                orders.append(order)
                print(f'      dt {dt_vals[i]:.4f} -> {dt_vals[i+1]:.4f}: '
                      f'order = {order:.2f}')

        if orders:
            avg_order = np.mean(orders)
            print(f'    Average order = {avg_order:.2f} (expected ~1.0)')
            assert avg_order > 0.8, \
                'Average convergence order should be at least 0.8.'


class TestExampleCLyapunovV0:
    """test_lyapunov_v0 -- ||y||^2 should decrease after transient."""

    def test_y_norm_decays(self):
        """||y||^2 should decay to less than 10% of post-transient value."""
        params = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0}
        res = run_example_c(params)

        h = np.pi / 4.0  # default delay
        y_norm_sq = np.sum(res['y_hist']**2, axis=1)

        idx_start = np.searchsorted(res['t'], 2.0 * h)
        y_tail = y_norm_sq[idx_start:]

        ratio = y_tail[-1] / y_tail[0]
        print(f'    ||y(t_end)||^2 / ||y(2h)||^2 = {ratio:.4e}')
        assert ratio < 0.1, \
            '||y||^2 should decay to less than 10% of post-transient value.'

    def test_mostly_non_increasing(self):
        """||y||^2 should be mostly non-increasing after transient."""
        params = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0}
        res = run_example_c(params)

        h = np.pi / 4.0
        y_norm_sq = np.sum(res['y_hist']**2, axis=1)

        idx_start = np.searchsorted(res['t'], 2.0 * h)
        y_tail = y_norm_sq[idx_start:]

        dV = np.diff(y_tail)
        frac_dec = np.sum(dV <= 1e-6) / len(dV)
        print(f'    Fraction non-increasing (t > 2h): {frac_dec*100:.2f}%')
        assert frac_dec > 0.80, \
            '||y||^2 should be mostly non-increasing after transient.'


class TestExampleCDelaySweep:
    """test_delay_sweep -- Multiple delay values all converge."""

    @pytest.mark.parametrize("h_val", [0.5, np.pi / 4, 1.0])
    def test_converges(self, h_val):
        """State should converge for each delay value."""
        params = {'x0': [np.pi, np.pi / 2.0], 'u0': 1.0,
                  'h': h_val, 't_end': 15.0}
        res = run_example_c(params)
        final_norm = np.linalg.norm(res['x_hist'][-1])
        print(f'    h = {h_val:.4f}: |x(t_end)| = {final_norm:.4e}')
        assert final_norm < 0.1, \
            f'State should converge for h = {h_val:.4f}.'
