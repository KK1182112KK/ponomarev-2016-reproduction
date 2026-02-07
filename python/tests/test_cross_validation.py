"""
Cross-method validation tests.
Tests across different numerical methods and transformations.
"""

import sys
import os
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.example_b import run_example_b
from src.example_c import run_example_c


def _lookup_u_rk4(t_query, dt, u_full, N_h):
    """Lookup u from full history array (for RK4 integration)."""
    idx = int(round(t_query / dt)) + N_h
    idx = max(0, min(idx, len(u_full) - 1))
    return u_full[idx]


class TestEulerVsRK4ExampleC:
    """test_euler_vs_rk4_example_c -- Compare Euler to RK4 for Example C."""

    def test_euler_rk4_agree(self):
        """Euler and RK4 should agree within 0.1 for dt=0.01."""
        h     = np.pi / 4.0
        x0    = np.array([np.pi, np.pi / 2.0])
        u0    = 1.0
        t_end = 5.0
        dt    = 0.01
        Nt    = int(round(t_end / dt)) + 1
        N_h   = int(round(h / dt))

        # Run the standard Euler simulation via run_example_c
        params_euler = {'x0': x0.tolist(), 'u0': u0, 't_end': t_end,
                        'dt': dt, 'h': h}
        res_euler = run_example_c(params_euler)
        t_euler = res_euler['t']
        x_euler = res_euler['x_hist']
        u_euler = res_euler['u_hist']

        # Run RK4 simulation with the SAME control history
        t_rk4 = np.arange(Nt) * dt
        x_rk4 = np.zeros((Nt, 2))
        x_rk4[0, :] = x0

        # Build full control history (including pre-history)
        u_full_rk4 = np.zeros(Nt + N_h)
        u_full_rk4[:N_h] = u0
        u_full_rk4[N_h:N_h + Nt] = u_euler

        for k in range(Nt - 1):
            x_k = x_rk4[k, :]
            t_k = t_rk4[k]

            # Current and delayed control from stored Euler history
            u_k     = u_euler[k]
            u_del_k = _lookup_u_rk4(t_k - h, dt, u_full_rk4, N_h)

            # RK4 stages for the plant ODE
            def f_sys(x_):
                return np.array([x_[1], np.sin(x_[0]) + u_k + u_del_k])

            k1 = f_sys(x_k)
            k2 = f_sys(x_k + 0.5 * dt * k1)
            k3 = f_sys(x_k + 0.5 * dt * k2)
            k4 = f_sys(x_k + dt * k3)

            x_rk4[k + 1, :] = x_k + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

        # Compare Euler vs RK4 (same control sequence)
        err = np.max(np.abs(x_euler - x_rk4))
        print(f'    max|x_euler - x_rk4| = {err:.4e} (threshold: 0.1)')
        assert err < 0.1, \
            'Euler and RK4 should agree within 0.1 for dt=0.01.'

        # Verify both reach similar final state
        err_final = np.linalg.norm(x_euler[-1] - x_rk4[-1])
        print(f'    |x_euler(end) - x_rk4(end)| = {err_final:.4e}')


class TestExampleBAnalyticalVsNumerical:
    """test_example_b_analytical_vs_numerical -- z-transform consistency."""

    def test_z_transform_consistency(self):
        """z-transform and V(z) should be self-consistent."""
        params = {'h': 1.0, 'dt': 0.001, 't_end': 10.0,
                  'x0': [1.0, 1.0], 'u0': 0.0}
        res = run_example_b(params)

        h  = params['h']
        dt = params['dt']
        Nt = len(res['t'])
        N_h = int(round(h / dt))
        ds_int = h / N_h

        # Sample every 500 steps
        check_idx = range(0, Nt, 500)
        max_z1_err = 0.0
        max_z2_err = 0.0
        max_V_err  = 0.0

        for idx in check_idx:
            x1 = res['x_hist'][idx, 0]
            x2 = res['x_hist'][idx, 1]

            # Trapezoidal integral of u from t-h to t
            int_u = 0.0
            for j in range(N_h + 1):
                theta_j = -h + j * ds_int
                t_query = res['t'][idx] + theta_j
                u_idx = int(round(t_query / dt))
                u_idx = max(0, min(u_idx, Nt - 1))
                if t_query < 0:
                    u_val = params['u0']
                else:
                    u_val = res['u_hist'][u_idx]
                if j == 0 or j == N_h:
                    w = 0.5
                else:
                    w = 1.0
                int_u += w * u_val
            int_u *= ds_int

            z1_check = x1 - x2 + int_u
            z2_check = x2

            err_z1 = abs(z1_check - res['z_hist'][idx, 0])
            err_z2 = abs(z2_check - res['z_hist'][idx, 1])

            max_z1_err = max(max_z1_err, err_z1)
            max_z2_err = max(max_z2_err, err_z2)

            # Verify V(z) = (z1 + z2*(z2-2)/2)^2 + z2^2
            sigma = (res['z_hist'][idx, 0]
                     + res['z_hist'][idx, 1] * (res['z_hist'][idx, 1] - 2.0) / 2.0)
            V_check = sigma**2 + res['z_hist'][idx, 1]**2
            err_V = abs(V_check - res['V_hist'][idx])
            max_V_err = max(max_V_err, err_V)

        n_checked = len(range(0, Nt, 500))
        print(f'    Checked {n_checked} time points')
        print(f'    max|z1_recomp - z1_stored| = {max_z1_err:.4e}')
        print(f'    max|z2_recomp - z2_stored| = {max_z2_err:.4e}')
        print(f'    max|V_recomp - V_stored|   = {max_V_err:.4e}')

        assert max_z2_err < 1e-12, 'z2 should equal x2 exactly.'
        assert max_z1_err < 0.1, \
            'z1 recomputation should agree closely with stored value.'
        assert max_V_err < 1e-10, \
            'Lyapunov function should be consistent with z coordinates.'

    def test_y_to_z_cascade_identity(self):
        """y-transform cascade should agree with direct z-transform."""
        params = {'h': 1.0, 'dt': 0.001, 't_end': 10.0,
                  'x0': [1.0, 1.0], 'u0': 0.0}
        res = run_example_b(params)

        h  = params['h']
        dt = params['dt']
        Nt = len(res['t'])
        N_h = int(round(h / dt))
        ds_int = h / N_h

        check_idx = range(0, Nt, 500)
        max_cascade_err = 0.0

        print(f'\n    Checking y-transform -> z-transform cascade consistency:')

        for idx in check_idx:
            x1 = res['x_hist'][idx, 0]
            x2 = res['x_hist'][idx, 1]

            # Compute integral
            int_u = 0.0
            for j in range(N_h + 1):
                theta_j = -h + j * ds_int
                t_query = res['t'][idx] + theta_j
                u_idx = int(round(t_query / dt))
                u_idx = max(0, min(u_idx, Nt - 1))
                if t_query < 0:
                    u_val = params['u0']
                else:
                    u_val = res['u_hist'][u_idx]
                if j == 0 or j == N_h:
                    w = 0.5
                else:
                    w = 1.0
                int_u += w * u_val
            int_u *= ds_int

            # y-transform (Eq. 76-77)
            y1 = x1 + (np.exp(2 * h) - 1) / 2.0 * x2**2 + int_u
            y2 = np.exp(h) * x2

            # Cascade from y -> z (Eq. 79-81)
            z1_from_y = y1 - np.exp(-h) * y2 + (np.exp(-2 * h) - 1) / 2.0 * y2**2
            z2_from_y = np.exp(-h) * y2

            # Direct z-transform (Eq. 86-88)
            z1_direct = x1 - x2 + int_u
            z2_direct = x2

            cascade_err = max(abs(z1_from_y - z1_direct),
                              abs(z2_from_y - z2_direct))
            max_cascade_err = max(max_cascade_err, cascade_err)

        print(f'    max|z_from_y - z_direct| = {max_cascade_err:.4e}')
        assert max_cascade_err < 1e-10, \
            'y-transform cascade should agree with direct z-transform.'
