"""
Validation tests for Example B (explicit prediction, cascade).
Ponomarev (2016), Section VI-B, Eq. 75-88.
"""

import sys
import os
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.example_b import run_example_b


class TestExampleBConvergence:
    """test_convergence -- Default params, state converges."""

    def test_state_converges(self):
        """State norm should converge."""
        res = run_example_b()
        final_norm = np.linalg.norm(res['x_hist'][-1])
        print(f'    |x(t_end)| = {final_norm:.4e} (threshold: 1e-2)')
        assert final_norm < 1e-2, 'State x did not converge to zero.'

    def test_lyapunov_small_at_end(self):
        """Lyapunov should be small at end."""
        res = run_example_b()
        print(f'    V(t_end) = {res["V_hist"][-1]:.4e}')
        assert res['V_hist'][-1] < 1e-3, \
            'Lyapunov function should be near zero at end.'

    def test_dimensions(self):
        """Dimensional checks."""
        res = run_example_b()
        assert len(res['t']) == res['x_hist'].shape[0]
        assert res['x_hist'].shape[1] == 2
        assert res['z_hist'].shape[1] == 2


class TestExampleBLyapunovDecay:
    """test_lyapunov_decay -- V(z) should be non-increasing after transient."""

    def test_mostly_decreasing(self):
        """V(z) should be nearly monotonically decreasing after transient."""
        params = {'h': 1.0}
        res = run_example_b(params)

        transient_time = 2.0 * params['h']
        idx_start = np.searchsorted(res['t'], transient_time)

        V_tail = res['V_hist'][idx_start:]
        dV = np.diff(V_tail)

        n_violations = np.sum(dV > 1e-10)
        frac_decreasing = 1.0 - n_violations / len(dV)
        print(f'    V values after t={transient_time:.1f}: {len(V_tail)} points')
        print(f'    Fraction non-increasing: {frac_decreasing:.4f}')

        assert frac_decreasing > 0.98, \
            'V(z) should be nearly monotonically decreasing after transient.'

    def test_significant_decay(self):
        """V should decrease significantly: V(end) < 0.1 * V(transient_end)."""
        params = {'h': 1.0}
        res = run_example_b(params)

        transient_time = 2.0 * params['h']
        idx_start = np.searchsorted(res['t'], transient_time)
        V_tail = res['V_hist'][idx_start:]

        V_ratio = res['V_hist'][-1] / V_tail[0]
        print(f'    V(t_end)/V(t_transient) = {V_ratio:.4e}')
        assert V_ratio < 0.1, \
            'V should decay to less than 10% of its post-transient value.'


class TestExampleBCascadeConvergence:
    """test_cascade_convergence -- z coordinates converge to zero."""

    def test_z_converges(self):
        """Cascade coordinates z should converge to zero."""
        res = run_example_b()
        final_z_norm = np.linalg.norm(res['z_hist'][-1])
        print(f'    |z(t_end)| = {final_z_norm:.4e} (threshold: 1e-2)')
        assert final_z_norm < 1e-2, \
            'Cascade coordinates z should converge to zero.'
        print(f'    z1(t_end) = {res["z_hist"][-1, 0]:.4e}, '
              f'z2(t_end) = {res["z_hist"][-1, 1]:.4e}')


class TestExampleBZTransformConsistency:
    """test_z_transform_consistency -- Verify z1 = x1 - x2 + int u at each step."""

    def test_z_recomputation(self):
        """z1 recomputation should match stored z1."""
        params = {'h': 1.0, 'dt': 0.001}
        res = run_example_b(params)

        h  = params['h']
        dt = params['dt']
        N_h = int(round(h / dt))
        Nt = len(res['t'])
        ds = h / N_h

        max_err_z1 = 0.0
        max_err_z2 = 0.0
        n_check = 0

        # Check at sampled time points (every 100 steps)
        check_indices = range(0, Nt, 100)
        for idx in check_indices:
            x1 = res['x_hist'][idx, 0]
            x2 = res['x_hist'][idx, 1]

            # Compute integral of u from t-h to t using trapezoidal rule
            int_u = 0.0
            for j in range(N_h + 1):
                theta_j = -h + j * ds
                t_query = res['t'][idx] + theta_j
                u_idx = int(round(t_query / dt))
                u_idx = max(0, min(u_idx, Nt - 1))
                if t_query < 0:
                    u_val = 0.0  # pre-history is 0
                else:
                    u_val = res['u_hist'][u_idx]
                if j == 0 or j == N_h:
                    w = 0.5
                else:
                    w = 1.0
                int_u += w * u_val
            int_u *= ds

            z1_recomp = x1 - x2 + int_u
            z2_recomp = x2

            err_z1 = abs(z1_recomp - res['z_hist'][idx, 0])
            err_z2 = abs(z2_recomp - res['z_hist'][idx, 1])

            max_err_z1 = max(max_err_z1, err_z1)
            max_err_z2 = max(max_err_z2, err_z2)
            n_check += 1

        print(f'    Checked {n_check} time points')
        print(f'    max|z1_recomp - z1_stored| = {max_err_z1:.4e}')
        print(f'    max|z2_recomp - z2_stored| = {max_err_z2:.4e}')

        # z2 = x2 should be exact
        assert max_err_z2 < 1e-12, 'z2 should equal x2 exactly.'

        # z1 recomputation should agree closely
        assert max_err_z1 < 0.1, 'z1 recomputation should match stored z1.'


class TestExampleBDtRefinement:
    """test_dt_refinement -- O(dt) convergence under step-size halving."""

    def test_convergence_order(self):
        """Convergence order should be approximately 1.0 (Euler)."""
        params_base = {'t_end': 5.0}

        # Coarse: dt = 0.002
        params_c = {**params_base, 'dt': 0.002}
        res_c = run_example_b(params_c)

        # Fine: dt = 0.001
        params_f = {**params_base, 'dt': 0.001}
        res_f = run_example_b(params_f)

        # Reference: dt = 0.0005
        params_r = {**params_base, 'dt': 0.0005}
        res_r = run_example_b(params_r)

        # Interpolate onto coarse grid
        x_f_interp = np.zeros_like(res_c['x_hist'])
        x_r_interp = np.zeros_like(res_c['x_hist'])
        for col in range(2):
            x_f_interp[:, col] = np.interp(res_c['t'], res_f['t'],
                                            res_f['x_hist'][:, col])
            x_r_interp[:, col] = np.interp(res_c['t'], res_r['t'],
                                            res_r['x_hist'][:, col])

        # Errors relative to reference
        err_coarse = np.max(np.abs(res_c['x_hist'] - x_r_interp))
        err_fine   = np.max(np.abs(x_f_interp - x_r_interp))

        print(f'    err(dt=0.002) = {err_coarse:.4e}')
        print(f'    err(dt=0.001) = {err_fine:.4e}')

        if err_fine > 1e-14:
            order = np.log2(err_coarse / err_fine)
            print(f'    Estimated order = {order:.2f} (expected ~1.0)')
            assert order > 0.7, 'Convergence order should be at least 0.7.'
        else:
            print('    Errors too small for order estimate.')
