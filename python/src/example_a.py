"""
Scalar predictor feedback (Ponomarev 2016, Eq. 73-74).

System:
    dx/dt = f(x) + b0*u(t) + b1*u(t-h) + int_{-h}^{0} bint*u(t+theta) dtheta

Predictor transformed system:
    dy/dt = f(y) + B(y, phi)*u(t)

Feedback (Eq. 74):
    kappa = -(f(y) + y) / B(y, phi)
"""

import numpy as np


def _lookup_u(t_query, dt, u_full, N_h):
    """Lookup u from full history array (including pre-history).

    u_full is indexed so that u_full[N_h] = u(0),
    u_full[0] = u(-h), etc.
    Index for time t_query: idx = round(t_query / dt) + N_h
    """
    idx = int(round(t_query / dt)) + N_h
    idx = max(0, min(idx, len(u_full) - 1))
    return u_full[idx]


def run_example_a(params=None):
    """Run Example A: scalar predictor feedback.

    Parameters
    ----------
    params : dict, optional
        Override default parameters. Keys:
            f     : callable, nonlinearity (default: np.sin)
            b0    : float, direct input gain (default: 1.0)
            b1    : float, delayed input gain (default: 0.5)
            bint  : float, distributed delay kernel (default: 1.0)
            h     : float, delay (default: 0.5)
            x0    : float, initial state (default: 1.0)
            u0    : float, initial control history value (default: 0.0)
            t_end : float, simulation end time (default: 10.0)
            dt    : float, time step (default: 0.001)

    Returns
    -------
    dict with keys: t, x_hist, u_hist
    """
    if params is None:
        params = {}

    f_func = params.get('f', np.sin)
    b0     = params.get('b0', 1.0)
    b1     = params.get('b1', 0.5)
    bint   = params.get('bint', 1.0)
    h      = params.get('h', 0.5)
    x0     = params.get('x0', 1.0)
    u0     = params.get('u0', 0.0)
    t_end  = params.get('t_end', 10.0)
    dt     = params.get('dt', 0.001)

    # Derived quantities
    Nt  = int(round(t_end / dt)) + 1
    N_h = int(round(h / dt))

    # Pre-allocate
    t      = np.arange(Nt) * dt
    x_hist = np.zeros(Nt)
    u_full = np.zeros(Nt + N_h)

    # Initial conditions
    x_hist[0] = x0
    u_full[:N_h] = u0

    # Jacobian of f (scalar): df/dx via finite differences
    def df(x):
        eps = 1e-7
        return (f_func(x + eps) - f_func(x - eps)) / (2.0 * eps)

    # Main time-stepping loop
    ds = h / N_h if N_h > 0 else h

    for k in range(Nt - 1):
        x_k = x_hist[k]

        # --- Solve predictor ODE and beta ODE over s in [0, h] ---
        xi   = x_k
        beta = b0

        # Pre-compute distributed delay integral at s=0:
        # I(0) = int_{-h}^{0} bint * u(t + alpha) d_alpha
        # Trapezoidal rule over stored history
        int_pred = 0.0
        for jj in range(N_h + 1):
            t_q = t[k] - h + jj * ds
            u_q = _lookup_u(t_q, dt, u_full, N_h)
            if jj == 0 or jj == N_h:
                ww = 0.5
            else:
                ww = 1.0
            int_pred += ww * bint * u_q
        int_pred *= ds

        for j in range(N_h):
            s_j = j * ds

            # u(t + s - h)
            t_query_pred = t[k] + s_j - h
            u_pred = _lookup_u(t_query_pred, dt, u_full, N_h)

            # A(xi) for scalar = df/dx evaluated at xi
            A_xi = df(xi)

            # Euler step for predictor
            xi = xi + ds * (f_func(xi) + b1 * u_pred + int_pred)

            # Euler step for beta
            beta = beta + ds * (A_xi * beta + bint)

            # Incremental update: I(j+1) = I(j) - ds * bint * u(t+s_j-h)
            int_pred = int_pred - ds * bint * u_pred

        y = xi                  # predicted state
        B = b1 + beta           # effective input gain

        # --- Feedback (Eq. 74) ---
        if abs(B) < 1e-10:
            u_k = 0.0           # safety: avoid division by zero
        else:
            u_k = -(f_func(y) + y) / B

        # Store control
        u_full[N_h + k] = u_k

        # --- Actual system dynamics ---
        # dx/dt = f(x) + b0*u(t) + b1*u(t-h) + int_{-h}^{0} bint*u(t+theta) dtheta

        # u(t-h)
        u_delayed = _lookup_u(t[k] - h, dt, u_full, N_h)

        # Distributed delay integral via trapezoidal rule
        int_val = 0.0
        for j in range(N_h + 1):
            theta_j = -h + j * ds
            t_query_int = t[k] + theta_j
            u_j = _lookup_u(t_query_int, dt, u_full, N_h)
            if j == 0 or j == N_h:
                w = 0.5
            else:
                w = 1.0
            int_val += w * bint * u_j
        int_val *= ds

        dxdt = f_func(x_k) + b0 * u_k + b1 * u_delayed + int_val

        # Euler step
        x_hist[k + 1] = x_k + dt * dxdt

    # Final control (hold last value)
    u_full[N_h + Nt - 1] = u_full[N_h + Nt - 2]

    # Extract control history aligned with t
    u_hist = u_full[N_h:N_h + Nt].copy()

    print(f'  Example A: |x(t_end)| = {abs(x_hist[-1]):.4e}')

    return {'t': t, 'x_hist': x_hist, 'u_hist': u_hist}
