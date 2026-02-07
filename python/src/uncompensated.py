"""
Example C without predictor feedback (no delay compensation).

System (same as Example C):
    dx1/dt = x2
    dx2/dt = sin(x1) + u(t) + u(t-h),   h = pi/4

Control: Simple proportional feedback WITHOUT predictor.
Three modes available via params['control_mode']:
    'zero'         -- u(t) = 0 (open-loop, default)
    'proportional' -- u(t) = -k1*x1 - k2*x2 (delay-ignorant)
    'naive_pred'   -- u(t) = -x1 - x2 (as if no delay)
"""

import numpy as np


def _lookup_u(t_query, dt, u_full, N_h):
    """Lookup u from full history array."""
    idx = int(round(t_query / dt)) + N_h
    idx = max(0, min(idx, len(u_full) - 1))
    return u_full[idx]


def run_uncompensated(params=None):
    """Run uncompensated baseline for the inverted pendulum.

    Parameters
    ----------
    params : dict, optional
        Override default parameters. Keys:
            h            : float, delay (default: pi/4)
            x0           : array-like, initial state [x1, x2] (default: [pi, pi/2])
            u0           : float, initial control history value (default: 1.0)
            t_end        : float, simulation end time (default: 10.0)
            dt           : float, time step (default: 0.01)
            control_mode : str, 'zero', 'proportional', or 'naive_pred' (default: 'zero')
            k1           : float, proportional gain 1 (default: 1.0)
            k2           : float, proportional gain 2 (default: 1.0)

    Returns
    -------
    dict with keys: t, x_hist, u_hist
    """
    if params is None:
        params = {}

    h     = params.get('h', np.pi / 4.0)
    x0    = np.asarray(params.get('x0', [np.pi, np.pi / 2.0]), dtype=float).ravel()
    u0    = params.get('u0', 1.0)
    t_end = params.get('t_end', 10.0)
    dt    = params.get('dt', 0.01)
    mode  = params.get('control_mode', 'zero')
    k1    = params.get('k1', 1.0)
    k2    = params.get('k2', 1.0)

    # Derived quantities
    Nt  = int(round(t_end / dt)) + 1
    N_h = int(round(h / dt))

    # Pre-allocate
    t      = np.arange(Nt) * dt
    x_hist = np.zeros((Nt, 2))
    u_full = np.zeros(Nt + N_h)

    # Initial conditions
    x_hist[0, :] = x0
    u_full[:N_h] = u0

    # Main time-stepping loop
    for k in range(Nt - 1):
        x1_k = x_hist[k, 0]
        x2_k = x_hist[k, 1]

        # --- Control law (no predictor) ---
        mode_lower = mode.lower()
        if mode_lower == 'zero':
            u_k = 0.0
        elif mode_lower == 'proportional':
            u_k = -k1 * x1_k - k2 * x2_k
        elif mode_lower == 'naive_pred':
            u_k = -x1_k - x2_k
        else:
            raise ValueError(f'Unknown control_mode: {mode}')

        # Store control
        u_full[N_h + k] = u_k

        # --- System dynamics ---
        u_delayed = _lookup_u(t[k] - h, dt, u_full, N_h)

        dx1 = x2_k
        dx2 = np.sin(x1_k) + u_k + u_delayed

        # Euler step
        x_hist[k + 1, 0] = x1_k + dt * dx1
        x_hist[k + 1, 1] = x2_k + dt * dx2

    # Hold last control value
    u_full[N_h + Nt - 1] = u_full[N_h + Nt - 2]
    u_hist = u_full[N_h:N_h + Nt].copy()

    print(f'  Uncompensated ({mode}): |x(t_end)| = {np.linalg.norm(x_hist[-1]):.4e}')

    return {'t': t, 'x_hist': x_hist, 'u_hist': u_hist}
