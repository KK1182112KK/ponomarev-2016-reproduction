"""
Explicit prediction (Ponomarev 2016, Eq. 75-88).

System (Eq. 75):
    dx1/dt = x2^2 + u(t-h)
    dx2/dt = x2 + u(t)

Simplified transformation (Eq. 86-88):
    z1 = x1 - x2 + int_{-h}^{0} u(t+theta) dtheta
    z2 = x2
    u  = -2*x2 + u_tilde

Cascade system (Eq. 82):
    dz1/dt = z2^2 - z2
    dz2/dt = -z2 + u_tilde

Lyapunov function (Eq. 83):
    V(z) = (z1 + z2*(z2-2)/2)^2 + z2^2

Feedback (Eq. 84-85):
    u_tilde = -dV/dz2 = -(2*(z1 + z2*(z2-2)/2)*(z2-1) + 2*z2)
"""

import numpy as np


def _lookup_u(t_query, dt, u_full, N_h):
    """Lookup u from full history array."""
    idx = int(round(t_query / dt)) + N_h
    idx = max(0, min(idx, len(u_full) - 1))
    return u_full[idx]


def run_example_b(params=None):
    """Run Example B: explicit prediction with cascade design.

    Parameters
    ----------
    params : dict, optional
        Override default parameters. Keys:
            h     : float, delay (default: 1.0)
            x0    : array-like, initial state [x1, x2] (default: [1.0, 1.0])
            u0    : float, initial control history value (default: 0.0)
            t_end : float, simulation end time (default: 15.0)
            dt    : float, time step (default: 0.001)

    Returns
    -------
    dict with keys: t, x_hist, u_hist, z_hist, V_hist
    """
    if params is None:
        params = {}

    h     = params.get('h', 1.0)
    x0    = np.asarray(params.get('x0', [1.0, 1.0]), dtype=float).ravel()
    u0    = params.get('u0', 0.0)
    t_end = params.get('t_end', 15.0)
    dt    = params.get('dt', 0.001)

    # Derived quantities
    Nt  = int(round(t_end / dt)) + 1
    N_h = int(round(h / dt))

    # Pre-allocate
    t      = np.arange(Nt) * dt
    x_hist = np.zeros((Nt, 2))
    z_hist = np.zeros((Nt, 2))
    V_hist = np.zeros(Nt)
    u_full = np.zeros(Nt + N_h)

    # Initial conditions
    x_hist[0, :] = x0
    u_full[:N_h] = u0

    ds = h / N_h if N_h > 0 else h

    # Main time-stepping loop
    for k in range(Nt - 1):
        x1 = x_hist[k, 0]
        x2 = x_hist[k, 1]

        # --- Compute z coordinates (Eq. 86-88) ---
        # z1 = x1 - x2 + int_{-h}^{0} u(t+theta) dtheta
        # z2 = x2

        # Distributed delay integral via trapezoidal rule
        int_u = 0.0
        for j in range(N_h + 1):
            theta_j = -h + j * ds
            t_query = t[k] + theta_j
            u_j = _lookup_u(t_query, dt, u_full, N_h)
            if j == 0 or j == N_h:
                w = 0.5
            else:
                w = 1.0
            int_u += w * u_j
        int_u *= ds

        z1 = x1 - x2 + int_u
        z2 = x2

        z_hist[k, :] = [z1, z2]

        # --- Lyapunov function (Eq. 83) ---
        sigma = z1 + z2 * (z2 - 2.0) / 2.0
        V_hist[k] = sigma**2 + z2**2

        # --- Feedback (Eq. 84-85) ---
        # dV/dz2 = 2*sigma*(z2-1) + 2*z2
        dVdz2 = 2.0 * sigma * (z2 - 1.0) + 2.0 * z2
        u_tilde = -dVdz2

        # u = -2*x2 + u_tilde (Eq. 88)
        u_k = -2.0 * x2 + u_tilde

        # Store control
        u_full[N_h + k] = u_k

        # --- Actual system dynamics (Eq. 75) ---
        # dx1/dt = x2^2 + u(t-h)
        # dx2/dt = x2 + u(t)
        u_delayed = _lookup_u(t[k] - h, dt, u_full, N_h)

        dx1 = x2**2 + u_delayed
        dx2 = x2 + u_k

        # Euler step
        x_hist[k + 1, 0] = x1 + dt * dx1
        x_hist[k + 1, 1] = x2 + dt * dx2

    # Final step: compute z and V for last time step
    x1 = x_hist[Nt - 1, 0]
    x2 = x_hist[Nt - 1, 1]
    int_u = 0.0
    for j in range(N_h + 1):
        theta_j = -h + j * ds
        t_query = t[Nt - 1] + theta_j
        u_j = _lookup_u(t_query, dt, u_full, N_h)
        if j == 0 or j == N_h:
            w = 0.5
        else:
            w = 1.0
        int_u += w * u_j
    int_u *= ds
    z1 = x1 - x2 + int_u
    z2 = x2
    z_hist[Nt - 1, :] = [z1, z2]
    sigma = z1 + z2 * (z2 - 2.0) / 2.0
    V_hist[Nt - 1] = sigma**2 + z2**2

    # Extract control history aligned with t
    u_full[N_h + Nt - 1] = u_full[N_h + Nt - 2]  # hold last
    u_hist = u_full[N_h:N_h + Nt].copy()

    print(f'  Example B: |x(t_end)| = {np.linalg.norm(x_hist[-1]):.4e}, '
          f'V(t_end) = {V_hist[-1]:.4e}')

    return {'t': t, 'x_hist': x_hist, 'u_hist': u_hist,
            'z_hist': z_hist, 'V_hist': V_hist}
