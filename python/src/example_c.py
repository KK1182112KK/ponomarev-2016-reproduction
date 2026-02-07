"""
Numerical predictor feedback for inverted pendulum
(Ponomarev 2016, Eq. 89-104). THE MAIN EXAMPLE.

System (Eq. 89):
    dx1/dt = x2
    dx2/dt = sin(x1) + u(t) + u(t-h),   h = pi/4

Predictor ODE (Eq. 46-47 specialized):
    xi'(s)  = [xi2; sin(xi1) + u(t+s-h)],   xi(0) = x(t)
    beta'(s) = A(xi(s)) * beta(s),            beta(0) = [0, 1]^T
    where A(x) = [[0, 1], [cos(x1), 0]]

Predicted state: y = xi(h)
Effective gain:  B = [0, 1]^T + beta(h)
Feedback (Eq. 104): u(t) = -B^T * y
"""

import numpy as np


def _lookup_u(t_query, dt, u_full, N_h):
    """Lookup u from full history array."""
    idx = int(round(t_query / dt)) + N_h
    idx = max(0, min(idx, len(u_full) - 1))
    return u_full[idx]


def run_example_c(params=None):
    """Run Example C: inverted pendulum with numerical predictor feedback.

    Parameters
    ----------
    params : dict, optional
        Override default parameters. Keys:
            h     : float, delay (default: pi/4)
            x0    : array-like, initial state [x1, x2] (default: [pi, pi/2])
            u0    : float, initial control history value (default: 1.0)
            t_end : float, simulation end time (default: 10.0)
            dt    : float, time step (default: 0.01)

    Returns
    -------
    dict with keys: t, x_hist, u_hist, y_hist
    """
    if params is None:
        params = {}

    h     = params.get('h', np.pi / 4.0)
    x0    = np.asarray(params.get('x0', [np.pi, np.pi / 2.0]), dtype=float).ravel()
    u0    = params.get('u0', 1.0)
    t_end = params.get('t_end', 10.0)
    dt    = params.get('dt', 0.01)

    # Derived quantities
    Nt  = int(round(t_end / dt)) + 1
    N_h = int(round(h / dt))
    ds  = h / N_h if N_h > 0 else h

    # Pre-allocate
    t      = np.arange(Nt) * dt
    x_hist = np.zeros((Nt, 2))
    y_hist = np.zeros((Nt, 2))
    u_full = np.zeros(Nt + N_h)

    # Initial conditions
    x_hist[0, :] = x0
    u_full[:N_h] = u0

    # Main time-stepping loop
    for k in range(Nt - 1):
        x1_k = x_hist[k, 0]
        x2_k = x_hist[k, 1]

        # ===== Step 1: Solve predictor ODE xi'(s) over s in [0, h] =====
        # Also solve beta ODE simultaneously
        xi   = np.array([x1_k, x2_k])
        beta = np.array([0.0, 1.0])

        for j in range(N_h):
            s_j = j * ds

            # u(t + s - h): query into history
            t_query = t[k] + s_j - h
            u_pred = _lookup_u(t_query, dt, u_full, N_h)

            # A(xi(s)) = [[0, 1], [cos(xi1), 0]]
            A_xi = np.array([[0.0, 1.0], [np.cos(xi[0]), 0.0]])

            # Euler step for xi
            dxi = np.array([xi[1], np.sin(xi[0]) + u_pred])
            xi = xi + ds * dxi

            # Euler step for beta
            dbeta = A_xi @ beta
            beta = beta + ds * dbeta

        # ===== Step 2: Predicted state =====
        y = xi.copy()
        y_hist[k, :] = y

        # ===== Step 3: Effective input gain =====
        B = np.array([0.0, 1.0]) + beta

        # ===== Step 4: Feedback (Eq. 104) =====
        u_k = -np.dot(B, y)

        # Store control
        u_full[N_h + k] = u_k

        # ===== Step 5: Actual system dynamics =====
        # dx1/dt = x2
        # dx2/dt = sin(x1) + u(t) + u(t-h)
        u_delayed = _lookup_u(t[k] - h, dt, u_full, N_h)

        dx1 = x2_k
        dx2 = np.sin(x1_k) + u_k + u_delayed

        # Euler step
        x_hist[k + 1, 0] = x1_k + dt * dx1
        x_hist[k + 1, 1] = x2_k + dt * dx2

    # Final step: compute predicted state for last time step
    xi   = x_hist[Nt - 1, :].copy()
    beta = np.array([0.0, 1.0])
    for j in range(N_h):
        s_j = j * ds
        t_query = t[Nt - 1] + s_j - h
        u_pred = _lookup_u(t_query, dt, u_full, N_h)
        A_xi = np.array([[0.0, 1.0], [np.cos(xi[0]), 0.0]])
        dxi = np.array([xi[1], np.sin(xi[0]) + u_pred])
        xi = xi + ds * dxi
        dbeta = A_xi @ beta
        beta = beta + ds * dbeta
    y_hist[Nt - 1, :] = xi

    # Hold last control value
    u_full[N_h + Nt - 1] = u_full[N_h + Nt - 2]

    # Extract control history aligned with t
    u_hist = u_full[N_h:N_h + Nt].copy()

    print(f'  Example C (x0=[{x0[0]:.2f}, {x0[1]:.2f}], h={h:.4f}): '
          f'|x(t_end)| = {np.linalg.norm(x_hist[-1]):.4e}')

    return {'t': t, 'x_hist': x_hist, 'u_hist': u_hist, 'y_hist': y_hist}
