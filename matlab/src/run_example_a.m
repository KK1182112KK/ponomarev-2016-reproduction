function [t, x_hist, u_hist] = run_example_a(params)
%RUN_EXAMPLE_A  Scalar predictor feedback (Ponomarev 2016, Eq. 73-74).
%
%  [t, x_hist, u_hist] = run_example_a()
%  [t, x_hist, u_hist] = run_example_a(params)
%
%  System:
%    dx/dt = f(x) + b0*u(t) + b1*u(t-h) + int_{-h}^{0} bint*u(t+theta) dtheta
%
%  Predictor transformed system:
%    dy/dt = f(y) + B(y,phi)*u(t)
%
%  Feedback (Eq. 74):
%    kappa = -(f(y) + y) / B(y, phi)
%
%  Returns:
%    t       — (Nt x 1)  time vector
%    x_hist  — (Nt x 1)  state trajectory
%    u_hist  — (Nt x 1)  control input trajectory

    %% Default parameters
    if nargin < 1, params = struct(); end

    f_func = get_field(params, 'f', @(x) sin(x));
    b0     = get_field(params, 'b0', 1.0);
    b1     = get_field(params, 'b1', 0.5);
    bint   = get_field(params, 'bint', 1.0);   % constant kernel
    h      = get_field(params, 'h', 0.5);
    x0     = get_field(params, 'x0', 1.0);
    u0     = get_field(params, 'u0', 0.0);      % initial control history value
    t_end  = get_field(params, 't_end', 10.0);
    dt     = get_field(params, 'dt', 0.001);

    %% Derived quantities
    Nt = round(t_end / dt) + 1;
    N_h = round(h / dt);           % number of steps in one delay window

    %% Pre-allocate
    t      = (0:Nt-1)' * dt;
    x_hist = zeros(Nt, 1);
    u_full = zeros(Nt + N_h, 1);   % includes pre-history indices

    % Initial conditions
    x_hist(1) = x0;
    u_full(1:N_h) = u0;            % u(theta) for theta in [-h, 0)

    %% Jacobian of f (scalar)
    % For scalar case: A(x) = df/dx. Use finite differences.
    df = @(x) (f_func(x + 1e-7) - f_func(x - 1e-7)) / 2e-7;

    %% Main time-stepping loop
    for k = 1:Nt-1
        x_k = x_hist(k);

        % --- Solve predictor ODE: xi'(s) = f(xi) + b1*u(t+s-h) + int ---
        % and beta ODE: beta'(s) = A(xi(s))*beta(s) + bint
        % simultaneously over s in [0, h]
        xi   = x_k;
        beta = b0;              % beta(0) = B0 = b0 for scalar case
        ds   = h / N_h;

        % Pre-compute distributed delay integral at s=0:
        % I(0) = int_{-h}^{0} bint * u(t + alpha) d_alpha
        % Uses trapezoidal rule over stored history.
        int_pred = 0;
        for jj = 0:N_h
            t_q = t(k) - h + jj * ds;
            u_q = lookup_u(t_q, dt, u_full, N_h, h);
            if jj == 0 || jj == N_h
                ww = 0.5;
            else
                ww = 1.0;
            end
            int_pred = int_pred + ww * bint * u_q;
        end
        int_pred = int_pred * ds;

        for j = 0:N_h-1
            s_j = j * ds;

            % u(t + s - h): time = t(k) + s_j - h
            t_query_pred = t(k) + s_j - h;
            u_pred = lookup_u(t_query_pred, dt, u_full, N_h, h);

            % A(xi) for scalar = df/dx evaluated at xi
            A_xi = df(xi);

            % Euler step for predictor
            % xi'(s) = f(xi) + b1*u(t+s-h) + int_{-h}^{-s} bint*u(t+s+theta) dtheta
            xi   = xi + ds * (f_func(xi) + b1 * u_pred + int_pred);

            % Euler step for beta
            beta = beta + ds * (A_xi * beta + bint);

            % Incremental update: I(j+1) = I(j) - ds * bint * u(t+s_j-h)
            int_pred = int_pred - ds * bint * u_pred;
        end

        y = xi;                 % predicted state
        B = b1 + beta;          % effective input gain

        % --- Feedback (Eq. 74) ---
        if abs(B) < 1e-10
            u_k = 0;            % safety: avoid division by zero
        else
            u_k = -(f_func(y) + y) / B;
        end

        % Store control
        u_full(N_h + k) = u_k;

        % --- Actual system dynamics ---
        % dx/dt = f(x) + b0*u(t) + b1*u(t-h) + int_{-h}^{0} bint*u(t+theta) dtheta

        % u(t-h)
        u_delayed = lookup_u(t(k) - h, dt, u_full, N_h, h);

        % Distributed delay integral via trapezoidal rule
        int_val = 0;
        for j = 0:N_h
            theta_j = -h + j * ds;
            t_query_int = t(k) + theta_j;
            u_j = lookup_u(t_query_int, dt, u_full, N_h, h);
            if j == 0 || j == N_h
                w = 0.5;
            else
                w = 1.0;
            end
            int_val = int_val + w * bint * u_j;
        end
        int_val = int_val * ds;

        dxdt = f_func(x_k) + b0 * u_k + b1 * u_delayed + int_val;

        % Euler step
        x_hist(k+1) = x_k + dt * dxdt;
    end

    % Final control (hold last value)
    u_full(N_h + Nt) = u_full(N_h + Nt - 1);

    % Extract control history aligned with t
    u_hist = u_full(N_h + (1:Nt));

    fprintf('  Example A: |x(t_end)| = %.4e\n', abs(x_hist(end)));
end

%% Helper: lookup u from full history array (including pre-history)
function u_val = lookup_u(t_query, dt, u_full, N_h, h)
    % u_full is indexed so that u_full(N_h + 1) = u(0),
    % u_full(1) = u(-h), etc.
    % Index for time t_query: idx = round(t_query / dt) + N_h + 1
    idx = round(t_query / dt) + N_h + 1;
    idx = max(1, min(idx, length(u_full)));
    u_val = u_full(idx);
end

%% Helper: get field with default
function val = get_field(s, name, default)
    if isfield(s, name)
        val = s.(name);
    else
        val = default;
    end
end
