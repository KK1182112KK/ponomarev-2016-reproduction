function [t, x_hist, u_hist, z_hist, V_hist] = run_example_b(params)
%RUN_EXAMPLE_B  Explicit prediction (Ponomarev 2016, Eq. 75-88).
%
%  [t, x_hist, u_hist, z_hist, V_hist] = run_example_b()
%  [t, x_hist, u_hist, z_hist, V_hist] = run_example_b(params)
%
%  System (Eq. 75):
%    dx1/dt = x2^2 + u(t-h)
%    dx2/dt = x2 + u(t)
%
%  Simplified transformation (Eq. 86-88):
%    z1 = x1 - x2 + int_{-h}^{0} u(t+theta) dtheta
%    z2 = x2
%    u  = -2*x2 + u_tilde
%
%  Cascade system (Eq. 82):
%    dz1/dt = z2^2 - z2
%    dz2/dt = -z2 + u_tilde
%
%  Lyapunov function (Eq. 83):
%    V(z) = (z1 + z2*(z2-2)/2)^2 + z2^2
%
%  Feedback (Eq. 84-85):
%    u_tilde = -dV/dz2 = -(2*(z1 + z2*(z2-2)/2)*(z2-1) + 2*z2)
%
%  Returns:
%    t      — (Nt x 1)  time vector
%    x_hist — (Nt x 2)  state trajectory [x1, x2]
%    u_hist — (Nt x 1)  control input u(t)
%    z_hist — (Nt x 2)  cascade coordinates [z1, z2]
%    V_hist — (Nt x 1)  Lyapunov function V(z(t))

    %% Default parameters
    if nargin < 1, params = struct(); end

    h     = get_field(params, 'h', 1.0);
    x0    = get_field(params, 'x0', [1; 1]);
    u0    = get_field(params, 'u0', 0.0);       % initial control history value
    t_end = get_field(params, 't_end', 15.0);
    dt    = get_field(params, 'dt', 0.001);

    %% Derived quantities
    Nt  = round(t_end / dt) + 1;
    N_h = round(h / dt);

    %% Pre-allocate
    t      = (0:Nt-1)' * dt;
    x_hist = zeros(Nt, 2);
    z_hist = zeros(Nt, 2);
    V_hist = zeros(Nt, 1);
    u_full = zeros(Nt + N_h, 1);    % includes pre-history

    % Initial conditions
    x_hist(1,:) = x0(:)';
    u_full(1:N_h) = u0;

    %% Main time-stepping loop
    for k = 1:Nt-1
        x1 = x_hist(k, 1);
        x2 = x_hist(k, 2);

        % --- Compute z coordinates (Eq. 86-88) ---
        % z1 = x1 - x2 + int_{-h}^{0} u(t+theta) dtheta
        % z2 = x2

        % Distributed delay integral via trapezoidal rule
        int_u = 0;
        ds = h / N_h;
        for j = 0:N_h
            theta_j = -h + j * ds;
            t_query = t(k) + theta_j;
            u_j = lookup_u(t_query, dt, u_full, N_h);
            if j == 0 || j == N_h
                w = 0.5;
            else
                w = 1.0;
            end
            int_u = int_u + w * u_j;
        end
        int_u = int_u * ds;

        z1 = x1 - x2 + int_u;
        z2 = x2;

        z_hist(k,:) = [z1, z2];

        % --- Lyapunov function (Eq. 83) ---
        sigma = z1 + z2 * (z2 - 2) / 2;
        V_hist(k) = sigma^2 + z2^2;

        % --- Feedback (Eq. 84-85) ---
        % dV/dz2 = 2*sigma*(z2-1) + 2*z2
        dVdz2 = 2 * sigma * (z2 - 1) + 2 * z2;
        u_tilde = -dVdz2;

        % u = -2*x2 + u_tilde (Eq. 88)
        u_k = -2 * x2 + u_tilde;

        % Store control
        u_full(N_h + k) = u_k;

        % --- Actual system dynamics (Eq. 75) ---
        % dx1/dt = x2^2 + u(t-h)
        % dx2/dt = x2 + u(t)
        u_delayed = lookup_u(t(k) - h, dt, u_full, N_h);

        dx1 = x2^2 + u_delayed;
        dx2 = x2 + u_k;

        % Euler step
        x_hist(k+1, 1) = x1 + dt * dx1;
        x_hist(k+1, 2) = x2 + dt * dx2;
    end

    % Final step: compute z and V for last time step
    x1 = x_hist(Nt, 1);
    x2 = x_hist(Nt, 2);
    int_u = 0;
    ds = h / N_h;
    for j = 0:N_h
        theta_j = -h + j * ds;
        t_query = t(Nt) + theta_j;
        u_j = lookup_u(t_query, dt, u_full, N_h);
        if j == 0 || j == N_h
            w = 0.5;
        else
            w = 1.0;
        end
        int_u = int_u + w * u_j;
    end
    int_u = int_u * ds;
    z1 = x1 - x2 + int_u;
    z2 = x2;
    z_hist(Nt,:) = [z1, z2];
    sigma = z1 + z2 * (z2 - 2) / 2;
    V_hist(Nt) = sigma^2 + z2^2;

    % Extract control history aligned with t
    u_full(N_h + Nt) = u_full(N_h + Nt - 1);  % hold last
    u_hist = u_full(N_h + (1:Nt));

    fprintf('  Example B: |x(t_end)| = %.4e, V(t_end) = %.4e\n', ...
            norm(x_hist(end,:)), V_hist(end));
end

%% Helper: lookup u from full history array
function u_val = lookup_u(t_query, dt, u_full, N_h)
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
