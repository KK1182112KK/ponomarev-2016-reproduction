function [t, x_hist, u_hist, y_hist] = run_example_c(params)
%RUN_EXAMPLE_C  Numerical predictor feedback for inverted pendulum
%  (Ponomarev 2016, Eq. 89-104). THE MAIN EXAMPLE.
%
%  [t, x_hist, u_hist, y_hist] = run_example_c()
%  [t, x_hist, u_hist, y_hist] = run_example_c(params)
%
%  System (Eq. 89):
%    dx1/dt = x2
%    dx2/dt = sin(x1) + u(t) + u(t-h),   h = pi/4
%
%  Predictor ODE (Eq. 46-47 specialized):
%    xi'(s)  = [xi2; sin(xi1) + u(t+s-h)],   xi(0) = x(t)
%    beta'(s) = A(xi(s)) * beta(s),            beta(0) = [0;1]
%    where A(x) = [0, 1; cos(x1), 0]
%
%  Predicted state: y = xi(h)
%  Effective gain:  B = [0;1] + beta(h)
%  Feedback (Eq. 104): u(t) = -B' * y
%
%  Returns:
%    t      — (Nt x 1)  time vector
%    x_hist — (Nt x 2)  state trajectory [x1, x2]
%    u_hist — (Nt x 1)  control input u(t)
%    y_hist — (Nt x 2)  predicted state y(t) = [y1, y2]

    %% Default parameters
    if nargin < 1, params = struct(); end

    h     = get_field(params, 'h', pi/4);
    x0    = get_field(params, 'x0', [pi; pi/2]);
    u0    = get_field(params, 'u0', 1.0);           % initial control u(theta)=1 for theta in [-h,0]
    t_end = get_field(params, 't_end', 10.0);
    dt    = get_field(params, 'dt', 0.01);

    %% Derived quantities
    Nt     = round(t_end / dt) + 1;
    N_h    = round(h / dt);            % sub-steps in [0, h]
    ds     = h / N_h;                  % predictor sub-step size

    %% Pre-allocate
    t      = (0:Nt-1)' * dt;
    x_hist = zeros(Nt, 2);
    y_hist = zeros(Nt, 2);
    u_full = zeros(Nt + N_h, 1);      % includes pre-history

    % Initial conditions
    x_hist(1,:) = x0(:)';
    u_full(1:N_h) = u0;               % u(theta) = 1 for theta in [-h, 0)

    %% Main time-stepping loop
    for k = 1:Nt-1
        x1_k = x_hist(k, 1);
        x2_k = x_hist(k, 2);

        % ===== Step 1: Solve predictor ODE xi'(s) over s in [0, h] =====
        % Also solve beta ODE simultaneously (they share the xi trajectory)
        xi = [x1_k; x2_k];            % xi(0) = x(t)
        beta = [0; 1];                % beta(0) = B0 = [0;1]

        % Store xi trajectory for beta ODE (both integrated simultaneously)
        for j = 0:N_h-1
            s_j = j * ds;

            % u(t + s - h): query into history
            t_query = t(k) + s_j - h;
            u_pred = lookup_u(t_query, dt, u_full, N_h);

            % A(xi(s)) = [0, 1; cos(xi1), 0]
            A_xi = [0, 1; cos(xi(1)), 0];

            % Euler step for xi
            dxi = [xi(2); sin(xi(1)) + u_pred];
            xi = xi + ds * dxi;

            % Euler step for beta
            dbeta = A_xi * beta;
            beta = beta + ds * dbeta;
        end

        % ===== Step 2: Predicted state =====
        y = xi;                        % y(t) = xi(h)
        y_hist(k,:) = y';

        % ===== Step 3: Effective input gain =====
        B = [0; 1] + beta;            % B(y, phi) = B1 + beta(h)

        % ===== Step 4: Feedback (Eq. 104) =====
        u_k = -B' * y;

        % Store control
        u_full(N_h + k) = u_k;

        % ===== Step 5: Actual system dynamics =====
        % dx1/dt = x2
        % dx2/dt = sin(x1) + u(t) + u(t-h)
        u_delayed = lookup_u(t(k) - h, dt, u_full, N_h);

        dx1 = x2_k;
        dx2 = sin(x1_k) + u_k + u_delayed;

        % Euler step
        x_hist(k+1, 1) = x1_k + dt * dx1;
        x_hist(k+1, 2) = x2_k + dt * dx2;
    end

    % Final step: compute predicted state for last time step
    % (run predictor one more time)
    xi = x_hist(Nt,:)';
    beta = [0; 1];
    for j = 0:N_h-1
        s_j = j * ds;
        t_query = t(Nt) + s_j - h;
        u_pred = lookup_u(t_query, dt, u_full, N_h);
        A_xi = [0, 1; cos(xi(1)), 0];
        dxi = [xi(2); sin(xi(1)) + u_pred];
        xi = xi + ds * dxi;
        dbeta = A_xi * beta;
        beta = beta + ds * dbeta;
    end
    y_hist(Nt,:) = xi';

    % Hold last control value
    u_full(N_h + Nt) = u_full(N_h + Nt - 1);

    % Extract control history aligned with t
    u_hist = u_full(N_h + (1:Nt));

    fprintf('  Example C (x0=[%.2f, %.2f], h=%.4f): |x(t_end)| = %.4e\n', ...
            x0(1), x0(2), h, norm(x_hist(end,:)));
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
