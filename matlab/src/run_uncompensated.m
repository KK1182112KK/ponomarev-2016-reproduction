function [t, x_hist, u_hist] = run_uncompensated(params)
%RUN_UNCOMPENSATED  Example C without predictor feedback (no delay compensation).
%
%  [t, x_hist, u_hist] = run_uncompensated()
%  [t, x_hist, u_hist] = run_uncompensated(params)
%
%  System (same as Example C):
%    dx1/dt = x2
%    dx2/dt = sin(x1) + u(t) + u(t-h),   h = pi/4
%
%  Control: Simple proportional feedback WITHOUT predictor:
%    u(t) = -k1*x1(t) - k2*x2(t)
%  This ignores the delay and typically fails or gives poor performance.
%
%  Three modes available via params.control_mode:
%    'zero'         — u(t) = 0 (open-loop, default)
%    'proportional' — u(t) = -k1*x1 - k2*x2  (delay-ignorant)
%    'naive_pred'   — u(t) = -x1 - x2  (as if no delay)
%
%  Returns:
%    t      — (Nt x 1)  time vector
%    x_hist — (Nt x 2)  state trajectory [x1, x2]
%    u_hist — (Nt x 1)  control input

    %% Default parameters
    if nargin < 1, params = struct(); end

    h      = get_field(params, 'h', pi/4);
    x0     = get_field(params, 'x0', [pi; pi/2]);
    u0     = get_field(params, 'u0', 1.0);
    t_end  = get_field(params, 't_end', 10.0);
    dt     = get_field(params, 'dt', 0.01);
    mode   = get_field(params, 'control_mode', 'zero');
    k1     = get_field(params, 'k1', 1.0);
    k2     = get_field(params, 'k2', 1.0);

    %% Derived quantities
    Nt  = round(t_end / dt) + 1;
    N_h = round(h / dt);

    %% Pre-allocate
    t      = (0:Nt-1)' * dt;
    x_hist = zeros(Nt, 2);
    u_full = zeros(Nt + N_h, 1);

    % Initial conditions
    x_hist(1,:) = x0(:)';
    u_full(1:N_h) = u0;

    %% Main time-stepping loop
    for k = 1:Nt-1
        x1_k = x_hist(k, 1);
        x2_k = x_hist(k, 2);

        % --- Control law (no predictor) ---
        switch lower(mode)
            case 'zero'
                u_k = 0;
            case 'proportional'
                u_k = -k1 * x1_k - k2 * x2_k;
            case 'naive_pred'
                u_k = -x1_k - x2_k;
            otherwise
                error('Unknown control_mode: %s', mode);
        end

        % Store control
        u_full(N_h + k) = u_k;

        % --- System dynamics ---
        u_delayed = lookup_u(t(k) - h, dt, u_full, N_h);

        dx1 = x2_k;
        dx2 = sin(x1_k) + u_k + u_delayed;

        % Euler step
        x_hist(k+1, 1) = x1_k + dt * dx1;
        x_hist(k+1, 2) = x2_k + dt * dx2;
    end

    % Hold last control value
    u_full(N_h + Nt) = u_full(N_h + Nt - 1);
    u_hist = u_full(N_h + (1:Nt));

    fprintf('  Uncompensated (%s): |x(t_end)| = %.4e\n', mode, norm(x_hist(end,:)));
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
