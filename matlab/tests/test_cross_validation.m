function tests = test_cross_validation
%TEST_CROSS_VALIDATION  Cross-method validation tests.
%
%  Tests across different numerical methods and transformations.
%
%  Tests:
%    test_euler_vs_rk4_example_c        — Euler vs RK4 for Example C
%    test_example_b_analytical_vs_numerical — analytical z-transform vs numerical
    tests = functiontests(localfunctions);
end

%% test_euler_vs_rk4_example_c — Compare Euler to RK4 for Example C
function test_euler_vs_rk4_example_c(testCase)
    fprintf('\n  [Cross] test_euler_vs_rk4_example_c\n');

    h     = pi/4;
    x0    = [pi; pi/2];
    u0    = 1.0;
    t_end = 5.0;
    dt    = 0.01;
    Nt    = round(t_end / dt) + 1;
    N_h   = round(h / dt);

    %% --- Run the standard Euler simulation via run_example_c ---
    params_euler = struct('x0', x0, 'u0', u0, 't_end', t_end, 'dt', dt, 'h', h);
    [t_euler, x_euler, u_euler] = run_example_c(params_euler);

    %% --- Run RK4 simulation with the SAME control history ---
    % We use the control sequence from the Euler run to avoid feedback
    % coupling differences. This isolates the integration method difference.
    t_rk4    = (0:Nt-1)' * dt;
    x_rk4    = zeros(Nt, 2);
    x_rk4(1,:) = x0';

    % Build full control history (including pre-history) matching Euler run
    u_full_rk4 = zeros(Nt + N_h, 1);
    u_full_rk4(1:N_h) = u0;
    u_full_rk4(N_h + (1:Nt)) = u_euler;

    for k = 1:Nt-1
        x_k = x_rk4(k,:)';
        t_k = t_rk4(k);

        % Current and delayed control from stored Euler history
        u_k     = u_euler(k);
        u_del_k = lookup_u_rk4(t_k - h, dt, u_full_rk4, N_h);

        % RK4 stages for the plant ODE:
        % dx1/dt = x2
        % dx2/dt = sin(x1) + u(t) + u(t-h)
        % Note: u(t) and u(t-h) are piecewise constant (zero-order hold)
        % over each Euler step, so we hold them constant in the RK4 sub-steps.
        f_sys = @(x_) [x_(2); sin(x_(1)) + u_k + u_del_k];

        k1 = f_sys(x_k);
        k2 = f_sys(x_k + 0.5*dt*k1);
        k3 = f_sys(x_k + 0.5*dt*k2);
        k4 = f_sys(x_k + dt*k3);

        x_rk4(k+1,:) = (x_k + dt/6 * (k1 + 2*k2 + 2*k3 + k4))';
    end

    % Compare Euler vs RK4 (same control sequence)
    err = max(max(abs(x_euler - x_rk4)));
    fprintf('    max|x_euler - x_rk4| = %.4e (threshold: 0.1)\n', err);

    % The difference is purely due to integration accuracy (Euler vs RK4).
    % With dt=0.01, Euler error is O(dt) ~ 0.01 per step, accumulated.
    verifyLessThan(testCase, err, 0.1, ...
        'Euler and RK4 should agree within 0.1 for dt=0.01.');

    % Verify both reach similar final state
    err_final = norm(x_euler(end,:) - x_rk4(end,:));
    fprintf('    |x_euler(end) - x_rk4(end)| = %.4e\n', err_final);
end

%% test_example_b_analytical_vs_numerical — z-transform consistency
function test_example_b_analytical_vs_numerical(testCase)
    fprintf('\n  [Cross] test_example_b_analytical_vs_numerical\n');
    fprintf('    Checking z-transform (Eq. 86-88) vs y-transform (Eq. 76-77)\n');

    % Run Example B to get trajectories
    params = struct();
    params.h = 1.0;
    params.dt = 0.001;
    params.t_end = 10.0;
    params.x0 = [1; 1];
    params.u0 = 0.0;

    [t, x_hist, u_hist, z_hist, V_hist] = run_example_b(params);

    h  = params.h;
    dt = params.dt;
    Nt = length(t);
    N_h = round(h / dt);

    % === Check 1: z-transform (Eq. 86-88) consistency ===
    % z1 = x1 - x2 + int_{-h}^{0} u(t+theta) dtheta
    % z2 = x2
    %
    % Also verify the Lyapunov function V(z) = (z1 + z2*(z2-2)/2)^2 + z2^2

    % Sample every 500 steps
    check_idx = 1:500:Nt;
    max_z1_err = 0;
    max_z2_err = 0;
    max_V_err  = 0;

    for ci = 1:length(check_idx)
        idx = check_idx(ci);
        x1 = x_hist(idx, 1);
        x2 = x_hist(idx, 2);

        % Trapezoidal integral of u from t-h to t
        int_u = 0;
        ds_int = h / N_h;
        for j = 0:N_h
            theta_j = -h + j * ds_int;
            t_query = t(idx) + theta_j;
            u_idx = round(t_query / dt) + 1;
            u_idx = max(1, min(u_idx, Nt));
            if t_query < 0
                u_val = params.u0;
            else
                u_val = u_hist(u_idx);
            end
            if j == 0 || j == N_h
                w = 0.5;
            else
                w = 1.0;
            end
            int_u = int_u + w * u_val;
        end
        int_u = int_u * ds_int;

        z1_check = x1 - x2 + int_u;
        z2_check = x2;

        err_z1 = abs(z1_check - z_hist(idx, 1));
        err_z2 = abs(z2_check - z_hist(idx, 2));

        max_z1_err = max(max_z1_err, err_z1);
        max_z2_err = max(max_z2_err, err_z2);

        % Verify V(z) = (z1 + z2*(z2-2)/2)^2 + z2^2
        sigma = z_hist(idx,1) + z_hist(idx,2) * (z_hist(idx,2) - 2) / 2;
        V_check = sigma^2 + z_hist(idx,2)^2;
        err_V = abs(V_check - V_hist(idx));
        max_V_err = max(max_V_err, err_V);
    end

    fprintf('    Checked %d time points\n', length(check_idx));
    fprintf('    max|z1_recomp - z1_stored| = %.4e\n', max_z1_err);
    fprintf('    max|z2_recomp - z2_stored| = %.4e\n', max_z2_err);
    fprintf('    max|V_recomp - V_stored|   = %.4e\n', max_V_err);

    verifyLessThan(testCase, max_z2_err, 1e-12, ...
        'z2 should equal x2 exactly.');

    verifyLessThan(testCase, max_z1_err, 0.1, ...
        'z1 recomputation should agree closely with stored value.');

    verifyLessThan(testCase, max_V_err, 1e-10, ...
        'Lyapunov function should be consistent with z coordinates.');

    % === Check 2: Analytical y-transform (Eq. 76-77) ===
    % y1 = x1 + (e^(2h)-1)/2 * x2^2 + int_{-h}^{0} u(t+theta) dtheta
    % y2 = e^h * x2
    % Cascade from y: z1 = y1 - e^(-h)*y2 + (e^(-2h)-1)/2 * y2^2
    %                 z2 = e^(-h) * y2
    % Simplification should give z1 = x1 - x2 + int u, z2 = x2 (Eq. 86-88)
    %
    % We verify that the y -> z cascade transformation is self-consistent.
    fprintf('\n    Checking y-transform -> z-transform cascade consistency:\n');

    max_cascade_err = 0;
    for ci = 1:length(check_idx)
        idx = check_idx(ci);
        x1 = x_hist(idx, 1);
        x2 = x_hist(idx, 2);

        % Compute integral
        int_u = 0;
        ds_int = h / N_h;
        for j = 0:N_h
            theta_j = -h + j * ds_int;
            t_query = t(idx) + theta_j;
            u_idx = round(t_query / dt) + 1;
            u_idx = max(1, min(u_idx, Nt));
            if t_query < 0
                u_val = params.u0;
            else
                u_val = u_hist(u_idx);
            end
            if j == 0 || j == N_h
                w = 0.5;
            else
                w = 1.0;
            end
            int_u = int_u + w * u_val;
        end
        int_u = int_u * ds_int;

        % y-transform (Eq. 76-77)
        y1 = x1 + (exp(2*h) - 1)/2 * x2^2 + int_u;
        y2 = exp(h) * x2;

        % Cascade from y -> z (Eq. 79-81)
        z1_from_y = y1 - exp(-h)*y2 + (exp(-2*h) - 1)/2 * y2^2;
        z2_from_y = exp(-h) * y2;

        % Direct z-transform (Eq. 86-88)
        z1_direct = x1 - x2 + int_u;
        z2_direct = x2;

        cascade_err = max(abs(z1_from_y - z1_direct), abs(z2_from_y - z2_direct));
        max_cascade_err = max(max_cascade_err, cascade_err);
    end

    fprintf('    max|z_from_y - z_direct| = %.4e\n', max_cascade_err);
    verifyLessThan(testCase, max_cascade_err, 1e-10, ...
        'y-transform cascade should agree with direct z-transform.');
end

%% ========================================================================
%  Local helper: lookup u from full history array (for RK4 integration)
%  ========================================================================
function u_val = lookup_u_rk4(t_query, dt, u_full, N_h)
    idx = round(t_query / dt) + N_h + 1;
    idx = max(1, min(idx, length(u_full)));
    u_val = u_full(idx);
end
