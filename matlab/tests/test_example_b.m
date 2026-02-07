function tests = test_example_b
%TEST_EXAMPLE_B  Validation tests for Example B (explicit prediction, cascade).
%
%  Ponomarev (2016), Section VI-B, Eq. 75-88.
%
%  Tests:
%    test_convergence           — default params converge |x(t_end)| < 1e-2
%    test_lyapunov_decay        — V(z) non-increasing after transient
%    test_cascade_convergence   — cascade coordinates z converge
%    test_z_transform_consistency — z1 = x1 - x2 + int u holds at each step
%    test_dt_refinement         — O(dt) convergence under step-size halving
    tests = functiontests(localfunctions);
end

%% test_convergence — Default params, state converges
function test_convergence(testCase)
    fprintf('\n  [Example B] test_convergence: default params\n');

    [t, x_hist, u_hist, z_hist, V_hist] = run_example_b();

    % State norm should converge
    final_norm = norm(x_hist(end,:));
    fprintf('    |x(t_end)| = %.4e (threshold: 1e-2)\n', final_norm);
    verifyLessThan(testCase, final_norm, 1e-2, ...
        'State x did not converge to zero.');

    % Lyapunov should be small at end
    fprintf('    V(t_end) = %.4e\n', V_hist(end));
    verifyLessThan(testCase, V_hist(end), 1e-3, ...
        'Lyapunov function should be near zero at end.');

    % Dimensional checks
    verifyEqual(testCase, length(t), size(x_hist, 1));
    verifyEqual(testCase, size(x_hist, 2), 2);
    verifyEqual(testCase, size(z_hist, 2), 2);
end

%% test_lyapunov_decay — V(z) should be non-increasing after transient
function test_lyapunov_decay(testCase)
    fprintf('\n  [Example B] test_lyapunov_decay: V(z) monotonicity\n');

    params = struct();
    params.h = 1.0;
    [t, ~, ~, ~, V_hist] = run_example_b(params);

    % Allow transient period of 2*h
    transient_time = 2 * params.h;
    idx_start = find(t >= transient_time, 1, 'first');

    V_tail = V_hist(idx_start:end);
    dV = diff(V_tail);

    % Count violations (V increasing)
    n_violations = sum(dV > 1e-10);
    frac_decreasing = 1 - n_violations / length(dV);
    fprintf('    V values after t=%.1f: %d points\n', transient_time, length(V_tail));
    fprintf('    Fraction non-increasing: %.4f\n', frac_decreasing);

    % Allow up to 2% violations (Euler discretization artifacts)
    verifyGreaterThan(testCase, frac_decreasing, 0.98, ...
        'V(z) should be nearly monotonically decreasing after transient.');

    % V should decrease significantly: V(end) < 0.1 * V(transient_end)
    V_ratio = V_hist(end) / V_tail(1);
    fprintf('    V(t_end)/V(t_transient) = %.4e\n', V_ratio);
    verifyLessThan(testCase, V_ratio, 0.1, ...
        'V should decay to less than 10%% of its post-transient value.');
end

%% test_cascade_convergence — z coordinates converge to zero
function test_cascade_convergence(testCase)
    fprintf('\n  [Example B] test_cascade_convergence: |z(t_end)| < 1e-2\n');

    [~, ~, ~, z_hist, ~] = run_example_b();

    final_z_norm = norm(z_hist(end,:));
    fprintf('    |z(t_end)| = %.4e (threshold: 1e-2)\n', final_z_norm);
    verifyLessThan(testCase, final_z_norm, 1e-2, ...
        'Cascade coordinates z should converge to zero.');

    % Individual components
    fprintf('    z1(t_end) = %.4e, z2(t_end) = %.4e\n', ...
        z_hist(end,1), z_hist(end,2));
end

%% test_z_transform_consistency — Verify z1 = x1 - x2 + int u at each step
function test_z_transform_consistency(testCase)
    fprintf('\n  [Example B] test_z_transform_consistency: checking Eq. 86-88\n');

    params = struct();
    params.h = 1.0;
    params.dt = 0.001;
    [t, x_hist, u_hist, z_hist, ~] = run_example_b(params);

    h = params.h;
    dt = params.dt;
    N_h = round(h / dt);
    Nt = length(t);

    % Recompute z1 = x1 - x2 + int_{-h}^{0} u(t+theta) dtheta independently
    % z2 = x2
    max_err_z1 = 0;
    max_err_z2 = 0;
    n_check = 0;

    % Check at sampled time points (every 100 steps to keep it fast)
    check_indices = 1:100:Nt;
    for idx = check_indices
        x1 = x_hist(idx, 1);
        x2 = x_hist(idx, 2);

        % Compute integral of u from t-h to t using trapezoidal rule
        int_u = 0;
        ds = h / N_h;
        for j = 0:N_h
            theta_j = -h + j * ds;
            t_query = t(idx) + theta_j;
            % Find the u value at t_query
            u_idx = round(t_query / dt) + 1;
            u_idx = max(1, min(u_idx, Nt));
            if t_query < 0
                u_val = 0;  % pre-history is 0
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
        int_u = int_u * ds;

        z1_recomp = x1 - x2 + int_u;
        z2_recomp = x2;

        err_z1 = abs(z1_recomp - z_hist(idx, 1));
        err_z2 = abs(z2_recomp - z_hist(idx, 2));

        max_err_z1 = max(max_err_z1, err_z1);
        max_err_z2 = max(max_err_z2, err_z2);
        n_check = n_check + 1;
    end

    fprintf('    Checked %d time points\n', n_check);
    fprintf('    max|z1_recomp - z1_stored| = %.4e\n', max_err_z1);
    fprintf('    max|z2_recomp - z2_stored| = %.4e\n', max_err_z2);

    % z2 = x2 should be exact
    verifyLessThan(testCase, max_err_z2, 1e-12, ...
        'z2 should equal x2 exactly.');

    % z1 recomputation should agree closely (trapezoidal rule differences)
    verifyLessThan(testCase, max_err_z1, 0.1, ...
        'z1 recomputation should match stored z1.');
end

%% test_dt_refinement — O(dt) convergence under step-size halving
function test_dt_refinement(testCase)
    fprintf('\n  [Example B] test_dt_refinement: checking O(dt) convergence\n');

    % Use shorter t_end for speed
    params_base = struct();
    params_base.t_end = 5.0;

    % Coarse: dt = 0.002
    params_c = params_base;
    params_c.dt = 0.002;
    [t_c, x_c, ~, ~, ~] = run_example_b(params_c);

    % Fine: dt = 0.001
    params_f = params_base;
    params_f.dt = 0.001;
    [t_f, x_f, ~, ~, ~] = run_example_b(params_f);

    % Reference: dt = 0.0005
    params_r = params_base;
    params_r.dt = 0.0005;
    [t_r, x_r, ~, ~, ~] = run_example_b(params_r);

    % Interpolate onto coarse grid
    x_f_interp = interp1(t_f, x_f, t_c, 'linear');
    x_r_interp = interp1(t_r, x_r, t_c, 'linear');

    % Errors relative to reference
    err_coarse = max(max(abs(x_c - x_r_interp)));
    err_fine   = max(max(abs(x_f_interp - x_r_interp)));

    fprintf('    err(dt=0.002) = %.4e\n', err_coarse);
    fprintf('    err(dt=0.001) = %.4e\n', err_fine);

    if err_fine > 1e-14
        order = log2(err_coarse / err_fine);
        fprintf('    Estimated order = %.2f (expected ~1.0)\n', order);
        verifyGreaterThan(testCase, order, 0.7, ...
            'Convergence order should be at least 0.7.');
    else
        fprintf('    Errors too small for order estimate.\n');
    end
end
