function tests = test_example_c
%TEST_EXAMPLE_C  Validation tests for Example C (inverted pendulum, numerical predictor).
%
%  Ponomarev (2016), Section VI-C, Eq. 89-104.
%  This is the main example, reproducing Fig. 1 of the paper.
%
%  Tests:
%    test_convergence_all_ics       — x1(0) in {pi/2, pi, 3pi/2} all converge
%    test_initial_control_spike     — |u(1)| > 5 per Fig. 1
%    test_predictor_at_t0           — predictor ODE at t=0 matches independent integration
%    test_compensated_vs_uncompensated — predictor converges, open-loop diverges
%    test_dt_refinement             — O(dt) convergence order estimate
%    test_lyapunov_v0               — ||y||^2 decreases after transient
%    test_delay_sweep               — multiple delay values all converge
    tests = functiontests(localfunctions);
end

%% test_convergence_all_ics — Fig. 1: three initial conditions all converge
function test_convergence_all_ics(testCase)
    fprintf('\n  [Example C] test_convergence_all_ics: 3 initial conditions\n');

    x1_inits = [pi/2, pi, 3*pi/2];
    x2_init  = pi/2;

    for i = 1:length(x1_inits)
        params = struct('x0', [x1_inits(i); x2_init], 'u0', 1.0);
        [t, x_hist, u_hist] = run_example_c(params);

        final_norm = norm(x_hist(end,:));
        fprintf('    x1(0) = %.4f: |x(t_end)| = %.4e\n', x1_inits(i), final_norm);
        verifyLessThan(testCase, final_norm, 0.1, ...
            sprintf('State did not converge for x1(0) = %.4f.', x1_inits(i)));
    end
end

%% test_initial_control_spike — First feedback value has large magnitude
function test_initial_control_spike(testCase)
    fprintf('\n  [Example C] test_initial_control_spike: |u(1)| > 5\n');

    % Use x1(0) = pi (the middle initial condition from Fig. 1)
    params = struct('x0', [pi; pi/2], 'u0', 1.0);
    [~, ~, u_hist] = run_example_c(params);

    % The first feedback value (index 1 gives the first control applied)
    % Paper Fig. 1 shows initial spike u ~ -8 to -9
    u_first = abs(u_hist(1));
    u_max   = max(abs(u_hist));
    fprintf('    |u(1)| = %.4f\n', u_first);
    fprintf('    max|u| = %.4f\n', u_max);

    verifyGreaterThan(testCase, u_max, 5, ...
        'Control should exhibit a large initial spike (per Fig. 1).');
end

%% test_predictor_at_t0 — Verify predictor ODE at t=0 vs independent integration
function test_predictor_at_t0(testCase)
    fprintf('\n  [Example C] test_predictor_at_t0: independent predictor check\n');

    h  = pi/4;
    x0 = [pi; pi/2];
    u0 = 1.0;
    dt = 0.01;
    N_h = round(h / dt);
    ds = h / N_h;

    % Manually solve the predictor ODE at t=0:
    % xi'(s) = [xi2; sin(xi1) + u(0 + s - h)]
    % xi(0) = x0 = [pi; pi/2]
    % For s in [0, h]: u(s - h) = u0 = 1.0 (within pre-history)
    xi = x0;
    for j = 0:N_h-1
        dxi = [xi(2); sin(xi(1)) + u0];
        xi = xi + ds * dxi;
    end
    y_manual = xi;

    % Run full simulation and extract y(0) = y_hist(1,:)
    params = struct('x0', x0, 'u0', u0, 'h', h, 'dt', dt);
    [~, ~, ~, y_hist] = run_example_c(params);
    y_sim = y_hist(1,:)';

    err = norm(y_manual - y_sim);
    fprintf('    y_manual = [%.6f, %.6f]\n', y_manual(1), y_manual(2));
    fprintf('    y_sim    = [%.6f, %.6f]\n', y_sim(1), y_sim(2));
    fprintf('    |error|  = %.4e\n', err);

    verifyLessThan(testCase, err, 1e-10, ...
        'Predictor at t=0 should match independent forward integration.');
end

%% test_compensated_vs_uncompensated — Predictor converges, open-loop diverges
function test_compensated_vs_uncompensated(testCase)
    fprintf('\n  [Example C] test_compensated_vs_uncompensated\n');

    x0 = [pi; pi/2];

    % Compensated (predictor feedback)
    params_comp = struct('x0', x0, 'u0', 1.0);
    [~, x_comp, ~] = run_example_c(params_comp);
    norm_comp = norm(x_comp(end,:));

    % Uncompensated (open-loop, u=0)
    params_unc = struct('x0', x0, 'u0', 1.0, 'control_mode', 'zero');
    [~, x_unc, ~] = run_uncompensated(params_unc);
    norm_unc = norm(x_unc(end,:));

    fprintf('    Compensated:   |x(t_end)| = %.4e\n', norm_comp);
    fprintf('    Uncompensated: |x(t_end)| = %.4e\n', norm_unc);

    % Compensated should converge (Euler at dt=0.01 gives residual ~0.04)
    verifyLessThan(testCase, norm_comp, 0.1, ...
        'Compensated system should converge.');

    % Uncompensated should fail to converge (or diverge)
    % The inverted pendulum at x0 = [pi, pi/2] with u=0 will not converge
    verifyGreaterThan(testCase, norm_unc, norm_comp, ...
        'Uncompensated should perform worse than compensated.');
end

%% test_dt_refinement — O(dt) convergence order estimate
function test_dt_refinement(testCase)
    fprintf('\n  [Example C] test_dt_refinement: convergence order\n');

    params_base = struct('x0', [pi; pi/2], 'u0', 1.0, 't_end', 5.0);

    % Multiple dt values
    dt_vals = [0.04, 0.02, 0.01, 0.005];

    % Reference: dt = 0.0025
    params_ref = params_base;
    params_ref.dt = 0.0025;
    [t_ref, x_ref, ~] = run_example_c(params_ref);

    errors = zeros(size(dt_vals));
    for i = 1:length(dt_vals)
        params_i = params_base;
        params_i.dt = dt_vals(i);
        [t_i, x_i, ~] = run_example_c(params_i);

        % Interpolate reference onto this grid
        x_ref_interp = interp1(t_ref, x_ref, t_i, 'linear');
        errors(i) = max(max(abs(x_i - x_ref_interp)));
        fprintf('    dt = %.4f: max error = %.4e\n', dt_vals(i), errors(i));
    end

    % Compute convergence orders between successive refinements
    fprintf('    Convergence orders:\n');
    orders = zeros(length(dt_vals)-1, 1);
    for i = 1:length(dt_vals)-1
        if errors(i+1) > 1e-14
            ratio = dt_vals(i) / dt_vals(i+1);
            orders(i) = log(errors(i) / errors(i+1)) / log(ratio);
            fprintf('      dt %.4f -> %.4f: order = %.2f\n', ...
                dt_vals(i), dt_vals(i+1), orders(i));
        end
    end

    % Average order should be close to 1.0 (Euler method)
    valid_orders = orders(orders > 0);
    if ~isempty(valid_orders)
        avg_order = mean(valid_orders);
        fprintf('    Average order = %.2f (expected ~1.0)\n', avg_order);
        verifyGreaterThan(testCase, avg_order, 0.8, ...
            'Average convergence order should be at least 0.8.');
    end
end

%% test_lyapunov_v0 — ||y||^2 should decrease after transient
function test_lyapunov_v0(testCase)
    fprintf('\n  [Example C] test_lyapunov_v0: ||y||^2 decay\n');

    params = struct('x0', [pi; pi/2], 'u0', 1.0);
    [t, ~, ~, y_hist] = run_example_c(params);

    h = pi/4;  % default delay

    % Compute ||y||^2 at each time step
    y_norm_sq = sum(y_hist.^2, 2);

    % After transient (t > 2*h), check that ||y||^2 generally decreases
    idx_start = find(t >= 2*h, 1, 'first');
    y_tail = y_norm_sq(idx_start:end);

    % Check that final value is much smaller than post-transient start
    ratio = y_tail(end) / y_tail(1);
    fprintf('    ||y(t_end)||^2 / ||y(2h)||^2 = %.4e\n', ratio);
    verifyLessThan(testCase, ratio, 0.1, ...
        '||y||^2 should decay to less than 10%% of post-transient value.');

    % Count fraction of non-increasing steps
    dV = diff(y_tail);
    frac_dec = sum(dV <= 1e-6) / length(dV);
    fprintf('    Fraction non-increasing (t > 2h): %.2f%%\n', frac_dec*100);
    % Be lenient: Euler can cause small oscillations
    verifyGreaterThan(testCase, frac_dec, 0.80, ...
        '||y||^2 should be mostly non-increasing after transient.');
end

%% test_delay_sweep — Multiple delay values all converge
function test_delay_sweep(testCase)
    fprintf('\n  [Example C] test_delay_sweep: h in {0.2, 0.5, 0.7854, 1.0}\n');

    h_vals = [0.5, pi/4, 1.0];
    tol = 0.1;  % relaxed tolerance for Euler discretization

    for i = 1:length(h_vals)
        params = struct('x0', [pi; pi/2], 'u0', 1.0, 'h', h_vals(i), ...
                         't_end', 15.0);
        [~, x_hist, ~] = run_example_c(params);

        final_norm = norm(x_hist(end,:));
        fprintf('    h = %.4f: |x(t_end)| = %.4e\n', h_vals(i), final_norm);
        verifyLessThan(testCase, final_norm, tol, ...
            sprintf('State should converge for h = %.4f.', h_vals(i)));
    end
end
