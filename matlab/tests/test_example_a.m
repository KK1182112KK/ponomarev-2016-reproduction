function tests = test_example_a
%TEST_EXAMPLE_A  Validation tests for Example A (scalar predictor feedback).
%
%  Ponomarev (2016), Section VI-A, Eq. 73-74.
%
%  Tests:
%    test_convergence          — default params converge |x(t_end)| < 1e-2
%    test_linear_predictor     — f(x)=0 linear case, faster convergence
%    test_exponential_decay    — f(x)=x linear case, qualitative convergence
%    test_dt_refinement        — O(dt) convergence under step-size halving
    tests = functiontests(localfunctions);
end

%% test_convergence — Run with default params, verify convergence
function test_convergence(testCase)
    fprintf('\n  [Example A] test_convergence: default params\n');

    [t, x_hist, u_hist] = run_example_a();

    % State should converge to zero
    final_x = abs(x_hist(end));
    fprintf('    |x(t_end)| = %.4e (threshold: 1e-2)\n', final_x);
    verifyLessThan(testCase, final_x, 1e-2, ...
        'State x did not converge to zero within tolerance.');

    % Control should also settle near zero
    final_u = abs(u_hist(end));
    fprintf('    |u(t_end)| = %.4e\n', final_u);
    verifyLessThan(testCase, final_u, 1.0, ...
        'Control u did not settle near zero.');

    % Basic dimensional checks
    verifyEqual(testCase, length(t), length(x_hist), ...
        'Dimension mismatch: t and x_hist.');
    verifyEqual(testCase, length(t), length(u_hist), ...
        'Dimension mismatch: t and u_hist.');
end

%% test_linear_predictor — f(x)=0, system is purely linear
function test_linear_predictor(testCase)
    fprintf('\n  [Example A] test_linear_predictor: f(x)=0 (linear case)\n');

    params = struct();
    params.f = @(x) 0;      % linear case: no nonlinearity
    params.b0 = 1.0;
    params.b1 = 0.5;
    params.bint = 1.0;
    params.h = 0.5;
    params.x0 = 1.0;
    params.u0 = 0.0;
    params.t_end = 10.0;
    params.dt = 0.001;

    [t, x_hist, u_hist] = run_example_a(params);

    % With f(x)=0, predictor should give exact prediction.
    % System should converge faster than nonlinear case.
    final_x = abs(x_hist(end));
    fprintf('    |x(t_end)| = %.4e (threshold: 1e-4)\n', final_x);
    verifyLessThan(testCase, final_x, 1e-4, ...
        'Linear case should converge more tightly.');

    % After transient (t > 2*h), x(t) should decay roughly exponentially
    % In the transformed coordinates dy/dt = -y, so after initial transient
    % the state should show monotone decay.
    idx_start = find(t >= 2*params.h, 1, 'first');
    x_tail = abs(x_hist(idx_start:end));
    % Check that the tail is mostly decreasing (allow small noise from Euler)
    diffs = diff(x_tail);
    frac_decreasing = sum(diffs <= 1e-10) / length(diffs);
    fprintf('    Fraction decreasing (t > 2h): %.2f%%\n', frac_decreasing*100);
    verifyGreaterThan(testCase, frac_decreasing, 0.90, ...
        'Linear case tail should be mostly monotonically decreasing.');
end

%% test_exponential_decay — f(x) = x (linear f), qualitative convergence
function test_exponential_decay(testCase)
    fprintf('\n  [Example A] test_exponential_decay: f(x) = x\n');

    params = struct();
    params.f = @(x) x;      % f(x) = x (linear, unstable open-loop)
    params.b0 = 1.0;
    params.b1 = 0.5;
    params.bint = 1.0;
    params.h = 0.5;
    params.x0 = 1.0;
    params.u0 = 0.0;
    params.t_end = 15.0;    % longer to allow convergence
    params.dt = 0.001;

    [t, x_hist, ~] = run_example_a(params);

    % Should converge despite unstable open-loop dynamics
    final_x = abs(x_hist(end));
    fprintf('    |x(t_end)| = %.4e (threshold: 1e-2)\n', final_x);
    verifyLessThan(testCase, final_x, 1e-2, ...
        'f(x)=x case should converge with predictor feedback.');

    % State should not blow up
    max_x = max(abs(x_hist));
    fprintf('    max|x| = %.4e\n', max_x);
    verifyLessThan(testCase, max_x, 100, ...
        'State should remain bounded.');
end

%% test_dt_refinement — O(dt) Euler convergence
function test_dt_refinement(testCase)
    fprintf('\n  [Example A] test_dt_refinement: checking O(dt) convergence\n');

    % Use a moderate t_end to keep tests fast
    params_base = struct();
    params_base.t_end = 5.0;

    % Coarse run: dt = 0.002
    params_coarse = params_base;
    params_coarse.dt = 0.002;
    [t_c, x_c, ~] = run_example_a(params_coarse);

    % Fine run: dt = 0.001
    params_fine = params_base;
    params_fine.dt = 0.001;
    [t_f, x_f, ~] = run_example_a(params_fine);

    % Reference run: dt = 0.0005
    params_ref = params_base;
    params_ref.dt = 0.0005;
    [t_r, x_r, ~] = run_example_a(params_ref);

    % Interpolate all onto the coarse time grid for comparison
    x_f_interp = interp1(t_f, x_f, t_c, 'linear');
    x_r_interp = interp1(t_r, x_r, t_c, 'linear');

    % Errors relative to reference
    err_coarse = max(abs(x_c - x_r_interp));
    err_fine   = max(abs(x_f_interp - x_r_interp));

    fprintf('    err(dt=0.002) = %.4e\n', err_coarse);
    fprintf('    err(dt=0.001) = %.4e\n', err_fine);

    % Convergence order: log2(err_coarse / err_fine) should be ~1 for O(dt)
    if err_fine > 1e-14
        order = log2(err_coarse / err_fine);
        fprintf('    Estimated order = %.2f (expected ~1.0)\n', order);
        verifyGreaterThan(testCase, order, 0.7, ...
            'Convergence order should be at least 0.7 (close to O(dt)).');
    else
        fprintf('    Error too small to estimate order; both are near reference.\n');
    end
end
