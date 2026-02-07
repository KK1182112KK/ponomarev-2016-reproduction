%% RUN_ALL  One-click reproduction of Ponomarev (2016) results
%
%  Usage:
%    run_all           — Run all 3 simulations (default 'sim' mode)
%    run_all('sim')    — Run simulations only (no figures)
%    run_all('fig')    — Generate figures only
%    run_all('test')   — Run validation tests only
%    run_all('all')    — Everything: simulation + figures + tests
%
%  Reference:
%    Ponomarev (2016), "Nonlinear Predictor Feedback for Input-Affine
%    Systems with Distributed Input Delays", IEEE TAC.
%
%  Designed for both interactive use and CI pipelines.

function run_all(mode)
    if nargin < 1, mode = 'sim'; end

    % Setup paths
    root = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(root, 'src')));

    print_banner();

    switch lower(mode)
        case 'sim'
            run_simulation();
        case 'fig'
            run_figures();
        case 'test'
            run_tests(root);
        case 'all'
            run_simulation();
            run_figures();
            run_tests(root);
        otherwise
            error('Unknown mode: %s. Use ''sim'', ''fig'', ''test'', or ''all''.', mode);
    end

    fprintf('\n=== Done ===\n');
end

%% ====================================================================
%  Banner
%  ====================================================================
function print_banner()
    fprintf('=============================================\n');
    fprintf('  Ponomarev (2016) — Reproduction\n');
    fprintf('  Nonlinear Predictor Feedback for\n');
    fprintf('  Input-Affine Systems with Distributed\n');
    fprintf('  Input Delays\n');
    fprintf('=============================================\n\n');
end

%% ====================================================================
%  Simulation
%  ====================================================================
function run_simulation()
    fprintf('--- Simulation ---\n\n');

    % Example A: Scalar case
    fprintf('[Example A] Scalar predictor feedback (Eq. 73-74)\n');
    [t_a, x_a, u_a] = run_example_a();
    fprintf('  Final: |x(%.0f)| = %.4e\n\n', t_a(end), abs(x_a(end)));

    % Example B: Explicit prediction
    fprintf('[Example B] Explicit prediction (Eq. 75-88)\n');
    [t_b, x_b, u_b, z_b, V_b] = run_example_b();
    fprintf('  Final: |x(%.0f)| = %.4e, V = %.4e\n\n', ...
            t_b(end), norm(x_b(end,:)), V_b(end));

    % Example C: Numerical prediction (inverted pendulum)
    fprintf('[Example C] Numerical prediction — Inverted pendulum (Eq. 89-104)\n');
    x1_inits = [pi/2, pi, 3*pi/2];
    for i = 1:length(x1_inits)
        p = struct('x0', [x1_inits(i); pi/2], 'u0', 1.0);
        [t_c, x_c, u_c] = run_example_c(p);
        fprintf('  x1(0) = %.4f: |x(%.0f)| = %.4e, max|u| = %.2f\n', ...
                x1_inits(i), t_c(end), norm(x_c(end,:)), max(abs(u_c)));
    end

    % Cross-validation: compensated vs uncompensated
    fprintf('\n[Comparison] Compensated vs Uncompensated (x0 = [pi, pi/2])\n');
    p_comp = struct('x0', [pi; pi/2], 'u0', 1.0);
    [~, x_comp] = run_example_c(p_comp);

    p_unc = struct('x0', [pi; pi/2], 'u0', 1.0, 'control_mode', 'zero');
    [~, x_unc] = run_uncompensated(p_unc);

    p_naive = struct('x0', [pi; pi/2], 'u0', 1.0, 'control_mode', 'proportional');
    [~, x_naive] = run_uncompensated(p_naive);

    fprintf('  Predictor:     |x(t_end)| = %.4e\n', norm(x_comp(end,:)));
    fprintf('  Open-loop:     |x(t_end)| = %.4e\n', norm(x_unc(end,:)));
    fprintf('  Proportional:  |x(t_end)| = %.4e\n', norm(x_naive(end,:)));
end

%% ====================================================================
%  Figures
%  ====================================================================
function run_figures()
    fprintf('--- Generating Figures ---\n\n');
    generate_figures(struct('set', 'all', 'dpi', 300));
end

%% ====================================================================
%  Tests
%  ====================================================================
function run_tests(root)
    fprintf('--- Running Tests ---\n\n');
    test_dir = fullfile(root, 'tests');
    if exist(test_dir, 'dir')
        results = runtests(test_dir, 'IncludeSubfolders', true);
        disp(results);
        n_failed = sum([results.Failed]);
        if n_failed > 0
            error('%d test(s) failed.', n_failed);
        end
    else
        fprintf('No tests directory found at: %s\n', test_dir);
    end
end
