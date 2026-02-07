function generate_figures(opts)
%GENERATE_FIGURES  Generate all figures for Ponomarev (2016) reproduction.
%
%  generate_figures()
%  generate_figures(opts)
%
%  Options (struct fields):
%    opts.set  — 'all' (default), or specific figure: 'fig1','fig2','fig3','fig4','fig5'
%    opts.dpi  — resolution for saving (default 300)
%
%  Figures:
%    Fig 1: Paper's Fig. 1 — u(t) and x1(t) for Example C, 3 initial conditions
%    Fig 2: Example A — scalar state and control trajectories
%    Fig 3: Example B — state, cascade coordinates, and Lyapunov function
%    Fig 4: Example C — compensated vs uncompensated comparison
%    Fig 5: Example C — delay sweep h in [0.1, pi/2]

    if nargin < 1, opts = struct(); end
    fig_set = get_field(opts, 'set', 'all');
    dpi     = get_field(opts, 'dpi', 300);

    results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    switch lower(fig_set)
        case 'all'
            make_fig1(results_dir, dpi);
            make_fig2(results_dir, dpi);
            make_fig3(results_dir, dpi);
            make_fig4(results_dir, dpi);
            make_fig5(results_dir, dpi);
        case 'fig1'
            make_fig1(results_dir, dpi);
        case 'fig2'
            make_fig2(results_dir, dpi);
        case 'fig3'
            make_fig3(results_dir, dpi);
        case 'fig4'
            make_fig4(results_dir, dpi);
        case 'fig5'
            make_fig5(results_dir, dpi);
        otherwise
            error('Unknown figure set: %s', fig_set);
    end
end

%% ====================================================================
%  Fig 1: Reproduce paper's Figure 1 (Example C, 3 initial conditions)
%  ====================================================================
function make_fig1(results_dir, dpi)
    fprintf('  Generating Fig 1 (Paper Fig. 1 — Example C, 3 ICs)...\n');

    x1_inits = [pi/2, pi, 3*pi/2];
    colors = {'b', 'r', 'k'};
    labels = {'x_1(0) = \pi/2', 'x_1(0) = \pi', 'x_1(0) = 3\pi/2'};

    fig = figure('Position', [100 100 900 350], 'Visible', 'off');

    % Left: control u(t)
    subplot(1,2,1); hold on; grid on;
    for i = 1:3
        p = struct('x0', [x1_inits(i); pi/2], 'u0', 1.0);
        [t, ~, u_hist] = run_example_c(p);
        plot(t, u_hist, colors{i}, 'LineWidth', 1.2);
    end
    xlabel('t'); ylabel('u(t)');
    title('Control Input');
    legend(labels, 'Location', 'best');

    % Right: state x1(t)
    subplot(1,2,2); hold on; grid on;
    for i = 1:3
        p = struct('x0', [x1_inits(i); pi/2], 'u0', 1.0);
        [t, x_hist] = run_example_c(p);
        plot(t, x_hist(:,1), colors{i}, 'LineWidth', 1.2);
    end
    xlabel('t'); ylabel('x_1(t)');
    title('Angular Position');
    legend(labels, 'Location', 'best');

    sgtitle('Figure 1: Example C — Predictor Feedback (Ponomarev 2016)');

    save_figure(fig, fullfile(results_dir, 'fig1_paper_figure1.png'), dpi);
    fprintf('    Saved fig1_paper_figure1.png\n');
end

%% ====================================================================
%  Fig 2: Example A — scalar state and control
%  ====================================================================
function make_fig2(results_dir, dpi)
    fprintf('  Generating Fig 2 (Example A — scalar case)...\n');

    [t, x_hist, u_hist] = run_example_a();

    fig = figure('Position', [100 100 900 350], 'Visible', 'off');

    subplot(1,2,1);
    plot(t, x_hist, 'b', 'LineWidth', 1.5);
    hold on; grid on;
    yline(0, 'k--', 'LineWidth', 0.5);
    xlabel('t'); ylabel('x(t)');
    title('State');

    subplot(1,2,2);
    plot(t, u_hist, 'r', 'LineWidth', 1.5);
    hold on; grid on;
    yline(0, 'k--', 'LineWidth', 0.5);
    xlabel('t'); ylabel('u(t)');
    title('Control Input');

    sgtitle('Figure 2: Example A — Scalar Predictor Feedback');

    save_figure(fig, fullfile(results_dir, 'fig2_example_a.png'), dpi);
    fprintf('    Saved fig2_example_a.png\n');
end

%% ====================================================================
%  Fig 3: Example B — state trajectories and Lyapunov function
%  ====================================================================
function make_fig3(results_dir, dpi)
    fprintf('  Generating Fig 3 (Example B — explicit prediction)...\n');

    [t, x_hist, u_hist, z_hist, V_hist] = run_example_b();

    fig = figure('Position', [100 100 1200 400], 'Visible', 'off');

    % States x1, x2
    subplot(1,3,1);
    plot(t, x_hist(:,1), 'b', 'LineWidth', 1.2); hold on;
    plot(t, x_hist(:,2), 'r', 'LineWidth', 1.2);
    grid on;
    xlabel('t'); ylabel('x(t)');
    title('State Trajectories');
    legend('x_1', 'x_2', 'Location', 'best');

    % Cascade coordinates z1, z2
    subplot(1,3,2);
    plot(t, z_hist(:,1), 'b', 'LineWidth', 1.2); hold on;
    plot(t, z_hist(:,2), 'r', 'LineWidth', 1.2);
    grid on;
    xlabel('t'); ylabel('z(t)');
    title('Cascade Coordinates');
    legend('z_1', 'z_2', 'Location', 'best');

    % Lyapunov function
    subplot(1,3,3);
    semilogy(t, max(V_hist, 1e-16), 'k', 'LineWidth', 1.5);
    grid on;
    xlabel('t'); ylabel('V(z)');
    title('Lyapunov Function');

    sgtitle('Figure 3: Example B — Explicit Prediction');

    save_figure(fig, fullfile(results_dir, 'fig3_example_b.png'), dpi);
    fprintf('    Saved fig3_example_b.png\n');
end

%% ====================================================================
%  Fig 4: Example C — compensated vs uncompensated
%  ====================================================================
function make_fig4(results_dir, dpi)
    fprintf('  Generating Fig 4 (Example C — compensated vs uncompensated)...\n');

    % Compensated (predictor feedback)
    p_comp = struct('x0', [pi; pi/2], 'u0', 1.0);
    [t_c, x_c, u_c] = run_example_c(p_comp);

    % Uncompensated (open-loop, u=0)
    p_unc = struct('x0', [pi; pi/2], 'u0', 1.0, 'control_mode', 'zero');
    [t_u, x_u, u_u] = run_uncompensated(p_unc);

    % Uncompensated (naive proportional)
    p_naive = struct('x0', [pi; pi/2], 'u0', 1.0, 'control_mode', 'proportional');
    [t_n, x_n, u_n] = run_uncompensated(p_naive);

    fig = figure('Position', [100 100 1200 400], 'Visible', 'off');

    % x1(t)
    subplot(1,3,1);
    plot(t_c, x_c(:,1), 'b', 'LineWidth', 1.5); hold on;
    plot(t_u, x_u(:,1), 'r--', 'LineWidth', 1.2);
    plot(t_n, x_n(:,1), 'g-.', 'LineWidth', 1.2);
    grid on;
    xlabel('t'); ylabel('x_1(t)');
    title('Angular Position');
    legend('Predictor', 'Open-loop', 'Proportional', 'Location', 'best');

    % x2(t)
    subplot(1,3,2);
    plot(t_c, x_c(:,2), 'b', 'LineWidth', 1.5); hold on;
    plot(t_u, x_u(:,2), 'r--', 'LineWidth', 1.2);
    plot(t_n, x_n(:,2), 'g-.', 'LineWidth', 1.2);
    grid on;
    xlabel('t'); ylabel('x_2(t)');
    title('Angular Velocity');
    legend('Predictor', 'Open-loop', 'Proportional', 'Location', 'best');

    % u(t)
    subplot(1,3,3);
    plot(t_c, u_c, 'b', 'LineWidth', 1.5); hold on;
    plot(t_n, u_n, 'g-.', 'LineWidth', 1.2);
    grid on;
    xlabel('t'); ylabel('u(t)');
    title('Control Input');
    legend('Predictor', 'Proportional', 'Location', 'best');

    sgtitle('Figure 4: Example C — Compensated vs Uncompensated');

    save_figure(fig, fullfile(results_dir, 'fig4_compensated_comparison.png'), dpi);
    fprintf('    Saved fig4_compensated_comparison.png\n');
end

%% ====================================================================
%  Fig 5: Example C — delay sweep
%  ====================================================================
function make_fig5(results_dir, dpi)
    fprintf('  Generating Fig 5 (Example C — delay sweep)...\n');

    h_vals = linspace(0.1, pi/2, 8);
    final_norms = zeros(size(h_vals));

    fig = figure('Position', [100 100 900 350], 'Visible', 'off');

    subplot(1,2,1); hold on; grid on;
    colors_map = lines(length(h_vals));

    for i = 1:length(h_vals)
        p = struct('x0', [pi; pi/2], 'u0', 1.0, 'h', h_vals(i));
        [t, x_hist] = run_example_c(p);
        final_norms(i) = norm(x_hist(end,:));
        plot(t, x_hist(:,1), 'Color', colors_map(i,:), 'LineWidth', 1.0);
    end
    xlabel('t'); ylabel('x_1(t)');
    title('State x_1 for various h');
    leg_labels = arrayfun(@(h) sprintf('h=%.2f', h), h_vals, 'UniformOutput', false);
    legend(leg_labels, 'Location', 'best', 'FontSize', 7);

    subplot(1,2,2);
    bar(h_vals, log10(max(final_norms, 1e-16)));
    grid on;
    xlabel('Delay h'); ylabel('log_{10}(|x(t_{end})|)');
    title('Final State Norm vs Delay');

    sgtitle('Figure 5: Example C — Delay Sweep');

    save_figure(fig, fullfile(results_dir, 'fig5_delay_sweep.png'), dpi);
    fprintf('    Saved fig5_delay_sweep.png\n');
end

%% ====================================================================
%  Utility functions
%  ====================================================================

function save_figure(fig, filepath, dpi)
    print(fig, filepath, '-dpng', sprintf('-r%d', dpi));
    close(fig);
end

function val = get_field(s, name, default)
    if isfield(s, name)
        val = s.(name);
    else
        val = default;
    end
end
