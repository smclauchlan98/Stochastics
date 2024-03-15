close all
clc

% Define parameters
t = -10.0;
l = 1.0;
% c = [0.0 0.1 1.0];
c = 0.0;
radius = 1.0;

% Set spatial step size k
N = 99;
k = 1 / (N + 1);

% Error for steady state solve
error = 10^(-6);

% Parameters important for SPDE defn and soln
kappa = 1.0;
% reg = [0.01 0.1 1.0]; % Controls decay rate
reg = 0.01;
sigma = 1.0; % Scale noise

% Time and time steps
T = 5.0;
M = 100000; % Number of time steps (if T = 2.0)
dtref = 2.0 / 100000; % Size of time steps
dt = dtref * kappa;

% Spatial grid
r = linspace(0, radius, 1 / k)';

for u = 1:3
    % Get value of s on boundary
    syms ss
    eqn =  ss^3 - (3.0 * sqrt(6.0) / 4.0) * ss^2 ...
        + ((3.0 * t / 4.0) - 3.0 * c(u)^2) * ss - (3.0 * c(u) / 4.0) == 0;
    Sol = vpasolve(eqn, ss);
    s_BC = 0.0;
    roots = zeros(length(Sol), 1);
    for i = 1:length(Sol)
        roots(i) = Sol(i);
        if isreal(roots(i)) && roots(i) > s_BC
            s_BC = roots(i);
        end
    end
    m_BC = sqrt(1.0 + (4.0 / 3.0) * c(u) * s_BC);
    % Initial condition
    % Scale to s_plus?
    s0 = rand((1 / k) - 2, 1);
    m0 = rand((1 / k) - 2, 1);

    % Values of r without BC
    r_ref = r(2:(1 / k) - 1);

    % Solve without noise
    it = 1;
    fn1 = sf_eqn_line(s0, m0, r_ref, t, l, c(u), s_BC, k);
    fn2 = mf_eqn_line(s0, m0, r_ref, l, c(u), m_BC, k);
    while sqrt(norm(fn1)^2 + norm(fn2)^2) * k^0.5 > error % steady solution solve
        fprintf('Iteration %d, error: %.7f\n', ...
            it, sqrt(norm(fn1)^2 + norm(fn2)^2) * k^0.5);

        k1s = sf_eqn_line(s0, m0, r_ref, t, l, c(u), s_BC, k);
        k1m = mf_eqn_line(s0, m0, r_ref, l, c(u), m_BC, k);
        k2s = sf_eqn_line(s0 + dt * k1s / 2.0, m0 + dt * k1m / 2.0, ...
            r_ref, t, l, c(u), s_BC, k);
        k2m = mf_eqn_line(s0 + dt * k1s / 2.0, m0 + dt * k1m / 2.0, ...
            r_ref, l, c(u), m_BC, k);
        k3s = sf_eqn_line(s0 + dt * k2s / 2.0, m0 + dt * k2m / 2.0, ...
            r_ref, t, l, c(u), s_BC, k);
        k3m = mf_eqn_line(s0 + dt * k2s / 2.0, m0 + dt * k2m / 2.0, ...
            r_ref, l, c(u), m_BC, k);
        k4s = sf_eqn_line(s0 + dt * k3s / 2.0, m0 + dt * k3m / 2.0, ...
            r_ref, t, l, c(u), s_BC, k);
        k4m = mf_eqn_line(s0 + dt * k3s / 2.0, m0 + dt * k3m / 2.0, ...
            r_ref, l, c(u), m_BC, k);

        sD = s0 + (dt / 6.0) * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
        mD = m0 + (dt / 6.0) * (k1m + 2.0 * k2m + 2.0 * k3m + k4m);
        s0 = sD;
        m0 = mD;

        fn1 = sf_eqn_line(s0, m0, r_ref, t, l, c(u), s_BC, k);
        fn2 = mf_eqn_line(s0, m0, r_ref, l, c(u), m_BC, k);

        it = it + 1;
    end

    n_iter = it - 1;
    PDEssol = [0; sD; s_BC];
    PDEmsol = [0; mD; m_BC];

    filename = sprintf(['./Output/FerroLine/Solutions' ...
        '/dets_t%.1fl%.1fc%.1f.csv'], ...
        t, l, c(u));
    TT = array2table(PDEssol);
    writetable(TT, filename, 'WriteVariableNames', false);

    filename = sprintf(['./Output/FerroLine/Solutions' ...
        '/detm_t%.1fl%.1fc%.1f.csv'], ...
        t, l, c(u));
    TT = array2table(PDEmsol);
    writetable(TT, filename, 'WriteVariableNames', false);

    figure
    plot(r(:, 1), PDEssol / s_BC, '-', 'LineWidth', 2)
    %title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
    ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
    xlabel('r', 'fontsize', 14, 'fontweight', 'b')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    set(gcf,'visible','off')
    ylim([0 1])
    plotname = sprintf(['./Output/FerroLine/Plots' ...
        '/detsprofile_t%.1fl%.1fc%.1f.png'], ...
        t, l, c(u));
    saveas(gcf, plotname, 'png')

    figure
    plot(r(:, 1), PDEmsol / m_BC, '-', 'LineWidth', 2)
    %title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
    ylabel('m / m_f', 'fontsize', 14, 'fontweight', 'b')
    xlabel('r', 'fontsize', 14, 'fontweight', 'b')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    set(gcf,'visible','off')
    ylim([0 1])
    plotname = sprintf(['./Output/FerroLine/Plots' ...
        '/detmprofile_t%.1fl%.1fc%.1f.png'], ...
        t, l, c(u));
    saveas(gcf, plotname, 'png')

    % Solve with noise
    % Set number of simulations
    n_sim = 100;
    MC = 1;
    for p = 1:3

        sm_mean = zeros(MC, 2);
        sm_var = zeros(MC, 2);
        min_profile = zeros(n_sim, 2);
        % Empty array to collect SPDE solutions for Monte-Carlo
        for w = 1:MC
            sMC = zeros(1 / k, 1);
            mMC = zeros(1 / k, 1);
            for i = 1:n_sim

                % Initial condition
                % Scale to s_plus?
                s0 = rand((1 / k) - 2, 1);
                m0 = rand((1 / k) - 2, 1);

                % Solve with noise
                bj = get_oned_bj(dtref, 98, 1.0, reg(p));

                it = 1;
                tt = 0.0;
                while tt < T
                    fprintf('Simulation %d of %d, reg = %.2f, iteration %d of %d\n', ...
                        i, n_sim, reg(p), it, vpa(T / dtref));

                    % Use real part of noise
                    dW = get_oned_dW(bj, kappa, 0, 1);

                    k1s = sf_eqn_line(s0, m0, r_ref, t, l, c(u), s_BC, k) ...
                        + sigma * dW / dt;
                    k1m = mf_eqn_line(s0, m0, r_ref, l, c(u), m_BC, k) + sigma * dW / dt;
                    k2s = sf_eqn_line(s0 + dt * k1s / 2.0, m0 + dt * k1m / 2.0, ...
                        r_ref, t, l, c(u), s_BC, k) + sigma * dW / dt;
                    k2m = mf_eqn_line(s0 + dt * k1s / 2.0, m0 + dt * k1m / 2.0, ...
                        r_ref, l, c(u), m_BC, k) + sigma * dW / dt;
                    k3s = sf_eqn_line(s0 + dt * k2s / 2.0, m0 + dt * k2m / 2.0, ...
                        r_ref, t, l, c(u), s_BC, k) + sigma * dW / dt;
                    k3m = mf_eqn_line(s0 + dt * k2s / 2.0, m0 + dt * k2m / 2.0, ...
                        r_ref, l, c(u), m_BC, k) + sigma * dW / dt;
                    k4s = sf_eqn_line(s0 + dt * k3s / 2.0, m0 + dt * k3m / 2.0, ...
                        r_ref, t, l, c(u), s_BC, k) + sigma * dW / dt;
                    k4m = mf_eqn_line(s0 + dt * k3s / 2.0, m0 + dt * k3m / 2.0, ...
                        r_ref, l, c(u), m_BC, k) + sigma * dW / dt;

                    sS = s0 + (dt / 6.0) * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
                    mS = m0 + (dt / 6.0) * (k1m + 2.0 * k2m + 2.0 * k3m + k4m);
                    s0 = sS;
                    m0 = mS;

                    it = it + 1;
                    tt = tt + dt;
                end

                SPDEssol = [0.0; sS; s_BC];
                sMC = sMC + SPDEssol;
                SPDEmsol = [0.0; mS; m_BC];
                mMC = mMC + SPDEmsol;

                min_profile(i, 1) = min(SPDEssol);
                min_profile(i, 2) = min(SPDEmsol);


                % figure
                % SPDEaverage = sum(SPDEsol, 2) / size(SPDEsol, 1);
                % plot(r(:, 1), SPDEaverage / s_BC, '-', 'LineWidth', 2)
                % %title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
                % ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
                % xlabel('r', 'fontsize', 14, 'fontweight', 'b')
                % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
                % %ylim([0 1])
                %
                % figure
                % hold on
                % plot(r(:, 1), PDEaverage / s_BC, '-black','LineWidth', 2)
                % plot(r(:, 1), SPDEaverage / s_BC, '-red', 'LineWidth', 2)
                % ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
                % xlabel('r', 'fontsize', 14, 'fontweight', 'b')
                % set(gca, 'FontSize', 14, 'FontWeight', 'bold')

            end

            filename = sprintf(['./Output/FerroLine/Solutions' ...
                '/min_values_t%.1fl%.1fc%.1fT%dreg%.2fnsim%d.csv'], ...
                t, l, c(u), vpa(T), reg(p), n_sim);
            TT = array2table(min_profile);
            writetable(TT, filename, 'WriteVariableNames', false);

            s_avg = sMC / n_sim;
            filename = sprintf(['./Output/FerroLine/Solutions' ...
                '/stochs_t%.1fl%.1fc%.1fT%dreg%.2fnsim%d.csv'], ...
                t, l, c(u), vpa(T), reg(p), n_sim);
            TT = array2table(s_avg);
            writetable(TT, filename, 'WriteVariableNames', false);

            m_avg = mMC / n_sim;
            filename = sprintf(['./Output/FerroLine/Solutions' ...
                '/stochm_t%.1fl%.1fc%.1fT%dreg%.2fnsim%d.csv'], ...
                t, l, c(u), vpa(T), reg(p), n_sim);
            TT = array2table(m_avg);
            writetable(TT, filename, 'WriteVariableNames', false);

            % absolute difference
            s_diff = abs(s_avg - PDEssol);
            m_diff = abs(m_avg - PDEmsol);
            % mean absolute difference
            sm_mean(w, :) = [sum(s_diff, 'all') / (98 * 100) ...
                sum(m_diff, 'all') / (98 * 100)];
            % variance absolute difference
            sm_var(w, :) = [sum((s_diff - sm_mean(w, 1)).^2, 'all') / (98 * 100)
                sum((m_diff - sm_mean(w, 2)).^2, 'all') / (98 * 100)];

            filename = sprintf(['./error/FerroLine_t%.1fl%.1fc%.1f' ...
                'T%dreg%.2fMC%dnsim%d.csv'], ...
                t, l, c(u), vpa(T), reg(p), MC, n_sim);
            TT = array2table([sm_mean sm_var]);
            TT.Properties.VariableNames(1:4) = {'Mean abs s difference', ...
                'Mean abs m difference', 'Var abs s difference', ...
                'Var abs m difference'};
            writetable(TT, filename);

            %         figure
            %         surf(x, y, s_avg / s_BC, EdgeColor = 'interp')
            %         title('s / s_f on a circle radius R=1', ...
            %             'fontsize', 14,'fontweight', 'b')
            %         ylabel('y', 'fontsize', 14, 'fontweight', 'b')
            %         xlabel('x', 'fontsize', 14, 'fontweight', 'b')
            %         set(gca, 'FontSize', 14, 'FontWeight', 'bold')
            %         view(2)
            %         colorbar('Ticks',[0,0.25,0.5,0.75,1],...
            %             'TickLabels',{'0','0.25','0.5','0.75','1'})

            % figure
            % s_average =  / size(s_avg, 1);
            % plot(r(:, 1), SPDEaverage / s_BC, '-', 'LineWidth', 2)
            % %title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
            % ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
            % xlabel('r', 'fontsize', 14, 'fontweight', 'b')
            % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
            % %ylim([0 1])

            figure
            hold on
            plot(r(:, 1), s_avg / s_BC, '-red', 'LineWidth', 2)
            ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
            xlabel('r', 'fontsize', 14, 'fontweight', 'b')
            set(gca, 'FontSize', 14, 'FontWeight', 'bold')
            set(gcf,'visible','off')
            plotname = sprintf(['./Output/FerroLine/Plots' ...
                '/stochsprofile_t%.1fl%.1fc%.1fT%dreg%.2fnsim%d.png'], ...
                t, l, c(u), vpa(T), reg(p), n_sim);
            saveas(gcf, plotname, 'png')

            figure
            hold on
            plot(r(:, 1), m_avg / m_BC, '-red', 'LineWidth', 2)
            ylabel('m / m_f', 'fontsize', 14, 'fontweight', 'b')
            xlabel('r', 'fontsize', 14, 'fontweight', 'b')
            set(gca, 'FontSize', 14, 'FontWeight', 'bold')
            set(gcf,'visible','off')
            plotname = sprintf(['./Output/FerroLine/Plots' ...
                '/stochmprofile_t%.1fl%.1fc%.1fT%dreg%.2fnsim%d.png'], ...
                t, l, c(u), vpa(T), reg(p), n_sim);
            saveas(gcf, plotname, 'png')


            %         figure
            %         hold on
            %         plot(r(:, 1), PDEaverage / s_BC, '-black','LineWidth', 2)
            %         plot(r(:, 1), s_avgaverage / s_BC, '-red', 'LineWidth', 2)
            %         ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
            %         xlabel('r', 'fontsize', 14, 'fontweight', 'b')
            %         set(gca, 'FontSize', 14, 'FontWeight', 'bold')

            %         % Compute difference between MC and deterministic
            %         figure
            %         surf(x, y, (s_avg - PDEsol) / s_BC, EdgeColor = 'interp')
            %         title('Difference between SPDE solution and PDE solution', ...
            %             'fontsize', 14,'fontweight', 'b')
            %         ylabel('y', 'fontsize', 14, 'fontweight', 'b')
            %         xlabel('x', 'fontsize', 14, 'fontweight', 'b')
            %         set(gca, 'FontSize', 14, 'FontWeight', 'bold')
            %         view(2)

        end
    end
end