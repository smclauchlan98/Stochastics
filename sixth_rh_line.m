close all
clc

% Define parameters
t = 0.0;
eps = 1.0;
d = 1.0;
e = 0.0;
f = 1.0;
radius = 1.0;

% Get value of s on boundary
syms ss
eqn = (2.0 * (f - e) / 81.0) * ss^4 + (4.0 * e / 27.0) * ss^4 ...
    + (4.0 * d / 45.0) * ss^3 + (2.0 / 3.0) * ss^2 ...
    - (2.0 * sqrt(6.0) / 3.0) * ss + t == 0;
Sol = vpasolve(eqn, ss);
s_BC = 0.0;
roots = zeros(length(Sol), 1);
for i = 1:length(Sol)
    roots(i) = Sol(i);
    if isreal(roots(i)) && roots(i) > s_BC
        s_BC = roots(i);
    end
end

% Set spatial step size k
N = 99;
k = 1 / (N + 1);

% Error for steady state solve
error = 10^(-6);

% Parameters important for SPDE defn and soln
kappa = 1.0;
reg = [0.01 0.1 1.0]; % Controls decay rate
sigma = 1.0; % Scale noise

% Time and time steps
T = 1.0;
M = 100000; % Number of time steps (if T = 2.0)
dtref = 2.0 / 100000; % Size of time steps
dt = dtref * kappa;

% Spatial grid
r = linspace(0, radius, 1 / k)';

% Initial condition
% Scale to s_plus?
s0 = rand((1 / k) - 2, 1);

% Values of r without BC
r_ref = r(2:(1 / k) - 1);

% Solve without noise
it = 1;
fn = s_eqn_line(s0, r_ref, t, eps, d, e, f, s_BC, k);
while norm(fn) * k^0.5 > error % steady solution solve
    fprintf('Iteration %d, error: %.7f\n', it, norm(fn) * k^0.5);

    k1 = s_eqn_line(s0, r_ref, t, eps, d, e, f, s_BC, k);
    k2 = s_eqn_line(s0 + dt * k1 / 2.0, r_ref, t, eps, d, e, f, s_BC, k);
    k3 = s_eqn_line(s0 + dt * k2 / 2.0, r_ref, t, eps, d, e, f, s_BC, k);
    k4 = s_eqn_line(s0 + dt * k3 / 2.0, r_ref, t, eps, d, e, f, s_BC, k);

    sD = s0 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    s0 = sD;

    fn = s_eqn_line(s0, r_ref, t, eps, d, e, f, s_BC, k);
    it = it + 1;
end

n_iter = it - 1;
PDEsol = [0; sD; s_BC];

filename = sprintf(['./Output/SixthLine/Solutions' ...
    '/det_t%.1feps%.1f.csv'], ...
    t, eps);
TT = array2table(PDEsol);
writetable(TT, filename, 'WriteVariableNames', false);

figure
plot(r(:, 1), PDEsol / s_BC, '-', 'LineWidth', 2)
%title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
xlabel('r', 'fontsize', 14, 'fontweight', 'b')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
set(gcf,'visible','off')
ylim([0 1])
plotname = sprintf(['./Output/SixthLine/Plots' ...
    '/detsprofile_t%.1feps%.1f.png'], ...
    t, eps);
saveas(gcf, plotname, 'png')

% Solve with noise
% Set number of simulations
n_sim = 100;
MC = 1;
for p = 1:3

    s_mean = zeros(MC, 1);
    s_var = zeros(MC, 1);
    min_profile = zeros(n_sim, 1);
    % Empty array to collect SPDE solutions for Monte-Carlo
    for w = 1:MC
        sMC = zeros(1 / k, 1);
        for i = 1:n_sim

            % Initial condition
            % Scale to s_plus?
            s0 = rand((1 / k) - 2, 1);

            % Solve with noise
            bj = get_oned_bj(dtref, 98, 1.0, reg(p));

            it = 1;
            tt = 0.0;
            while tt < T
                fprintf('Simulation %d of %d, reg = %.2f, iteration %d of %d\n', ...
                    i, n_sim, reg(p), it, vpa(T / dtref));

                % Use real part of noise
                dW = get_oned_dW(bj, kappa, 0, 1);

                k1 = s_eqn_line(s0, r_ref, t, eps, d, e, f, s_BC, k) + sigma * dW / dt;
                k2 = s_eqn_line(s0 + dt * k1 / 2.0, r_ref, t, eps, d, e, f, s_BC, k) ...
                    + sigma * dW / dt;
                k3 = s_eqn_line(s0 + dt * k2 / 2.0, r_ref, t, eps, d, e, f, s_BC, k) ...
                    + sigma * dW / dt;
                k4 = s_eqn_line(s0 + dt * k3 / 2.0, r_ref, t, eps, d, e, f, s_BC, k) ...
                    + sigma * dW / dt;

                sS = s0 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                s0 = sS;

                it = it + 1;
                tt = tt + dt;
            end

            SPDEsol = [0.0; sS; s_BC];
            sMC = sMC + SPDEsol;

            min_profile(i) = min(SPDEsol);

           
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

        filename = sprintf(['./Output/SixthLine/Solutions' ...
            '/min_values_t%.1feps%.1fT%dreg%.2fnsim%d.csv'], ...
            t, eps, vpa(T), reg(p), n_sim);
        TT = array2table(min_profile);
        writetable(TT, filename, 'WriteVariableNames', false);

        s_avg = sMC / n_sim;
        filename = sprintf(['./Output/SixthLine/Solutions' ...
            '/stoch_t%.1feps%.1fT%dreg%.2fnsim%d.csv'], ...
            t, eps, vpa(T), reg(p), n_sim);
        TT = array2table(s_avg);
        writetable(TT, filename, 'WriteVariableNames', false);

        % absolute difference
        s_diff = abs(s_avg - PDEsol);
        % mean absolute difference
        s_mean(w, :) = sum(s_diff, 'all') / 98;
        % variance absolute difference
        s_var(w, :) = sum((s_diff - s_mean(w, :)).^2, 'all') / 98;

        filename = sprintf(['./error/sixthLine_t%.1feps%.1f' ...
            'T%dreg%.2fMC%dnsim%d.csv'], ...
            t, eps, vpa(T), reg(p), MC, n_sim);
        TT = array2table([s_mean s_var]);
        TT.Properties.VariableNames(1:2) = {'Mean abs difference', ...
            'Var abs difference'};
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
        plotname = sprintf(['./Output/SixthLine/Plots' ...
            '/stochsprofile_t%.1feps%.1fT%dalpha%.2fnsim%d.png'], ...
            t, eps, vpa(T), alpha(p), n_sim);
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




