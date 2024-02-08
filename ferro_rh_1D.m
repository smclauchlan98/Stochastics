close all
clc

% Define parameters
t = -1.0;
l = 0.1;
c = 0.1;
radius = 1.0;

% Get value of s on boundary
syms ss
eqn =  ss^3 - (3.0 * sqrt(6.0) / 4.0) * ss^2 ...
    + ((3.0 * t / 4.0) - 3.0 * c^2) * ss - (3.0 * c / 4.0) == 0;
Sol = vpasolve(eqn, ss);
s_BC = 0.0;
roots = zeros(length(Sol), 1);
for i = 1:length(Sol)
    roots(i) = Sol(i);
    if isreal(roots(i)) && roots(i) > s_BC
        s_BC = roots(i);
    end
end
m_BC = sqrt(1.0 + (4.0 / 3.0) * c * s_BC);

% Set spatial step size k
N = 99;
k = 1 / (N + 1);

% Tolerance for steady state solve
tol = 10^(-6);

% Parameters important for SPDE defn and soln
kappa = 1.0;
alpha = 0.01; % Controls decay rate
sigma = 1.0; % Scale noise

% Time and time steps
T = 6;
M = 300000; % Number of time steps
dtref = T / M; % Size of time steps
dt = dtref * kappa;

% Spatial grid
r = linspace(0, radius, 1 / k)';
theta = (pi / 2.0) * linspace(0, 1, 1 / k); % quarter circle
[theta, r] = meshgrid(theta, r);
x = r .* cos(theta);
y = r .* sin(theta);

% Initial condition
% Scale to s_plus?
s0 = rand((1 / k) - 2, 1 / k);
m0 = rand((1 / k) - 2, 1 / k);

% Values of r without BC
r_ref = r(2:(1 / k) - 1, :);

% Solve without noise
it = 1;
fn1 = s_eqn_ferro(s0, m0, r_ref, t, l, c, s_BC, k);
fn2 = m_eqn_ferro(s0, m0, r_ref, l, c, m_BC, k);
while sqrt(norm(fn1)^2 + norm(fn2)^2) * k^0.5 > tol % steady solution solve
    fprintf('Iteration %d, error: %.7f\n', ...
        it, sqrt(norm(fn1)^2 + norm(fn2)^2) * k^0.5);

    k1s = s_eqn_ferro(s0, m0, r_ref, t, l, c, s_BC, k);
    k1m = m_eqn_ferro(s0, m0, r_ref, l, c, m_BC, k);
    k2s = s_eqn_ferro(s0 + dt * k1s / 2.0, m0 + dt * k1m / 2.0, ...
        r_ref, t, l, c, s_BC, k);
    k2m = m_eqn_ferro(s0 + dt * k1s / 2.0, m0 + dt * k1m / 2.0, ...
        r_ref, l, c, m_BC, k);
    k3s = s_eqn_ferro(s0 + dt * k2s / 2.0, m0 + dt * k2m / 2.0, ...
        r_ref, t, l, c, s_BC, k);
    k3m = m_eqn_ferro(s0 + dt * k2s / 2.0, m0 + dt * k2m / 2.0, ...
        r_ref, l, c, m_BC, k);
    k4s = s_eqn_ferro(s0 + dt * k3s / 2.0, m0 + dt * k3m / 2.0, ...
        r_ref, t, l, c, s_BC, k);
    k4m = m_eqn_ferro(s0 + dt * k3s / 2.0, m0 + dt * k3m / 2.0, ...
        r_ref, l, c, m_BC, k);

    sD = s0 + (dt / 6.0) * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
    mD = m0 + (dt / 6.0) * (k1m + 2.0 * k2m + 2.0 * k3m + k4m);
    s0 = sD;
    m0 = mD;

    fn1 = s_eqn_ferro(s0, m0, r_ref, t, l, c, s_BC, k);
    fn2 = m_eqn_ferro(s0, m0, r_ref, l, c, m_BC, k);

    it = it + 1;
end

n_iter = it - 1;
PDEssol = [zeros(1, 1 / k); sD; s_BC * ones(1, 1 / k)];
PDEmsol = [zeros(1, 1 / k); mD; m_BC * ones(1, 1 / k)];

figure
surf(x, y, PDEssol / s_BC, EdgeColor = 'interp')
title('s / s_f on a circle radius R=1','fontsize', 14,'fontweight', 'b')
ylabel('y', 'fontsize', 14, 'fontweight', 'b')
xlabel('x', 'fontsize', 14, 'fontweight', 'b')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
view(2)

figure
surf(x, y, PDEmsol / m_BC, EdgeColor = 'interp')
title('m / m_f on a circle radius R=1','fontsize', 14,'fontweight', 'b')
ylabel('y', 'fontsize', 14, 'fontweight', 'b')
xlabel('x', 'fontsize', 14, 'fontweight', 'b')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
view(2)

figure
PDEsaverage = sum(PDEssol, 2) / size(PDEssol, 1);
plot(r(:, 1), PDEsaverage / s_BC, '-', 'LineWidth', 2)
%title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
xlabel('r', 'fontsize', 14, 'fontweight', 'b')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
%ylim([0 1])

figure
PDEmaverage = sum(PDEmsol, 2) / size(PDEmsol, 1);
plot(r(:, 1), PDEmaverage / m_BC, '-', 'LineWidth', 2)
%title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
ylabel('m / m_f', 'fontsize', 14, 'fontweight', 'b')
xlabel('r', 'fontsize', 14, 'fontweight', 'b')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
%ylim([0 1])

%%
% Solve with noise
% Set number of simulations
% Something in literature does 10 simulations and averages
n_sim = 10;
MC = 2;
sm_error = zeros(MC, 2);
for w = 1:MC
    % Empty array to collect SPDE solutions for Monte-Carlo
    sMC = zeros(1 / k, 1 / k);
    mMC = zeros(1 / k, 1 / k);
    for i = 1:n_sim

        % Initial condition
        % Scale to s_plus?
        s0 = rand((1 / k) - 2, 1 / k);
        m0 = rand((1 / k) - 2, 1 / k);

        bj = get_twod_bj_1D(dtref, [(N - 1) / 2 - 1, (N - 1) / 2], ...
            alpha, radius);

        it = 1;
        for j = 1:M / kappa
            fprintf('Simulation %d, iteration %d of %d\n', i, it, M / kappa);

            % Use real part of noise
            dW = get_twod_dW(bj, kappa, 1);

            k1s = s_eqn_ferro(s0, m0, r_ref, t, l, c, s_BC, k) ...
                + sigma * dW / dt;
            k1m = m_eqn_ferro(s0, m0, r_ref, l, c, m_BC, k) + sigma * dW / dt;
            k2s = s_eqn_ferro(s0 + dt * k1s / 2.0, m0 + dt * k1m / 2.0, ...
                r_ref, t, l, c, s_BC, k) + sigma * dW / dt;
            k2m = m_eqn_ferro(s0 + dt * k1s / 2.0, m0 + dt * k1m / 2.0, ...
                r_ref, l, c, m_BC, k) + sigma * dW / dt;
            k3s = s_eqn_ferro(s0 + dt * k2s / 2.0, m0 + dt * k2m / 2.0, ...
                r_ref, t, l, c, s_BC, k) + sigma * dW / dt;
            k3m = m_eqn_ferro(s0 + dt * k2s / 2.0, m0 + dt * k2m / 2.0, ...
                r_ref, l, c, m_BC, k) + sigma * dW / dt;
            k4s = s_eqn_ferro(s0 + dt * k3s / 2.0, m0 + dt * k3m / 2.0, ...
                r_ref, t, l, c, s_BC, k) + sigma * dW / dt;
            k4m = m_eqn_ferro(s0 + dt * k3s / 2.0, m0 + dt * k3m / 2.0, ...
                r_ref, l, c, m_BC, k) + sigma * dW / dt;

            sS = s0 + (dt / 6.0) * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
            mS = m0 + (dt / 6.0) * (k1m + 2.0 * k2m + 2.0 * k3m + k4m);
            s0 = sS;
            m0 = mS;

            it = it + 1;
        end

        SPDEssol = [zeros(1, 1 / k); sS; s_BC * ones(1, 1 / k)];
        SPDEmsol = [zeros(1, 1 / k); mS; m_BC * ones(1, 1 / k)];
        sMC = sMC + SPDEssol;
        mMC = mMC + SPDEmsol;

        % figure
        % surf(x, y, SPDEssol / s_BC, EdgeColor = 'interp')
        % title('s / s_f on a circle radius R=1','fontsize', 14,'fontweight', 'b')
        % ylabel('y', 'fontsize', 14, 'fontweight', 'b')
        % xlabel('x', 'fontsize', 14, 'fontweight', 'b')
        % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
        % view(2)
        % %colorbar('Ticks',[0,0.25,0.5,0.75,1],...
        % %        'TickLabels',{'0','0.25','0.5','0.75','1'})
        %
        % figure
        % surf(x, y, SPDEmsol / m_BC, EdgeColor = 'interp')
        % title('m / m_f on a circle radius R=1','fontsize', 14,'fontweight', 'b')
        % ylabel('y', 'fontsize', 14, 'fontweight', 'b')
        % xlabel('x', 'fontsize', 14, 'fontweight', 'b')
        % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
        % view(2)

        % figure
        % SPDEsaverage = sum(SPDEssol, 2) / size(SPDEssol, 1);
        % plot(r(:, 1), SPDEsaverage / s_BC, '-', 'LineWidth', 2)
        % %title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
        % ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
        % xlabel('r', 'fontsize', 14, 'fontweight', 'b')
        % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
        % %ylim([0 1])
        %
        % figure
        % SPDEmaverage = sum(SPDEmsol, 2) / size(SPDEmsol, 1);
        % plot(r(:, 1), SPDEmaverage / m_BC, '-', 'LineWidth', 2)
        % %title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
        % ylabel('m / m_f', 'fontsize', 14, 'fontweight', 'b')
        % xlabel('r', 'fontsize', 14, 'fontweight', 'b')
        % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
        % %ylim([0 1])
        %
        % figure
        % hold on
        % plot(r(:, 1), PDEsaverage / s_BC, '-black','LineWidth', 2)
        % plot(r(:, 1), SPDEsaverage / s_BC, '-red', 'LineWidth', 2)
        % ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
        % xlabel('r', 'fontsize', 14, 'fontweight', 'b')
        % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
        %
        % figure
        % hold on
        % plot(r(:, 1), PDEmaverage / m_BC, '-black','LineWidth', 2)
        % plot(r(:, 1), SPDEmaverage / m_BC, '-red', 'LineWidth', 2)
        % ylabel('m / m_f', 'fontsize', 14, 'fontweight', 'b')
        % xlabel('r', 'fontsize', 14, 'fontweight', 'b')
        % set(gca, 'FontSize', 14, 'FontWeight', 'bold')

    end

    s_avg = sMC / n_sim;
    m_avg = mMC / n_sim;

    figure
    surf(x, y, s_avg / s_BC, EdgeColor = 'interp')
    title('s / s_f on a circle radius R=1','fontsize', 14,'fontweight', 'b')
    ylabel('y', 'fontsize', 14, 'fontweight', 'b')
    xlabel('x', 'fontsize', 14, 'fontweight', 'b')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    view(2)
    colorbar('Ticks',[0,0.25,0.5,0.75,1],...
        'TickLabels',{'0','0.25','0.5','0.75','1'})

    figure
    surf(x, y, m_avg / m_BC, EdgeColor = 'interp')
    title('m / m_f on a circle radius R=1','fontsize', 14,'fontweight', 'b')
    ylabel('y', 'fontsize', 14, 'fontweight', 'b')
    xlabel('x', 'fontsize', 14, 'fontweight', 'b')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    view(2)
    colorbar('Ticks',[0,0.25,0.5,0.75,1],...
        'TickLabels',{'0','0.25','0.5','0.75','1'})

    % figure
    s_avgaverage = sum(s_avg, 2) / size(s_avg, 1);
    % plot(r(:, 1), SPDEsaverage / s_BC, '-', 'LineWidth', 2)
    % %title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
    % ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
    % xlabel('r', 'fontsize', 14, 'fontweight', 'b')
    % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    % %ylim([0 1])
    %
    % figure
    m_avgaverage = sum(m_avg, 2) / size(m_avg, 1);
    % plot(r(:, 1), SPDEmaverage / m_BC, '-', 'LineWidth', 2)
    % %title('Average of 100 h realisations','fontsize', 14,'fontweight','b')
    % ylabel('m / m_f', 'fontsize', 14, 'fontweight', 'b')
    % xlabel('r', 'fontsize', 14, 'fontweight', 'b')
    % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    % %ylim([0 1])

    sm_error(w, :) = [sum(s_avg, 'all') / (98 * 100) ...
        sum(m_avg, 'all') / (98 * 100)];

    filename = sprintf(['./error/ferro1D_t%.1fl%.1fc%.1f' ...
        'M%dsigma%.2fMC%dnsim%d.csv'], t, l, c, M, sigma, MC, n_sim);
    writematrix(sm_error, filename);

    figure
    hold on
    plot(r(:, 1), PDEsaverage / s_BC, '-black','LineWidth', 2)
    plot(r(:, 1), s_avgaverage / s_BC, '-red', 'LineWidth', 2)
    ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
    xlabel('r', 'fontsize', 14, 'fontweight', 'b')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')

    figure
    hold on
    plot(r(:, 1), PDEmaverage / m_BC, '-black','LineWidth', 2)
    plot(r(:, 1), m_avgaverage / m_BC, '-red', 'LineWidth', 2)
    ylabel('m / m_f', 'fontsize', 14, 'fontweight', 'b')
    xlabel('r', 'fontsize', 14, 'fontweight', 'b')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')

    % % Compute difference between MC and deterministic
    % figure
    % surf(x, y, (s_avg - PDEssol) / s_BC, EdgeColor = 'interp')
    % title('Difference between SPDE s solution and PDE solution', ...
    %     'fontsize', 14,'fontweight', 'b')
    % ylabel('y', 'fontsize', 14, 'fontweight', 'b')
    % xlabel('x', 'fontsize', 14, 'fontweight', 'b')
    % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    % view(2)
    %
    % figure
    % surf(x, y, (m_avg - PDEmsol) / m_BC, EdgeColor = 'interp')
    % title('Difference between SPDE m solution and PDE solution', ...
    %     'fontsize', 14,'fontweight', 'b')
    % ylabel('y', 'fontsize', 14, 'fontweight', 'b')
    % xlabel('x', 'fontsize', 14, 'fontweight', 'b')
    % set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    % view(2)

end



