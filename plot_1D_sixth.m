clc
close all

% Define parameters
t = 0.0;
eps = 1.0;
d = 1.0;
e = 0.0;
f = 1.0;

r = linspace(0, 1, 100);

alpha = 0.01;
M = 100000;
n_sim = 100;

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

filename = sprintf(['./Output/Solutions' ...
    '/det_t%.1feps%.1f%de%df%d.csv'], ...
    t, eps, d, e, f);
PDEsol = table2array(readtable(filename));

% filename = sprintf(['./Output/Solutions' ...
%     '/stoch_t%.1feps%.1fd%de%df%dM%ddalpha%.2fnsim%d.csv'], ...
%     t, eps, d, e, f, M, alpha, n_sim');
% SPDEsol = table2array(readtable(filename));

PDEaverage = sum(PDEsol, 2) / size(PDEsol, 1);
% SPDEaverage = sum(SPDEsol, 2) / size(SPDEsol, 1);

% PDE profile
figure
plot(r, PDEaverage / s_BC, '-', 'LineWidth', 2)
ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
xlabel('r', 'fontsize', 14, 'fontweight', 'b')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

% SPDE profile
figure
plot(r(:, 1), SPDEaverage / s_BC, '-red', 'LineWidth', 2)
ylabel('s / s_f', 'fontsize', 14, 'fontweight', 'b')
xlabel('r', 'fontsize', 14, 'fontweight', 'b')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')


