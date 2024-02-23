close all
clc

% Define parameters
t = -20.0;
eps = 0.1;
d = 1.0;
e = 0.0;
f = 1.0;
radius = 1.0;

% Set spatial step size k
N = 99;
k = 1 / (N + 1);

% Parameters important for SPDE defn and soln
kappa = 1.0;
alpha = 0.01; % Controls decay rate
sigma = 1.0; % Scale noise

% Time and time steps
T = 2.0;
M = 100000; % Number of time steps
dtref = T / M; % Size of time steps
dt = dtref * kappa;

% Spatial grid
r = linspace(0, radius, 1 / k)';
theta_1 = (pi / 2.0) * linspace(0, 1, 1 / k); % quarter circle
[theta, r] = meshgrid(theta_1, r);
x = r .* cos(theta);
y = r .* sin(theta);

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

% Define vector of q BCs.
q_BC = [sqrt(2.0 / 3.0) * (1.0 - (3.0 / 2.0) * sin(theta_1).^2) * s_BC;
    (sqrt(2.0) / 2.0) * sin(theta_1).^2 * s_BC;
    sqrt(2.0) * sin(theta_1) .* cos(theta_1) * s_BC];

% Initial condition
% Scale to s_plus?
q1_0 = s_BC * rand((1 / k) - 2, 1 / k);
q2_0 = s_BC * [rand((1 / k) - 2, (1 / k) - 1) zeros((1 / k) - 2, 1)];
q3_0 = s_BC * [zeros((1 / k) - 2, 1) rand((1 / k) - 2, (1 / k) - 2) ...
    zeros((1 / k) - 2, 1)];

% Values of r without BC
r_ref = r(2:(1 / k) - 1, :);

% Tolerance for steady state solve
tol = 10^(-6);
% tol = 10;

% Solve without noise
it = 1;
fn1 = q1_eqn_sixth(q1_0, q2_0, q3_0, r_ref, t, eps, d, e, f, q_BC, k);
fn2 = q2_eqn_sixth(q1_0, q2_0, q3_0, r_ref, t, eps, d, e, f, q_BC, k);
fn3 = q3_eqn_sixth(q1_0, q2_0, q3_0, r_ref, t, eps, d, e, f, q_BC, k);
while sqrt(norm(fn1)^2 + norm(fn2)^2 + norm(fn3)^2) * k^0.5 > tol % steady solution solve
    fprintf('Iteration %d, error: %.7f\n', ...
        it, sqrt(norm(fn1)^2 + norm(fn2)^2 + norm(fn3)^2) * k^0.5); 

    k1_1 = q1_eqn_sixth(q1_0, q2_0, q3_0, r_ref, t, eps, d, e, f, q_BC, k);
    k1_2 = q2_eqn_sixth(q1_0, q2_0, q3_0, r_ref, t, eps, d, e, f, q_BC, k);
    k1_3 = q3_eqn_sixth(q1_0, q2_0, q3_0, r_ref, t, eps, d, e, f, q_BC, k);
    k2_1 = q1_eqn_sixth(q1_0 + dt * k1_1 / 2.0, q2_0 + dt * k1_2 / 2.0, ...
        q3_0 + dt * k1_3 / 2.0, r_ref, t, eps, d, e, f, q_BC, k);
    k2_2 = q2_eqn_sixth(q1_0 + dt * k1_1 / 2.0, q2_0 + dt * k1_2 / 2.0, ...
        q3_0 + dt * k1_3 / 2.0, r_ref, t, eps, d, e, f, q_BC, k);
    k2_3 = q3_eqn_sixth(q1_0 + dt * k1_1 / 2.0, q2_0 + dt * k1_2 / 2.0, ...
        q3_0 + dt * k1_3 / 2.0, r_ref, t, eps, d, e, f, q_BC, k);
    k3_1 = q1_eqn_sixth(q1_0 + dt * k2_1 / 2.0, q2_0 + dt * k2_2 / 2.0, ...
        q3_0 + dt * k2_3 / 2.0, r_ref, t, eps, d, e, f, q_BC, k);
    k3_2 = q2_eqn_sixth(q1_0 + dt * k2_1 / 2.0, q2_0 + dt * k2_2 / 2.0, ...
        q3_0 + dt * k2_3 / 2.0, r_ref, t, eps, d, e, f, q_BC, k);
    k3_3 = q3_eqn_sixth(q1_0 + dt * k2_1 / 2.0, q2_0 + dt * k2_2 / 2.0, ...
        q3_0 + dt * k2_3 / 2.0, r_ref, t, eps, d, e, f, q_BC, k);
    k4_1 = q1_eqn_sixth(q1_0 + dt * k3_1 / 2.0, q2_0 + dt * k3_2 / 2.0, ...
        q3_0 + dt * k3_3 / 2.0, r_ref, t, eps, d, e, f, q_BC, k);
    k4_2 = q2_eqn_sixth(q1_0 + dt * k3_1 / 2.0, q2_0 + dt * k3_2 / 2.0, ...
        q3_0 + dt * k3_3 / 2.0, r_ref, t, eps, d, e, f, q_BC, k);
    k4_3 = q3_eqn_sixth(q1_0 + dt * k3_1 / 2.0, q2_0 + dt * k3_2 / 2.0, ...
        q3_0 + dt * k3_3 / 2.0, r_ref, t, eps, d, e, f, q_BC, k);

    q1D = q1_0 + (dt / 6.0) * (k1_1 + 2.0 * k2_1 + 2.0 * k3_1 + k4_1);
    q2D = q2_0 + (dt / 6.0) * (k1_2 + 2.0 * k2_2 + 2.0 * k3_2 + k4_2);
    q3D = q3_0 + (dt / 6.0) * (k1_3 + 2.0 * k2_3 + 2.0 * k3_3 + k4_3);
    q1_0 = q1D;
    q2_0 = q2D;
    q3_0 = q3D;

    fn1 = q1_eqn_sixth(q1_0, q2_0, q3_0, r_ref, t, eps, d, e, f, q_BC, k);
    fn2 = q2_eqn_sixth(q1_0, q2_0, q3_0, r_ref, t, eps, d, e, f, q_BC, k);
    fn3 = q3_eqn_sixth(q1_0, q2_0, q3_0, r_ref, t, eps, d, e, f, q_BC, k);

    it = it + 1;

end

n_iter = it - 1;

q1PDEsol = [zeros(1, 1 / k); q1D; q_BC(1, :)];
q2PDEsol = [zeros(1, 1 / k); q2D; q_BC(2, :)];
q3PDEsol = [zeros(1, 1 / k); q3D; q_BC(3, :)];
betaD = beta_biax(q1PDEsol, q2PDEsol, q3PDEsol);

figure
surf(x, y, betaD, EdgeColor = 'interp')
% title('s / s_f on a circle radius R=1','fontsize', 14,'fontweight', 'b')
ylabel('y', 'fontsize', 14, 'fontweight', 'b')
xlabel('x', 'fontsize', 14, 'fontweight', 'b')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
view(2)

filename = sprintf(['./Output/SixthqSolutions' ...
    '/detq1_t%.1feps%.1f%de%df%dniter%d.csv'], ...
    t, eps, d, e, f, n_iter);
T = array2table(q1PDEsol);
writetable(T, filename, 'WriteVariableNames', false);

filename = sprintf(['./Output/SixthqSolutions' ...
    '/detq2_t%.1feps%.1f%de%df%dniter%d.csv'], ...
    t, eps, d, e, f, n_iter);
T = array2table(q2PDEsol);
writetable(T, filename, 'WriteVariableNames', false);

filename = sprintf(['./Output/SixthqSolutions' ...
    '/detq3_t%.1feps%.1f%de%df%dniter%d.csv'], ...
    t, eps, d, e, f, n_iter);
T = array2table(q3PDEsol);
writetable(T, filename, 'WriteVariableNames', false);

filename = sprintf(['./Output/SixthqSolutions' ...
    '/detbeta_t%.1feps%.1f%de%df%dniter%d.csv'], ...
    t, eps, d, e, f, n_iter);
T = array2table(betaD);
writetable(T, filename, 'WriteVariableNames', false);


