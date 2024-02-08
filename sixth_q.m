close all
clc

% Define parameters
t = -1.0;
eps = 1.0;
d = 1.0;
e = 0.0;
f = 1.0;
radius = 1.0;

% Set spatial step size k
N = 99;
k = 1 / (N + 1);

% Error for steady state solve
error = 10^(-6);

% Parameters important for SPDE defn and soln
kappa = 1.0;
alpha = 0.01; % Controls decay rate
sigma = 1.0; % Scale noise

% Time and time steps
T = 0.1;
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
q1_0 = rand((1 / k) - 2, 1 / k);
q2_0 = rand((1 / k) - 2, 1 / k);
q3_0 = rand((1 / k) - 2, 1 / k);

% Values of r without BC
r_ref = r(2:(1 / k) - 1, :);

% Tolerance for steady state solve
tol = 10^(-6);

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

betaD = beta(q1D, q2D, q3D);

figure
surf(x, y, beta, EdgeColor = 'interp')
% title('s / s_f on a circle radius R=1','fontsize', 14,'fontweight', 'b')
ylabel('y', 'fontsize', 14, 'fontweight', 'b')
xlabel('x', 'fontsize', 14, 'fontweight', 'b')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
view(2)


