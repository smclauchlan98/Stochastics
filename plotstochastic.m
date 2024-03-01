clc
close all

t = 0.0;
eps = 1.0;
r = linspace(0, 1, 100);

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

filename = sprintf(['./Output/Sixth/Solutions' ...
    '/det_t%.1feps%.1f.csv'], ...
    t, eps);
detprofile = table2array(readtable(filename));
det_avg = sum(detprofile, 2) /  size(detprofile, 1);

T = 1;
n_sim = 100;
alpha = 0.01;

filename = sprintf(['./Output/Sixth/Solutions' ...
            '/stoch_t%.1feps%.1fT%dalpha%.2fnsim%d.csv'], ...
            t, eps, T, alpha, n_sim);
stochprofile1 = table2array(readtable(filename));
stoch1_avg = sum(stochprofile1, 2) /  size(stochprofile1, 1);

alpha = 0.1;

filename = sprintf(['./Output/Sixth/Solutions' ...
            '/stoch_t%.1feps%.1fT%dalpha%.2fnsim%d.csv'], ...
            t, eps, T, alpha, n_sim);
stochprofile2 = table2array(readtable(filename));
stoch2_avg = sum(stochprofile2, 2) /  size(stochprofile2, 1);

alpha = 1.0;

filename = sprintf(['./Output/Sixth/Solutions' ...
            '/stoch_t%.1feps%.1fT%dalpha%.2fnsim%d.csv'], ...
            t, eps, T, alpha, n_sim);
stochprofile3 = table2array(readtable(filename));
stoch3_avg = sum(stochprofile3, 2) /  size(stochprofile3, 1);


figure
hold on
plot(r(:), det_avg / s_BC, 'black');
plot(r(:), stoch1_avg / s_BC, 'red');
plot(r(:), stoch2_avg / s_BC, 'blue');
plot(r(:), stoch3_avg / s_BC, 'green');
ylim([-0.05 1]);




