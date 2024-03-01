function bj = get_twod_bj_1D(dtref, J, alpha, radius)

% We normalise lambda here, instead of just using j1^2 + j2^2
lambdax = 2.0 * pi * [0:J(1) + 1 -J(1):-1]' / (pi / 2.0); % Divided J by 2 in function argument
lambday = 2.0 * pi * [0:J(2) + 1 -J(2):-1]' / radius;
[lambdaxx, lambdayy] = meshgrid(lambday, lambdax);
root_qj = exp(-alpha * (lambdaxx.^2 + lambdayy.^2) / 2.0); % set decay rate noise
bj = root_qj * sqrt(dtref) * (2 * J(1) + 1) * (2 * J(2) + 1) ...
    / sqrt(radius * pi / 2.0);