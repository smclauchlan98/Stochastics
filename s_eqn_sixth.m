function answer = s_eqn_sixth(sval, r, t, eps, d, e, f, s_BC, k)

% Change this to sixth-order
eqn = eps * eps * (([zeros(1, 1 / k); sval(1:(1 / k) - 3, :)] ...
    + [sval(2:(1 / k) - 2, :); s_BC * ones(1, 1 / k)] ...
    - 2.0 * sval) / (k^2) ... % Second deriv
    + 2.0 * r.^(-1) .* ([sval(2:(1 / k) - 2, :); s_BC * ones(1, 1 / k)] ... 
    - [zeros(1, 1 / k); sval(1:(1 / k) - 3, :)]) / (2.0 * k) ... % 1st deriv
    - (6.0 * sval) .* r.^(-2)) ...
    - t * sval + sqrt(6.0) * sval.^2 - (4.0 / 3.0) * sval.^3 ...
    - (2.0 * d / 9.0) * sval.^4 - (4.0 * e / 9.0) * sval.^5 ...
    - (2.0 * (f - e) / 27.0) * sval.^5;

answer = eqn;
end