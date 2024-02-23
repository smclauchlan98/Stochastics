function eqn2 = q2_eqn_sixth(q1val, q2val, q3val, r, t, eps, d, e, f, q_BC, k)

eqn2 = eps * eps * (([zeros(1, 1 / k); q2val(1:(1 / k) - 3, :)] ...
    + [q2val(2:(1 / k) - 2, :); q_BC(2, :)] ...
    - 2.0 * q2val) / (k^2) ... % Second r deriv
    + r.^(-1) .* ([q2val(:, 2:(1 / k)) zeros((1 / k) - 2, 1)] ...
    + [zeros((1 / k) - 2, 1) q2val(:, 1:(1 / k) - 1)] ...
    - 2.0 * q2val) / (k^2) ... % Second theta deriv
    + r.^(-2) .* ([q2val(2:(1 / k) - 2, :); q_BC(2, :)] ...
    - [zeros(1, 1 / k); q2val(1:(1 / k) - 3, :)]) / (2.0 * k)) ... % First r deriv
    - (4.0 * eps * eps * q2val .* r.^(-2) + t * q2val ...
    + 6.0 * q1val .* q2val - 3.0 * sqrt(3.0) * q3val.^2 / 2.0 ...
    + 2.0 * q2val.^3 + 2.0 * q1val.^2 .* q2val ...
    + 2.0 * q2val .* q3val.^2 ...
    + (d / 5.0) * (- 2.0 * sqrt(6.0) * q1val.^3 .* q2val / 3.0 ...
    - 2.0 * sqrt(6.0) * q1val .* q2val.^3 ...
    - sqrt(6.0) * q1val .* q2val .* q3val.^2 / 2.0 ...
    + 3.0 * sqrt(2.0) * q1val.^2 .* q3val.^2 / 4.0 ...
    + 9.0 * sqrt(2.0) * q2val.^2 .* q3val.^2 / 4.0 ...
    + 3.0 * sqrt(2.0) * q3val.^4 / 4.0) ...
    + e * (q2val.^5 + q1val.^4 .* q2val + 2.0 * q1val.^2 .* q2val.^3 ...
    + 2.0 * q2val.^3 .* q3val.^2 + q2val .* q3val.^4 ...
    + 2.0 * q1val.^2 .* q2val .* q3val.^2) ...
    + ((f - e) / 6.0) * (- 2.0 * q1val.^4 .* q2val ...
    + sqrt(3.0) * q1val.^3 .* q3val.^2 / 2.0 ...
    + 6.0 * q1val.^2 .* q2val.^3 - 3.0 * q1val.^2 .* q2val .* q3val.^2 ...
    - 9.0 * sqrt(3.0) * q1val .* q2val.^2 .* q3val.^2 / 2.0 ...
    + 3.0 * sqrt(3.0) * q1val .* q3val.^4 / 4.0 ...
    + 9.0 * q2val .* q3val.^4 / 4.0));

eqn2(:, end) = zeros((1 / k) - 2, 1);

end