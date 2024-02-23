function b = beta_biax(q1val, q2val, q3val)

trQ2 = q1val.^2 + q2val.^2 + q3val.^2;
trQ3 = sqrt(6.0) * q1val.^3 / 6.0 - sqrt(6.0) * q1val .* q2val.^2 / 2.0 ...
    + sqrt(6.0) * q1val .* q3val.^2 / 4.0 ...
    + 3.0 * sqrt(2.0) * q2val .* q3val.^2 / 4.0;

b =  1.0 - 6.0 * trQ3.^2 .* trQ2.^(-3);

end