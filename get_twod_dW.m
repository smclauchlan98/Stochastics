function [dW1, dW2] = get_twod_dW(bj, kappa, MM) % bj coefficeints kappa = dt 
J = size(bj);
if (kappa == 1) % get xi_j
    nnr = randn(J(1), J(2), MM); 
    nnc = randn(J(1), J(2), MM);
else 
    nnr = squeeze(sum(randn(J(1), J(2), MM, kappa), 4));
    nnc = squeeze(sum(randn(J(1), J(2), MM, kappa), 4));
end
nn2 = nnr + sqrt(-1) * nnc; 
tmphat = bsxfun(@times, bj, nn2);
tmp = ifft2(tmphat);
dW1 = real(tmp);
dW2 = imag(tmp);