function F = initial_field_guess(P, par)
A = par.u - par.sigma1;
B = par.v + par.sigma2;
x0 = sqrt(max(A ./ B, eps));
xL = -x0;
xR = x0;
[valid, x1, ~, x2] = extrema_from_F(par.bias, par);
if valid
    xL = x1;
    xR = x2;
end
M1 = P .* xL + (1 - P) .* xR;
M2 = P .* xL.^2 + (1 - P) .* xR.^2;
F = par.sigma1 .* M1 + par.sigma2 .* M2 .* M1 + par.bias;
end

%% ============================================================
