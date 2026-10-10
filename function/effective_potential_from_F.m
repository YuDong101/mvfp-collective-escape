function U = effective_potential_from_F(x, par, F)
A = par.u - par.sigma1;
B = par.v + par.sigma2;
U = -0.5 .* A .* x.^2 + 0.25 .* B .* x.^4 - F .* x;
end

%% ============================================================
