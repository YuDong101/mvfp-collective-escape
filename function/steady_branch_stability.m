function stable = steady_branch_stability(F, par, x)
dF = max(1e-5, 1e-4 * max(1, abs(F)));
[~, m1p, m2p] = mvfp_steady_residual(F + dF, par, x);
[~, m1m, m2m] = mvfp_steady_residual(F - dF, par, x);
Gp = par.sigma1 .* m1p + par.sigma2 .* m2p .* m1p + par.bias;
Gm = par.sigma1 .* m1m + par.sigma2 .* m2m .* m1m + par.bias;
dGdF = (Gp - Gm) ./ (2 .* dF);
stable = double(isfinite(dGdF) && dGdF < 1);
end

%% ============================================================
