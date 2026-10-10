function [residual, m1, m2, Pss] = mvfp_steady_residual(F, par, x)
[m1, m2, Pss] = mvfp_steady_moments_from_F(F, par, x);
feedback = par.sigma1 .* m1 + par.sigma2 .* m2 .* m1 + par.bias;
residual = F - feedback;
end

%% ============================================================
