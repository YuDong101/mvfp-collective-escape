function [m1, m2, Pss] = mvfp_steady_moments_from_F(F, par, x)
[rho, m1, m2, Pss] = mvfp_steady_density_from_F(F, par, x);
if isempty(rho)
    m1 = nan;
    m2 = nan;
    Pss = nan;
end
end

%% ============================================================
