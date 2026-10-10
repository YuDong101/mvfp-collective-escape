function U2 = Ueff2_fc(x, par)
U2 = -(par.u - par.sigma1) + 3 .* (par.v + par.sigma2) .* x.^2;
end

%% ============================================================
