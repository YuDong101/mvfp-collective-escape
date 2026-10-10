function [Pmean, Pstd, vals] = langevin_final_left_fraction(par, P0, T, cfg, seed)
rng(seed, 'twister');
N = cfg.N;
vals = nan(1, cfg.nRealizations);
info0 = frozenP_initial_info(P0, par);
[xL, xR] = initial_representatives(info0, par);
for r = 1:cfg.nRealizations
    nL = round(P0 * N);
    x = [xL + 0.02 .* randn(nL, 1); xR + 0.02 .* randn(N - nL, 1)];
    x = x(randperm(N));
    nSteps = max(1, ceil(T / cfg.dt));
    dt = T / nSteps;
    for step = 1:nSteps
        M1 = mean(x);
        M2 = mean(x.^2);
        drift = (par.u - par.sigma1) .* x - (par.v + par.sigma2) .* x.^3 ...
            + par.sigma1 .* M1 + par.sigma2 .* M2 .* M1 + par.bias;
        x = x + drift .* dt + sqrt(2 .* par.D .* dt) .* randn(size(x));
        x = min(max(x, -cfg.xClip), cfg.xClip);
    end
    vals(r) = classify_left_fraction(x, par);
end
Pmean = mean(vals, 'omitnan');
Pstd = std(vals, 0, 'omitnan');
end

%% ============================================================
