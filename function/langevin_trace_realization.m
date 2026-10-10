function Ptrace = langevin_trace_realization(par, P0, cfg, t, nSteps, dt, ...
    recordSteps, xLinit, xRinit, seed)
rng(seed, 'twister');
nL = round(P0 .* cfg.N);
x = [xLinit + 0.02 .* randn(nL, 1); ...
    xRinit + 0.02 .* randn(cfg.N - nL, 1)];
x = x(randperm(cfg.N));
Ptrace = nan(1, numel(t));
rec = 1;
Ptrace(rec) = classify_left_fraction(x, par);
for step = 1:nSteps
    M1 = mean(x);
    M2 = mean(x.^2);
    drift = (par.u - par.sigma1) .* x - (par.v + par.sigma2) .* x.^3 ...
        + par.sigma1 .* M1 + par.sigma2 .* M2 .* M1 + par.bias;
    x = x + drift .* dt + sqrt(2 .* par.D .* dt) .* randn(size(x));
    x = min(max(x, -cfg.xClip), cfg.xClip);
    if rec < numel(recordSteps) && step == recordSteps(rec + 1)
        rec = rec + 1;
        Ptrace(rec) = classify_left_fraction(x, par);
    end
end
end

%% ============================================================
