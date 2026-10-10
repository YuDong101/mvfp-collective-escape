function readouts = mvfp_trace_moment_readouts(trace, par)
x = trace.x(:);
t = trace.t(:).';
nTimes = numel(t);
M1 = nan(1, nTimes);
M2 = nan(1, nTimes);
F = nan(1, nTimes);
PFixed = nan(1, nTimes);

for k = 1:nTimes
    rho = trace.rho(:, k);
    M1(k) = trapz(x, x .* rho);
    M2(k) = trapz(x, (x .^ 2) .* rho);
    F(k) = par.sigma1 .* M1(k) + par.sigma2 .* M2(k) .* M1(k) + par.bias;
    left = x <= 0;
    PFixed(k) = min(max(trapz(x(left), rho(left)), 0), 1);
end

readouts = struct('t', t, 'M1', M1, 'M2', M2, 'F', F, ...
    'PFixed', PFixed);
end

%% ============================================================
