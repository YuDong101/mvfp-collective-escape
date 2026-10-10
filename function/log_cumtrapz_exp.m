function logCum = log_cumtrapz_exp(x, logf)
x = x(:);
logf = logf(:);
n = numel(x);
logCum = -inf(n, 1);
acc = -inf;
for i = 1:(n - 1)
    dx = x(i + 1) - x(i);
    if ~(isfinite(dx) && dx > 0)
        continue;
    end
    term = log(dx ./ 2) + log_add_exp(logf(i), logf(i + 1));
    acc = log_add_exp(acc, term);
    logCum(i + 1) = acc;
end
end

%% ============================================================
