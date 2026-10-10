function logCum = log_cumtrapz_exp_reverse(x, logf)
x = x(:);
logf = logf(:);
n = numel(x);
logCum = -inf(n, 1);
acc = -inf;
for i = (n - 1):-1:1
    dx = x(i + 1) - x(i);
    if ~(isfinite(dx) && dx > 0)
        continue;
    end
    term = log(dx ./ 2) + log_add_exp(logf(i), logf(i + 1));
    acc = log_add_exp(acc, term);
    logCum(i) = acc;
end
end

%% ============================================================
