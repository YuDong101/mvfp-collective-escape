function logIntegral = log_trapz_exp(x, logf)
x = x(:);
logf = logf(:);
if numel(x) < 2
    logIntegral = -inf;
else
    logCum = log_cumtrapz_exp(x, logf);
    logIntegral = logCum(end);
end
end

%% ============================================================
