function y = log_add_exp(a, b)
if isinf(a) && a < 0
    y = b;
    return;
end
if isinf(b) && b < 0
    y = a;
    return;
end
m = max(a, b);
if isinf(m)
    y = m;
else
    y = m + log(exp(a - m) + exp(b - m));
end
end

%% ============================================================
