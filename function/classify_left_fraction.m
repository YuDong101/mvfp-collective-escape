function Pleft = classify_left_fraction(x, par) %#ok<INUSD>
finite = isfinite(x);
if ~any(finite)
    Pleft = nan;
    return;
end
x = x(finite);
Pleft = mean(x <= 0);
end

%% ============================================================
