function y = exp_if_representable(logValue)
if ~isfinite(logValue)
    y = nan;
elseif logValue > log(realmax)
    y = inf;
else
    y = exp(logValue);
end
end

%% ============================================================
