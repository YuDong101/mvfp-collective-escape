function Pleft = integrate_density_left_of_threshold(x, rho, threshold)
if threshold <= x(1)
    Pleft = 0;
    return;
end
if threshold >= x(end)
    Pleft = 1;
    return;
end
xAug = [x(:); threshold];
rhoAtThreshold = interp1(x, rho, threshold, 'linear', nan);
rhoAug = [rho(:); rhoAtThreshold];
[xAug, idx] = sort(xAug);
rhoAug = rhoAug(idx);
left = xAug <= threshold;
Pleft = trapz(xAug(left), rhoAug(left));
Pleft = min(max(Pleft, 0), 1);
end

%% ============================================================
