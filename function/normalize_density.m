function rho = normalize_density(x, rho)
rho = max(rho(:), 0);
mass = trapz(x, rho);
if mass > 0 && isfinite(mass)
    rho = rho ./ mass;
    return;
end
[~, idx] = min(abs(x));
rho(:) = 0;
rho(idx) = 1;
rho = rho ./ trapz(x, rho);
end

%% ============================================================
