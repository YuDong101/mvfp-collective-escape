function mass = boundary_density_mass(x, rho)
n = numel(x);
if n < 4
    mass = trapz(x, rho);
    return;
end
nEdge = max(2, ceil(0.05 * n));
nEdge = min(nEdge, floor(n / 2));
left = 1:nEdge;
right = (n - nEdge + 1):n;
mass = trapz(x(left), rho(left)) + trapz(x(right), rho(right));
end

%% ============================================================
