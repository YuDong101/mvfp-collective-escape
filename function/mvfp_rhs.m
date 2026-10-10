function drho = mvfp_rhs(xOrGrid, rho, par)
if isstruct(xOrGrid)
    grid = xOrGrid;
else
    grid = mvfp_grid_cache_from_x(xOrGrid);
end

x = grid.x;
rho = rho(:);
mass = grid.w.' * rho;
if mass > 0 && isfinite(mass)
    M1 = (grid.wx.' * rho) ./ mass;
    M2 = (grid.wx2.' * rho) ./ mass;
else
    M1 = grid.wx.' * rho;
    M2 = grid.wx2.' * rho;
end

dx = grid.dx;
if isfield(grid,'fluxScheme') && strcmp(grid.fluxScheme,'exponential-fit')
    F = par.sigma1 .* M1 + par.sigma2 .* M2 .* M1 + par.bias;
    J = mvfp_exponential_fit_flux(x,rho,par,F);
else
    b = (par.u - par.sigma1) .* x - (par.v + par.sigma2) .* grid.x3 ...
        + par.sigma1 .* M1 + par.sigma2 .* M2 .* M1 + par.bias;
    bFace = 0.5 .* (b(1:end-1) + b(2:end));
    rhoLeft = rho(1:end-1);
    rhoRight = rho(2:end);
    rhoUpwind = rhoLeft;
    rhoUpwind(bFace < 0) = rhoRight(bFace < 0);
    Jadv = bFace .* rhoUpwind;
    Jdiff = -par.D .* (rhoRight - rhoLeft) ./ dx;
    J = Jadv + Jdiff;
end
Jfull = [0; J; 0];
if isfield(grid,'volumeWeights')
    drho = -(Jfull(2:end) - Jfull(1:end-1)) ./ grid.volumeWeights;
else
    drho = -(Jfull(2:end) - Jfull(1:end-1)) ./ dx;
end
end

%% ============================================================
