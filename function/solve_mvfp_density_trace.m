function trace = solve_mvfp_density_trace(par, P0, times, cfg)
times = unique(times(:).');
if isempty(times)
    times = 0;
end
if any(times < 0)
    error('MVFP output times must be nonnegative.');
end

grid = mvfp_grid_cache(cfg);
x = grid.x;
rho0 = initial_density(x, par, P0, cfg.initStd, cfg);
tSolve = unique([0, times]);
if max(tSolve) <= 0
    tOut = 0;
    y = rho0.';
else
    odeOpts = odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol, ...
        'NonNegative', 1:numel(x));
    rhs = @(~, rho) mvfp_rhs(grid, rho, par);
    if isfield(grid,'fluxScheme') && strcmp(grid.fluxScheme,'exponential-fit') && par.sigma2==0
        [tOut,y]=mvfp_bordered_integrate(grid,par,tSolve,rho0,cfg,[]);
    else
        [tOut, y] = ode15s(rhs, tSolve, rho0, odeOpts);
    end
end

rhoOut = nan(numel(x), numel(times));
storeRaw = option_enabled(cfg,'storeRawDensity',false);
if storeRaw, rhoRaw = nan(numel(x),numel(times)); end
rawMinDensity = nan(numel(times),1);
rawBoundaryMass = nan(numel(times),1);
P = nan(numel(times), 1);
mass = nan(numel(times), 1);
rawMass = nan(numel(times), 1);
boundaryMass = nan(numel(times), 1);
threshold = nan(numel(times), 1);

for k = 1:numel(times)
    [~, idx] = min(abs(tOut - times(k)));
    rho = y(idx, :).';
    if storeRaw, rhoRaw(:,k) = rho; end
    rawMinDensity(k) = min(rho);
    rawBoundaryMass(k) = boundary_density_mass(x,rho);
    rawMass(k) = trapz(x, rho);
    rho = normalize_density(x, rho);
    threshold(k) = density_threshold(x, rho, par, cfg.thresholdMode);
    left = x <= threshold(k);
    P(k) = min(max(trapz(x(left), rho(left)), 0), 1);
    mass(k) = trapz(x, rho);
    boundaryMass(k) = boundary_density_mass(x, rho);
    rhoOut(:, k) = rho;
end

trace = struct('x', x, 't', times, 'rho', rhoOut, 'P', P, ...
    'mass', mass, 'rawMass', rawMass, 'boundaryMass', boundaryMass, ...
    'threshold', threshold);
trace.rawMinDensity = rawMinDensity;
trace.rawBoundaryMass = rawBoundaryMass;
trace.massDiscretization = 'legacy-uniform-divergence';
if isfield(grid,'volumeWeights')
    trace.massDiscretization = 'nodal-control-volume';
end
if storeRaw, trace.rhoRaw = rhoRaw; end
end

%% ============================================================
