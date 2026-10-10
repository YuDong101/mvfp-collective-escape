function free = compute_mvfp_free_energy_trace(x, rho, t, par)
%COMPUTE_MVFP_FREE_ENERGY_TRACE Pairwise MVFP free-energy post-processing.

if ~isfield(par, 'sigma2') || ~isfinite(par.sigma2) || par.sigma2 ~= 0
    error('MVFP:FreeEnergy:UnsupportedSigma2', ...
        'The free-energy diagnostic currently supports only sigma2 = 0.');
end

x = x(:);
if isvector(rho)
    rho = rho(:);
end
t = t(:).';
if size(rho, 1) ~= numel(x)
    error('MVFP:FreeEnergy:DensitySize', ...
        'rho must have one row per spatial grid point.');
end
if size(rho, 2) ~= numel(t)
    error('MVFP:FreeEnergy:TimeSize', ...
        'rho must have one column per requested time.');
end

nT = numel(t);
mass = nan(1, nT);
M1 = nan(1, nT);
M2 = nan(1, nT);
entropy = nan(1, nT);
localPotential = nan(1, nT);
interaction = nan(1, nT);
currentDissipation = nan(1, nT);
rhsFreeEnergyDecay = nan(1, nT);
V0 = -0.5 .* par.u .* x.^2 + 0.25 .* par.v .* x.^4 ...
    - par.bias .* x;
dxFace = diff(x);
% Match the production nodal-control-volume discretization when evaluating
% the numerical RHS used in the energy-rate check.  Calling mvfp_rhs with a
% bare x vector would silently fall back to the legacy uniform divergence.
rhsGrid = mvfp_grid_cache_from_x(x);
rhsGrid.volumeWeights = rhsGrid.w;

for k = 1:nT
    rhoK = rho(:, k);
    mass(k) = trapz(x, rhoK);
    if ~(isfinite(mass(k)) && mass(k) > 0)
        continue;
    end
    rhoK = rhoK ./ mass(k);
    M1(k) = trapz(x, x .* rhoK);
    M2(k) = trapz(x, x.^2 .* rhoK);
    entropyIntegrand = zeros(size(rhoK));
    positive = rhoK > 0;
    entropyIntegrand(positive) = rhoK(positive) .* log(rhoK(positive));
    entropy(k) = par.D .* trapz(x, entropyIntegrand);
    localPotential(k) = trapz(x, V0 .* rhoK);
    interaction(k) = 0.5 .* par.sigma1 .* (M2(k) - M1(k).^2);
    chemicalPotential = par.D .* (log(max(rhoK, realmin('double'))) + 1) ...
        + V0 + 0.5 .* par.sigma1 .* ...
        (x.^2 - 2 .* M1(k) .* x + M2(k));
    rhoFace = max(0.5 .* (rhoK(1:end-1) + rhoK(2:end)), ...
        realmin('double'));
    current = -rhoFace .* diff(chemicalPotential) ./ dxFace;
    currentDissipation(k) = sum((current.^2 ./ rhoFace) .* dxFace);
    numericalRhs = mvfp_rhs(rhsGrid, rhoK, par);
    rhsFreeEnergyDecay(k) = -trapz(x, chemicalPotential .* numericalRhs);
end

total = entropy + localPotential + interaction;
if nT >= 2
    minusDGdt = -gradient(total, t);
else
    minusDGdt = nan(size(total));
end
dissipationResidual = minusDGdt - currentDissipation;
relativeScale = max([abs(minusDGdt); currentDissipation; ...
    eps .* ones(size(total))], [], 1);
rhsResidual = minusDGdt - rhsFreeEnergyDecay;
rhsScale = max([abs(minusDGdt); abs(rhsFreeEnergyDecay); ...
    eps .* ones(size(total))], [], 1);
comparisonMask = isfinite(minusDGdt) & isfinite(rhsFreeEnergyDecay) ...
    & minusDGdt > 0 & rhsFreeEnergyDecay > 0;
if nT >= 2
    comparisonMask([1 end]) = false;
end

free = struct('method', 'pairwise-mckean-vlasov-free-energy', ...
    't', t, 'mass', mass, 'massError', abs(mass - 1), ...
    'M1', M1, 'M2', M2, ...
    'entropy', entropy, 'localPotential', localPotential, ...
    'interaction', interaction, ...
    'total', total, 'totalChange', total - total(1), ...
    'currentDissipation', currentDissipation, ...
    'rhsFreeEnergyDecay', rhsFreeEnergyDecay, ...
    'minusDGdt', minusDGdt, ...
    'dissipationResidual', dissipationResidual, ...
    'relativeDissipationResidual', abs(dissipationResidual) ./ relativeScale, ...
    'rhsEnergyRateResidual', rhsResidual, ...
    'relativeRhsEnergyRateResidual', abs(rhsResidual) ./ rhsScale, ...
    'dissipationComparisonMask', comparisonMask, ...
    'currentMethod', 'continuum-chemical-potential-centered-face-current', ...
    'rhsRateMethod', 'minus-integral-mu-times-upwind-finite-volume-rhs');
end
