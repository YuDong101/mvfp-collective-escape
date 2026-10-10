function diagnostics = compute_mvfp_free_energy_basin_diagnostics(par, opts)
%COMPUTE_MVFP_FREE_ENERGY_BASIN_DIAGNOSTICS Fixed-functional basin analysis.

p = basin_parameters(par, opts);
if p.sigma2 ~= 0
    error('MVFP:FreeEnergy:UnsupportedSigma2', ...
        'Basin diagnostics currently support only sigma2 = 0.');
end

times = basin_times(opts);
P0Values = basin_p0_values(opts);
x = linspace(-opts.mvfp.xMax, opts.mvfp.xMax, opts.mvfp.xGridN).';
branches = stationary_branch_diagnostics(p, opts, x);

cases = repmat(empty_basin_case(), 1, numel(P0Values));
for k = 1:numel(P0Values)
    cases(k) = trajectory_basin_case( ...
        p, opts, P0Values(k), times, branches);
end

criticalIndex = closest_unstable_case(cases);
dissipationTimes = basin_dissipation_times(opts, times);
if isequal(dissipationTimes, times)
    dissipationCase = cases(criticalIndex);
else
    dissipationCase = trajectory_basin_case( ...
        p, opts, P0Values(criticalIndex), dissipationTimes, branches);
end

selected = [cases.selectedStableIndex];
converged = logical([cases.convergedToStableBranch]);
PFixedFinal = arrayfun(@(c) c.PFixed(end), cases);
readoutThreshold = get_option(opts, 'readoutThreshold', 0.5);

diagnostics = struct('method', 'fixed-functional-basin-selection', ...
    'fixedFunctional', true, 'par', p, 'x', x, 't', times, ...
    'P0', P0Values, 'branches', branches, 'cases', cases, ...
    'dissipationCase', dissipationCase, ...
    'dissipationCaseIndex', criticalIndex, ...
    'nStableBranches', nnz(branches.stable), ...
    'nUnstableBranches', nnz(~branches.stable), ...
    'unstableResidenceDistanceTol', get_option(opts, ...
        'freeEnergyBasinUnstableDistanceTol', 0.1), ...
    'selectedStableIndex', selected, 'convergedMask', converged, ...
    'asymptoticBoundary', basin_boundary_summary( ...
        P0Values, selected, converged), ...
    'finiteTimeReadoutSummary', readout_curve_transition_summary( ...
        P0Values, PFixedFinal, readoutThreshold));
if isfield(opts, 'freeEnergyBasinPbInterval') && numel(opts.freeEnergyBasinPbInterval) == 2
    diagnostics.PbInterval = sort(opts.freeEnergyBasinPbInterval(:).');
elseif ~isempty(diagnostics.asymptoticBoundary.intervals)
    diagnostics.PbInterval = diagnostics.asymptoticBoundary.intervals(1, :);
else
    diagnostics.PbInterval = [nan nan];
end
end


%% ============================================================
function times = basin_dissipation_times(opts, fallback)
if isfield(opts, 'freeEnergyBasinDissipationTimes') ...
        && ~isempty(opts.freeEnergyBasinDissipationTimes)
    times = opts.freeEnergyBasinDissipationTimes(:).';
    finalTime = max(fallback);
    times = unique(times(isfinite(times) & times >= 0 ...
        & times <= finalTime), 'stable');
    if isempty(times)
        times = fallback;
    else
        times = unique([0, times, finalTime], 'stable');
    end
else
    times = fallback;
end
end

%% ============================================================
function index = closest_unstable_case(cases)
closestTime = [cases.closestUnstableTime];
[~, index] = max(closestTime, [], 'omitnan');
if isempty(index) || ~isfinite(index)
    distance = [cases.minDistanceToUnstable];
    [~, index] = min(distance, [], 'omitnan');
    if isempty(index) || ~isfinite(index)
        index = 1;
    end
end
end

%% ============================================================
function p = basin_parameters(par, opts)
p = par;
if isfield(opts, 'freeEnergyBasinD') && ~isempty(opts.freeEnergyBasinD)
    p.D = opts.freeEnergyBasinD;
end
if isfield(opts, 'freeEnergyBasinSigma1') ...
        && ~isempty(opts.freeEnergyBasinSigma1)
    p.sigma1 = opts.freeEnergyBasinSigma1;
end
end

%% ============================================================
function times = basin_times(opts)
if isfield(opts, 'freeEnergyBasinTimes') && ~isempty(opts.freeEnergyBasinTimes)
    times = opts.freeEnergyBasinTimes(:).';
else
    times = linspace(0, opts.T, 201);
end
valid = isfinite(times) & times >= 0;
hasFinalTime = isfield(opts, 'T') && isscalar(opts.T) ...
    && isfinite(opts.T) && opts.T >= 0;
if hasFinalTime
    valid = valid & times <= opts.T;
end
times = unique(times(valid), 'stable');
if isempty(times)
    error('MVFP:FreeEnergy:BasinTimes', ...
        'At least one nonnegative basin diagnostic time is required.');
end
if hasFinalTime
    times = unique([0, times, opts.T], 'stable');
elseif times(1) ~= 0
    times = [0, times];
end
end

%% ============================================================
function values = basin_p0_values(opts)
if isfield(opts, 'freeEnergyBasinP0List') ...
        && ~isempty(opts.freeEnergyBasinP0List)
    values = opts.freeEnergyBasinP0List(:).';
else
    values = opts.P0Scan(:).';
end
values = unique(values(isfinite(values) & values >= 0 & values <= 1), 'stable');
if isempty(values)
    error('MVFP:FreeEnergy:BasinP0', ...
        'At least one P0 value in [0,1] is required.');
end
end

%% ============================================================
function branches = stationary_branch_diagnostics(par, opts, x)
branches = compute_mvfp_steady_branches_for_param(par, opts, nan);
branches.stable = logical(branches.stable);
nBranches = numel(branches.F);
nX = numel(x);
rho = nan(nX, nBranches);
freeEnergyTotal = nan(1, nBranches);
freeEnergyEntropy = nan(1, nBranches);
freeEnergyLocal = nan(1, nBranches);
freeEnergyInteraction = nan(1, nBranches);
stationaryDissipation = nan(1, nBranches);
for j = 1:nBranches
    rhoJ = mvfp_steady_density_from_F(branches.F(j), par, x);
    if isempty(rhoJ)
        continue;
    end
    rho(:, j) = rhoJ;
    free = compute_mvfp_free_energy_trace(x, rhoJ, 0, par);
    freeEnergyTotal(j) = free.total;
    freeEnergyEntropy(j) = free.entropy;
    freeEnergyLocal(j) = free.localPotential;
    freeEnergyInteraction(j) = free.interaction;
    stationaryDissipation(j) = free.currentDissipation;
end
branches.x = x;
branches.rho = rho;
branches.freeEnergyTotal = freeEnergyTotal;
branches.freeEnergyEntropy = freeEnergyEntropy;
branches.freeEnergyLocal = freeEnergyLocal;
branches.freeEnergyInteraction = freeEnergyInteraction;
branches.stationaryDissipation = stationaryDissipation;
branches.globalStableIndex = global_stable_index( ...
    freeEnergyTotal, branches.stable);
end

%% ============================================================
function idx = global_stable_index(freeEnergy, stable)
idx = nan;
stableIdx = find(stable & isfinite(freeEnergy));
if isempty(stableIdx)
    return;
end
[~, local] = min(freeEnergy(stableIdx));
idx = stableIdx(local);
end

%% ============================================================
function c = trajectory_basin_case(par, opts, P0, times, branches)
trace = solve_mvfp_density_trace(par, P0, times, opts.mvfp);
x = trace.x(:);
nT = numel(times);
nBranches = numel(branches.F);
rhoNormalized = nan(size(trace.rho));
PFixed = nan(1, nT);

for k = 1:nT
    rhoK = trace.rho(:, k);
    mass = trapz(x, rhoK);
    if isfinite(mass) && mass > 0
        rhoK = rhoK ./ mass;
    end
    rhoNormalized(:, k) = rhoK;
    left = x <= 0;
    PFixed(k) = min(max(trapz(x(left), rhoK(left)), 0), 1);
end

free = compute_mvfp_free_energy_trace(x, rhoNormalized, times, par);
F = par.sigma1 .* free.M1 + par.bias;
distanceL1 = nan(nBranches, nT);
for j = 1:nBranches
    if all(isfinite(branches.rho(:, j)))
        distanceL1(j, :) = trapz(x, ...
            abs(rhoNormalized - branches.rho(:, j)), 1);
    end
end

stableIdx = find(branches.stable);
[densityDistance, densityIndex] = nearest_branch( ...
    distanceL1(:, end), stableIdx);
[feedbackDistance, feedbackIndex] = nearest_branch( ...
    abs(F(end) - branches.F(:)), stableIdx);
agreement = isfinite(densityIndex) && densityIndex == feedbackIndex;
densityTol = get_option(opts, 'freeEnergyBasinDensityTol', 0.05);
feedbackTol = get_option(opts, 'freeEnergyBasinFeedbackTol', 0.02);
converged = agreement && densityDistance <= densityTol ...
    && feedbackDistance <= feedbackTol;
selectedIndex = nan;
selectedF = nan;
if converged
    selectedIndex = densityIndex;
    selectedF = branches.F(selectedIndex);
end

[minUnstableDistance, closestUnstableIndex, closestUnstableTime] = ...
    closest_unstable_approach(distanceL1, times, branches.stable);
[distanceToClosestUnstable, unstableResidenceTime, unstableResidenceMask] = ...
    unstable_residence(distanceL1, times, branches.stable, ...
        get_option(opts, 'freeEnergyBasinUnstableDistanceTol', 0.1));

c = struct('P0', P0, 't', times, 'trace', trace, ...
    'PFixed', PFixed, 'M1', free.M1, 'M2', free.M2, 'F', F, ...
    'freeEnergy', free, 'distanceToBranchL1', distanceL1, ...
    'densitySelectedIndex', densityIndex, ...
    'feedbackSelectedIndex', feedbackIndex, ...
    'densityDistanceToSelected', densityDistance, ...
    'feedbackDistanceToSelected', feedbackDistance, ...
    'classificationAgreement', agreement, ...
    'convergedToStableBranch', converged, ...
    'selectedStableIndex', selectedIndex, ...
    'selectedStableF', selectedF, ...
    'minDistanceToUnstable', minUnstableDistance, ...
    'closestUnstableIndex', closestUnstableIndex, ...
    'closestUnstableTime', closestUnstableTime, ...
    'distanceToClosestUnstable', distanceToClosestUnstable, ...
    'unstableResidenceTime', unstableResidenceTime, ...
    'unstableResidenceMask', unstableResidenceMask, ...
    'massError', abs(trace.rawMass(:).' - 1), ...
    'boundaryMass', reshape(trace.boundaryMass, 1, []));
end

%% ============================================================
function [distance, residenceTime, residenceMask] = unstable_residence( ...
        distanceL1, times, stable, tolerance)
distance = nan(1, numel(times));
residenceTime = nan;
residenceMask = false(1, numel(times));
unstableIdx = find(~stable);
if isempty(unstableIdx)
    return;
end
distance = min(distanceL1(unstableIdx, :), [], 1, 'omitnan');
residenceMask = isfinite(distance) & distance <= tolerance;
if numel(times) == 1
    residenceTime = 0;
else
    residenceTime = trapz(times, double(residenceMask));
end
end

%% ============================================================
function [distance, index] = nearest_branch(allDistance, allowedIndices)
distance = nan;
index = nan;
if isempty(allowedIndices)
    return;
end
values = allDistance(allowedIndices);
valid = isfinite(values);
if ~any(valid)
    return;
end
validIndices = allowedIndices(valid);
[distance, local] = min(values(valid));
index = validIndices(local);
end

%% ============================================================
function [distance, branchIndex, timeValue] = closest_unstable_approach( ...
        distanceL1, times, stable)
distance = nan;
branchIndex = nan;
timeValue = nan;
unstableIdx = find(~stable);
if isempty(unstableIdx)
    return;
end
values = distanceL1(unstableIdx, :);
[distance, linearIndex] = min(values(:), [], 'omitnan');
if isempty(linearIndex) || ~isfinite(distance)
    distance = nan;
    return;
end
[row, timeIndex] = ind2sub(size(values), linearIndex);
branchIndex = unstableIdx(row);
timeValue = times(timeIndex);
end

%% ============================================================
function summary = basin_boundary_summary(P0, selected, converged)
[P0Sorted, order] = sort(P0);
selected = selected(order);
converged = converged(order);
intervals = nan(0, 2);
branchPairs = nan(0, 2);
for k = 1:(numel(P0Sorted) - 1)
    if converged(k) && converged(k + 1) ...
            && selected(k) ~= selected(k + 1)
        intervals(end + 1, :) = P0Sorted(k:k + 1); %#ok<AGROW>
        branchPairs(end + 1, :) = selected(k:k + 1); %#ok<AGROW>
    end
end
summary = struct('method', 'adjacent-converged-branch-change', ...
    'intervals', intervals, 'midpoints', mean(intervals, 2).', ...
    'branchPairs', branchPairs, 'nBoundaries', size(intervals, 1));
end

%% ============================================================
function value = get_option(opts, fieldName, fallback)
if isfield(opts, fieldName) && ~isempty(opts.(fieldName))
    value = opts.(fieldName);
else
    value = fallback;
end
end

%% ============================================================
function c = empty_basin_case()
c = struct('P0', nan, 't', [], 'trace', [], 'PFixed', [], ...
    'M1', [], 'M2', [], 'F', [], 'freeEnergy', [], ...
    'distanceToBranchL1', [], 'densitySelectedIndex', nan, ...
    'feedbackSelectedIndex', nan, 'densityDistanceToSelected', nan, ...
    'feedbackDistanceToSelected', nan, 'classificationAgreement', false, ...
    'convergedToStableBranch', false, 'selectedStableIndex', nan, ...
    'selectedStableF', nan, 'minDistanceToUnstable', nan, ...
    'closestUnstableIndex', nan, 'closestUnstableTime', nan, ...
    'distanceToClosestUnstable', [], 'unstableResidenceTime', nan, ...
    'unstableResidenceMask', [], ...
    'massError', [], 'boundaryMass', []);
end
