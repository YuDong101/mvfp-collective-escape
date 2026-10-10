function bif = compute_mvfp_x_bifurcation_vs_sigma(par, opts, steadySlice)
if nargin < 3 || isempty(steadySlice)
    sigmaVals = x_sigma_bifurcation_scan(par, opts);
    steadySlice = compute_mvfp_steady_bifurcation_slice(par, opts, ...
        'sigma1', sigmaVals, '\sigma_1');
end

bif = x_sigma_extrema_from_steady_slice(par, steadySlice, ...
    'steady-self-consistent-x-branch');
bif.frozenInitial = compute_initial_frozen_x_sigma(par, opts);
if can_compute_focused_x_sigma(opts)
    focusVals = x_sigma_focus_scan(opts);
    focusSlice = compute_mvfp_steady_bifurcation_slice(par, opts, ...
        'sigma1', focusVals, '\sigma_1');
    bif.focused = x_sigma_extrema_from_steady_slice(par, focusSlice, ...
        'focused-steady-self-consistent-x-branch');
    bif.focused.sigma1Range = [min(focusVals), max(focusVals)];
    bif.focused.nRequested = numel(focusVals);
end
end

%% ============================================================
function out = x_sigma_extrema_from_steady_slice(par, steadySlice, method)
branchValue = steadySlice.branchValue(:).';
F = steadySlice.F(:).';
xL = nan(size(branchValue));
xs = nan(size(branchValue));
xR = nan(size(branchValue));
xSingle = nan(size(branchValue));
A = nan(size(branchValue));

for k = 1:numel(branchValue)
    p = par;
    p.sigma1 = branchValue(k);
    A(k) = p.u - p.sigma1;
    extrema = effective_potential_extrema_from_F(F(k), p);
    xL(k) = extrema.xL;
    xs(k) = extrema.xs;
    xR(k) = extrema.xR;
    xSingle(k) = extrema.xSingle;
end

out = struct('fieldName', 'sigma1', 'label', '\sigma_1', ...
    'values', steadySlice.values(:).', ...
    'rootCount', steadySlice.rootCount(:).', ...
    'branchValue', branchValue, 'F', F, ...
    'xBranch', steadySlice.m1(:).', ...
    'm1', steadySlice.m1(:).', 'm2', steadySlice.m2(:).', ...
    'Pss', steadySlice.Pss(:).', ...
    'stable', logical(steadySlice.stable(:).'), ...
    'residualAtRoot', steadySlice.residualAtRoot(:).', ...
    'A', A, 'u', par.u, 'xL', xL, 'xs', xs, 'xR', xR, ...
    'xSingle', xSingle, 'method', method);
end

%% ============================================================
function frozen = compute_initial_frozen_x_sigma(par, opts)
sigmaVals = x_sigma_bifurcation_scan(par, opts);
xL = nan(size(sigmaVals));
xs = nan(size(sigmaVals));
xR = nan(size(sigmaVals));
xSingle = nan(size(sigmaVals));
F = nan(size(sigmaVals));
A = nan(size(sigmaVals));

for k = 1:numel(sigmaVals)
    p = par;
    p.sigma1 = sigmaVals(k);
    F(k) = initial_field_guess(opts.P0, p);
    A(k) = p.u - p.sigma1;
    extrema = effective_potential_extrema_from_F(F(k), p);
    xL(k) = extrema.xL;
    xs(k) = extrema.xs;
    xR(k) = extrema.xR;
    xSingle(k) = extrema.xSingle;
end

frozen = struct('fieldName', 'sigma1', 'label', '\sigma_1', ...
    'P0', opts.P0, 'sigma1', sigmaVals, 'values', sigmaVals, ...
    'F', F, 'A', A, 'u', par.u, 'xL', xL, 'xs', xs, ...
    'xR', xR, 'xSingle', xSingle, ...
    'method', 'initial-frozen-effective-potential');
end

%% ============================================================
function tf = can_compute_focused_x_sigma(opts)
tf = isfield(opts, 'mvfp') && isfield(opts.mvfp, 'xMax') ...
    && isfield(opts.mvfp, 'xGridN') && isfield(opts, 'steadyFGrid') ...
    && isfield(opts, 'steadyRootTol');
end

%% ============================================================
function sigmaVals = x_sigma_focus_scan(opts)
if isfield(opts, 'xSigmaBifurcationFocusScan')
    sigmaVals = opts.xSigmaBifurcationFocusScan(:).';
elseif isfield(opts, 'xSigmaFocusedScan')
    sigmaVals = opts.xSigmaFocusedScan(:).';
else
    sigmaVals = linspace(1, 4, 41);
end
sigmaVals = unique(sigmaVals(isfinite(sigmaVals) & sigmaVals > 0));
if isempty(sigmaVals)
    error('focused x-sigma bifurcation scan must contain positive finite sigma1 values.');
end
end

%% ============================================================
function sigmaVals = x_sigma_bifurcation_scan(par, opts)
if isfield(opts, 'xSigmaBifurcationScan')
    sigmaVals = opts.xSigmaBifurcationScan(:).';
elseif isfield(opts, 'sigma1BifurcationScan')
    localVals = linspace(0.5 .* par.u, 1.5 .* par.u, 81);
    sigmaVals = [opts.sigma1BifurcationScan(:).', localVals];
else
    sigmaVals = [10.^linspace(-2, 1, 64), linspace(0.5 .* par.u, 1.5 .* par.u, 81)];
end
sigmaVals = unique(sigmaVals(isfinite(sigmaVals) & sigmaVals > 0));
if isempty(sigmaVals)
    error('x-sigma bifurcation scan must contain positive finite sigma1 values.');
end
end

%% ============================================================
