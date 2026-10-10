function diagnostics = compute_mvfp_feedback_control_diagnostics(par, opts)
fprintf('    feedback controls: full MVFP versus frozen feedback and no coupling\n');
P0Vals = numeric_option(opts, 'feedbackControlP0List', [0.060, 0.0647, 0.070, 0.080]);
times = numeric_option(opts, 'feedbackControlTimes', linspace(0, opts.T, 201));
readoutThreshold = scalar_option(opts, 'readoutThreshold', 0.5);
readoutLevels = numeric_option(opts, 'dynamicFeedbackReadoutLevels', 0.2:0.2:0.8);
P0Vals = unique(P0Vals(isfinite(P0Vals) & P0Vals >= 0 & P0Vals <= 1));
times = unique([0, times(isfinite(times)), opts.T]);

for k = 1:numel(P0Vals)
    P0 = P0Vals(k);
    caseOpts = opts;
    if option_enabled(opts, 'feedbackControlsCommonRho0', false)
        grid = mvfp_grid_cache(opts.mvfp);
        caseOpts.mvfp.rho0X = grid.x;
        caseOpts.mvfp.rho0 = initial_density( ...
            grid.x, par, P0, opts.mvfp.initStd, opts.mvfp);
        initialProtocol = 'common-rho0';
    else
        initialProtocol = 'model-specific-density';
    end
    fullCase = full_mvfp_case(par, caseOpts, P0, times, readoutLevels, readoutThreshold);
    frozenCase = frozen_feedback_case(par, caseOpts, P0, times, readoutLevels, readoutThreshold, ...
        fullCase.F(1));
    noCouplingCase = no_coupling_case(par, caseOpts, P0, times, readoutLevels, readoutThreshold);

    rhoFull0 = fullCase.rho(:, 1);
    rhoFrozen0 = frozenCase.rho(:, 1);
    rhoNoCoupling0 = noCouplingCase.rho(:, 1);
    x = fullCase.x;
    initialL1 = struct( ...
        'fullVsFrozen', trapz(x, abs(rhoFull0 - rhoFrozen0)), ...
        'fullVsNoCoupling', trapz(x, abs(rhoFull0 - rhoNoCoupling0)), ...
        'frozenVsNoCoupling', trapz(x, abs(rhoFrozen0 - rhoNoCoupling0)));
    bFull0 = (par.u - par.sigma1) .* x - (par.v + par.sigma2) .* x.^3 + fullCase.F(1);
    bFrozen0 = (par.u - par.sigma1) .* x - (par.v + par.sigma2) .* x.^3 + frozenCase.Ffixed;
    initialDriftL1FullVsFrozen = trapz(x, abs(bFull0 - bFrozen0) .* rhoFull0);

    caseOut = struct('P0', P0, 'full', fullCase, ...
        'frozenFeedback', frozenCase, 'noCoupling', noCouplingCase, ...
        'initialProtocol', initialProtocol, 'initialL1', initialL1, ...
        'initialDriftL1FullVsFrozen', initialDriftL1FullVsFrozen, ...
        'deltaFinalPFixedFrozenMinusFull', ...
            frozenCase.finalPFixed - fullCase.finalPFixed, ...
        'deltaFinalPFixedNoCouplingMinusFull', ...
            noCouplingCase.finalPFixed - fullCase.finalPFixed);
    if k == 1
        cases = caseOut;
    else
        cases(k) = caseOut; %#ok<AGROW>
    end
end

diagnostics = struct('method', 'full-vs-frozen-feedback-vs-no-coupling', ...
    'P0', P0Vals, 't', times, 'readoutThreshold', readoutThreshold, ...
    'readoutLevels', readoutLevels, 'cases', cases);
end

%% ============================================================
function out = full_mvfp_case(par, opts, P0, times, readoutLevels, readoutThreshold)
trace = solve_mvfp_density_trace(par, P0, times, opts.mvfp);
out = package_trace_case('full-MVFP', trace, par, readoutLevels, readoutThreshold);
end

%% ============================================================
function out = frozen_feedback_case(par, opts, P0, times, readoutLevels, readoutThreshold, Ffixed)
trace = solve_mvfp_density_trace_fixed_feedback(par, P0, times, opts.mvfp, Ffixed);
out = package_trace_case('frozen-feedback', trace, par, readoutLevels, readoutThreshold, Ffixed);
out.Ffixed = Ffixed;
end

%% ============================================================
function out = no_coupling_case(par, opts, P0, times, readoutLevels, readoutThreshold)
p = par;
p.sigma1 = 0;
p.sigma2 = 0;
trace = solve_mvfp_density_trace(p, P0, times, opts.mvfp);
out = package_trace_case('no-coupling', trace, p, readoutLevels, readoutThreshold);
end

%% ============================================================
function out = package_trace_case(label, trace, par, readoutLevels, readoutThreshold, Foverride)
if nargin < 6
    Foverride = [];
end
readouts = mvfp_trace_moment_readouts(trace, par);
if ~isempty(Foverride)
    readouts.F = Foverride .* ones(size(readouts.F));
end
landscape = mvfp_trace_dynamic_landscape(trace, readouts, par, readoutLevels);
transitionTime = mvfp_readout_arrival_times( ...
    trace.t, readouts.PFixed, readoutThreshold);

out = struct('label', label, 'x', trace.x, 't', trace.t, ...
    'rho', trace.rho, 'PSaddle', saddle_readout(trace, landscape.xs), ...
    'PFixed', readouts.PFixed, 'M1', readouts.M1, ...
    'M2', readouts.M2, 'F', readouts.F, ...
    'xL', landscape.xL, 'xs', landscape.xs, 'xR', landscape.xR, ...
    'xSingle', landscape.xSingle, ...
    'barrierRL', landscape.barrierRL, 'barrierLR', landscape.barrierLR, ...
    'barrierRLOverD', landscape.barrierRLOverD, ...
    'barrierLROverD', landscape.barrierLROverD, ...
    'readoutLevels', landscape.readoutLevels, ...
    'arrivalTimes', landscape.arrivalTimes, ...
    'transitionTime', transitionTime, ...
    'finalPFixed', readouts.PFixed(end), ...
    'finalPSaddle', saddle_readout_at(trace.x, trace.rho(:, end), landscape.xs(end)), ...
    'readoutDefinition', struct('PFixed', 'x<=0', ...
        'PSaddle', 'x<=instantaneous-saddle'), ...
    'massError', abs(trace.rawMass(:).' - 1), ...
    'boundaryMass', trace.boundaryMass(:).');
if isfield(trace, 'rhoRaw')
    out.rhoRaw = trace.rhoRaw;
end
if isfield(trace, 'rawMinDensity')
    out.rawMinDensity = trace.rawMinDensity(:).';
end
end

%% ============================================================
function values = saddle_readout(trace, xs)
values = nan(size(trace.t(:).'));
for k = 1:numel(values)
    values(k) = saddle_readout_at(trace.x, trace.rho(:, k), xs(k));
end
end

%% ============================================================
function value = saddle_readout_at(x, rho, xs)
if ~isfinite(xs)
    value = nan;
    return;
end
value = integrate_density_left_of_threshold(x, rho, xs);
end

%% ============================================================
function value = scalar_option(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
value = value(1);
end

%% ============================================================
function values = numeric_option(opts, name, defaultValues)
if isfield(opts, name) && ~isempty(opts.(name))
    values = opts.(name)(:).';
else
    values = defaultValues;
end
end

%% ============================================================
