function landscape = mvfp_trace_dynamic_landscape(trace, readouts, par, levels)
if nargin < 4 || isempty(levels)
    levels = 0.2:0.2:0.8;
end

t = trace.t(:).';
F = readouts.F(:).';
nTimes = numel(t);
xL = nan(1, nTimes);
xs = nan(1, nTimes);
xR = nan(1, nTimes);
xSingle = nan(1, nTimes);
barrierLR = nan(1, nTimes);
barrierRL = nan(1, nTimes);

for k = 1:nTimes
    extrema = effective_potential_extrema_from_F(F(k), par);
    xL(k) = extrema.xL;
    xs(k) = extrema.xs;
    xR(k) = extrema.xR;
    xSingle(k) = extrema.xSingle;
    if all(isfinite([xL(k), xs(k), xR(k)]))
        U = effective_potential_from_F([xL(k), xs(k), xR(k)], par, F(k));
        barrierLR(k) = U(2) - U(1);
        barrierRL(k) = U(2) - U(3);
    end
end

arrivalTimes = mvfp_readout_arrival_times(t, readouts.PFixed, levels);
landscape = struct('xL', xL, 'xs', xs, 'xR', xR, ...
    'xSingle', xSingle, 'barrierLR', barrierLR, ...
    'barrierRL', barrierRL, ...
    'barrierLROverD', barrierLR ./ par.D, ...
    'barrierRLOverD', barrierRL ./ par.D, ...
    'readoutLevels', levels(:).', ...
    'arrivalTimes', arrivalTimes);
end

%% ============================================================
