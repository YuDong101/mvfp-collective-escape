function out = readout_curve_transition_summary(x, y, threshold)
if nargin < 3 || isempty(threshold)
    threshold = 0.5;
end
x = x(:).';
y = y(:).';
valid = isfinite(x) & isfinite(y);
xv = x(valid);
yv = y(valid);

out = struct('threshold', threshold, 'thresholdX', nan, ...
    'maxSlopeX', nan, 'maxAbsSlope', nan, 'sharpness', nan, ...
    'maxSlopeIndexLower', nan, 'maxSlopeIndexUpper', nan, ...
    'thresholdState', 'unresolved', 'jumpX', nan, 'jump', nan, ...
    'jumpIndexLower', nan, 'jumpIndexUpper', nan, ...
    'isSaturated', false, 'saturationState', 'none', ...
    'minY', nan, 'maxY', nan, 'n', numel(xv));
if numel(xv) < 1
    return;
end

out.minY = min(yv);
out.maxY = max(yv);
[out.isSaturated, out.saturationState] = saturation_state(yv);
out.thresholdX = threshold_crossing_x(xv, yv, threshold);
if isnan(out.thresholdX)
    if all(yv >= threshold)
        out.thresholdState = 'above';
        out.thresholdX = xv(1);
    elseif all(yv < threshold)
        out.thresholdState = 'below';
    end
else
    out.thresholdState = 'crossing';
end

if numel(xv) >= 2
    [out.jump, localIdx] = max(abs(diff(yv)));
    validIdx = find(valid);
    out.jumpIndexLower = validIdx(localIdx);
    out.jumpIndexUpper = validIdx(localIdx + 1);
    out.jumpX = 0.5 .* (x(out.jumpIndexLower) + x(out.jumpIndexUpper));
    dx = diff(xv);
    dy = diff(yv);
    slope = dy ./ dx;
    slope(~isfinite(slope) | dx == 0) = nan;
    [out.maxAbsSlope, slopeIdx] = finite_max_with_index(abs(slope));
    out.sharpness = out.maxAbsSlope;
    if isfinite(out.maxAbsSlope)
        out.maxSlopeX = 0.5 .* (xv(slopeIdx) + xv(slopeIdx + 1));
        validIdx = find(valid);
        out.maxSlopeIndexLower = validIdx(slopeIdx);
        out.maxSlopeIndexUpper = validIdx(slopeIdx + 1);
    end
end
end

%% ============================================================
function [tf, state] = saturation_state(y)
tol = 1e-6;
tf = false;
state = 'none';
if all(y >= 1 - tol)
    tf = true;
    state = 'upper';
elseif all(y <= tol)
    tf = true;
    state = 'lower';
end
end

%% ============================================================
function xCross = threshold_crossing_x(x, y, threshold)
xCross = nan;
if isempty(x)
    return;
end
exactIdx = find(y == threshold, 1, 'first');
if ~isempty(exactIdx)
    xCross = x(exactIdx);
    return;
end
offset = y - threshold;
crossIdx = find(offset(1:end-1) .* offset(2:end) < 0, 1, 'first');
if isempty(crossIdx)
    return;
end
x1 = x(crossIdx);
x2 = x(crossIdx + 1);
y1 = y(crossIdx);
y2 = y(crossIdx + 1);
if y2 == y1
    xCross = 0.5 .* (x1 + x2);
else
    xCross = x1 + (threshold - y1) .* (x2 - x1) ./ (y2 - y1);
end
end

%% ============================================================
function [m, idx] = finite_max_with_index(values)
idx = nan;
m = nan;
finiteMask = isfinite(values);
if ~any(finiteMask)
    return;
end
finiteIdx = find(finiteMask);
[m, localIdx] = max(values(finiteMask));
idx = finiteIdx(localIdx);
end

%% ============================================================
