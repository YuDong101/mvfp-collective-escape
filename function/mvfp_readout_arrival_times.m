function arrivalTimes = mvfp_readout_arrival_times(t, PFixed, levels)
t = t(:).';
PFixed = PFixed(:).';
levels = levels(:).';
arrivalTimes = nan(size(levels));
if numel(t) ~= numel(PFixed)
    return;
end

for i = 1:numel(levels)
    level = levels(i);
    diffVals = PFixed - level;
    hitIdx = find(diffVals == 0, 1, 'first');
    if ~isempty(hitIdx)
        arrivalTimes(i) = t(hitIdx);
        continue;
    end

    crossIdx = find(diffVals(1:end-1) .* diffVals(2:end) < 0, 1, 'first');
    if isempty(crossIdx)
        continue;
    end
    p1 = PFixed(crossIdx);
    p2 = PFixed(crossIdx + 1);
    if p2 == p1
        arrivalTimes(i) = t(crossIdx);
    else
        arrivalTimes(i) = t(crossIdx) ...
            + (level - p1) .* (t(crossIdx + 1) - t(crossIdx)) ./ (p2 - p1);
    end
end
end

%% ============================================================
