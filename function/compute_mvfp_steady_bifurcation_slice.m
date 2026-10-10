function out = compute_mvfp_steady_bifurcation_slice(par, opts, fieldName, values, label)
values = values(:).';
slices = cell(size(values));
if parallel_loop_enabled(opts)
    parfor k = 1:numel(values)
        p = set_param_value(par, fieldName, values(k));
        slices{k} = compute_mvfp_steady_branches_for_param(p, opts, values(k));
    end
else
    for k = 1:numel(values)
        p = set_param_value(par, fieldName, values(k));
        slices{k} = compute_mvfp_steady_branches_for_param(p, opts, values(k));
    end
end

rootCount = zeros(size(values));
branchValueCell = cell(size(values));
FCell = cell(size(values));
m1Cell = cell(size(values));
m2Cell = cell(size(values));
PssCell = cell(size(values));
stableCell = cell(size(values));
residualCell = cell(size(values));
for k = 1:numel(values)
    s = slices{k};
    rootCount(k) = numel(s.F);
    branchValueCell{k} = repmat(values(k), 1, numel(s.F));
    FCell{k} = s.F;
    m1Cell{k} = s.m1;
    m2Cell{k} = s.m2;
    PssCell{k} = s.Pss;
    stableCell{k} = s.stable;
    residualCell{k} = s.residualAtRoot;
end
branchValue = [branchValueCell{:}];
F = [FCell{:}];
m1 = [m1Cell{:}];
m2 = [m2Cell{:}];
Pss = [PssCell{:}];
stable = [stableCell{:}];
residualAtRoot = [residualCell{:}];

out = struct('fieldName', fieldName, 'label', label, 'values', values, ...
    'method', 'steady-self-consistent', 'rootCount', rootCount, ...
    'branchValue', branchValue, 'F', F, 'm1', m1, 'm2', m2, ...
    'Pss', Pss, 'stable', stable, 'residualAtRoot', residualAtRoot);
end

%% ============================================================
