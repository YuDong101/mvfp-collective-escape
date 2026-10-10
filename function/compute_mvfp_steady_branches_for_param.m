function branches = compute_mvfp_steady_branches_for_param(par, opts, paramValue)
x = linspace(-opts.mvfp.xMax, opts.mvfp.xMax, opts.mvfp.xGridN).';
FGrid = opts.steadyFGrid(:).';
residualGrid = nan(size(FGrid));
for k = 1:numel(FGrid)
    residualGrid(k) = mvfp_steady_residual(FGrid(k), par, x);
end

tol = opts.steadyRootTol;
rootCandidates = nan(1, 2 .* numel(FGrid));
nRoots = 0;
nearMask = isfinite(residualGrid) & abs(residualGrid) <= tol;
nNear = nnz(nearMask);
if nNear > 0
    rootCandidates(1:nNear) = FGrid(nearMask);
    nRoots = nNear;
end
for k = 1:(numel(FGrid) - 1)
    r1 = residualGrid(k);
    r2 = residualGrid(k + 1);
    if ~(isfinite(r1) && isfinite(r2))
        continue;
    end
    if r1 == 0 || r2 == 0
        continue;
    end
    if r1 * r2 < 0
        try
            rootF = fzero(@(F) mvfp_steady_residual(F, par, x), ...
                [FGrid(k), FGrid(k + 1)]);
            nRoots = nRoots + 1;
            rootCandidates(nRoots) = rootF;
        catch
        end
    end
end
rootsF = rootCandidates(1:nRoots);
rootsF = unique_roots(rootsF, max(10 .* tol, 1e-7));

n = numel(rootsF);
F = nan(1, n);
m1 = nan(1, n);
m2 = nan(1, n);
Pss = nan(1, n);
PssFixed = nan(1, n);
PssSaddle = nan(1, n);
thresholdSaddle = nan(1, n);
hasSaddle = false(1, n);
chi = nan(1, n);
stable = zeros(1, n);
residualAtRoot = nan(1, n);
for k = 1:n
    F(k) = rootsF(k);
    [residualAtRoot(k), m1(k), m2(k), Pss(k)] = ...
        mvfp_steady_residual(F(k), par, x);
    [~, ~, ~, ~, ~, readout] = mvfp_steady_density_from_F(F(k), par, x);
    PssFixed(k) = readout.PssFixed;
    PssSaddle(k) = readout.PssSaddle;
    thresholdSaddle(k) = readout.thresholdSaddle;
    hasSaddle(k) = readout.hasSaddle;
    if par.sigma2 == 0 && isfinite(par.D) && par.D > 0
        chi(k) = par.sigma1 .* max(m2(k) - m1(k).^2, 0) ./ par.D;
    end
    stable(k) = steady_branch_stability(F(k), par, x);
end
branches = struct('paramValue', paramValue, 'F', F, 'm1', m1, ...
    'm2', m2, 'Pss', Pss, 'PssFixed', PssFixed, ...
    'PssSaddle', PssSaddle, 'thresholdFixed', 0, ...
    'thresholdSaddle', thresholdSaddle, 'hasSaddle', hasSaddle, ...
    'legacyPssDefinition', 'saddle-if-present-otherwise-fixed-zero', ...
    'fixedReadoutDefinition', 'integral-x-less-than-or-equal-to-zero', ...
    'saddleReadoutDefinition', 'integral-x-less-than-or-equal-to-saddle', ...
    'chi', chi, 'chiDefinition', 'sigma1-times-variance-over-D', ...
    'chiInterpretation', 'static-free-energy-second-variation-diagnostic-not-time-growth-rate', ...
    'stable', stable, 'residualAtRoot', residualAtRoot);
end

%% ============================================================
