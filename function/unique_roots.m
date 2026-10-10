function rootsOut = unique_roots(rootsIn, tol)
rootsIn = sort(rootsIn(isfinite(rootsIn)));
nRoots = 0;
rootsWork = nan(size(rootsIn));
for k = 1:numel(rootsIn)
    if nRoots == 0 || abs(rootsIn(k) - rootsWork(nRoots)) > tol
        nRoots = nRoots + 1;
        rootsWork(nRoots) = rootsIn(k);
    end
end
rootsOut = rootsWork(1:nRoots);
end

%% ============================================================
