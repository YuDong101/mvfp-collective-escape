function extrema = effective_potential_extrema_from_F(F, par)
A = par.u - par.sigma1;
B = par.v + par.sigma2;

extrema = struct('roots', [], 'curvature', [], ...
    'xL', nan, 'xs', nan, 'xR', nan, 'xSingle', nan);
if ~(isfinite(A) && isfinite(B) && isfinite(F)) || B <= 0
    return;
end

r = roots([B, 0, -A, -F]);
r = sort(real(r(abs(imag(r)) < 1e-10)));
if isempty(r)
    return;
end

curv = Ueff2_fc(r, par);
minima = r(curv > -1e-10);
saddles = r(curv < 1e-10);

extrema.roots = r(:).';
extrema.curvature = curv(:).';
if numel(minima) >= 2
    extrema.xL = minima(1);
    extrema.xR = minima(end);
    if ~isempty(saddles)
        extrema.xs = saddles(1);
    end
elseif numel(minima) == 1
    extrema.xSingle = minima(1);
end
end

%% ============================================================
