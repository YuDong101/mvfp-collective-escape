function [xLinit, xRinit] = initial_representatives(info, par)
if info.valid && all(isfinite([info.xLstar, info.xRstar]))
    xLinit = info.xLstar;
    xRinit = info.xRstar;
    return;
end

[ok, xL, ~, xR] = extrema_from_F(par.bias, par);
if ok
    xLinit = xL;
    xRinit = xR;
    return;
end

A = par.u - par.sigma1;
B = par.v + par.sigma2;
if A > 0 && B > 0
    x0 = sqrt(max(A ./ B, eps));
else
    x0 = 1;
end
xLinit = -x0;
xRinit = x0;
end

%% ============================================================
