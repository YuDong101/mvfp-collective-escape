function [ok, xL, xs, xR] = extrema_from_F(F, par)
A = par.u - par.sigma1;
B = par.v + par.sigma2;
ok = 0;
xL = nan;
xs = nan;
xR = nan;
if ~(isfinite(A) && isfinite(B) && isfinite(F)) || B <= 0
    return;
end
r = roots([B, 0, -A, -F]);
r = sort(real(r(abs(imag(r)) < 1e-10)));
if numel(r) < 3
    return;
end
curv = Ueff2_fc(r, par);
if curv(1) > 0 && curv(2) < 0 && curv(3) > 0
    ok = 1;
    xL = r(1);
    xs = r(2);
    xR = r(3);
end
end

%% ============================================================
