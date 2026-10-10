function [tauLR, tauRL, logTauLR, logTauRL, localDxLR, localDxRL, ...
    nGridLR, nGridRL, barrierLR, barrierRL, xL, xs, xR] = ...
    compute_mvfp_mfpt_point_from_F(F, par, x, xMin, xMax, domainGridN)
tauLR = nan;
tauRL = nan;
logTauLR = nan;
logTauRL = nan;
localDxLR = nan;
localDxRL = nan;
nGridLR = 0;
nGridRL = 0;
barrierLR = nan;
barrierRL = nan;
xL = nan;
xs = nan;
xR = nan;

U = effective_potential_from_F(x, par, F);
[ok, xLk, xsk, xRk] = extrema_from_F(F, par);
if ~ok
    return;
end
xL = xLk;
xs = xsk;
xR = xRk;

if xsk <= xMin || xsk >= xMax
    return;
end

[tauLR, logTauLR, localDxLR, nGridLR] = ...
    mfpt_logquad_average_to_saddle(xMin, xsk, par, F, domainGridN, 'right');
[tauRL, logTauRL, localDxRL, nGridRL] = ...
    mfpt_logquad_average_to_saddle(xsk, xMax, par, F, domainGridN, 'left');

UL = interp1(x, U, xLk, 'linear', nan);
US = interp1(x, U, xsk, 'linear', nan);
UR = interp1(x, U, xRk, 'linear', nan);
barrierLR = US - UL;
barrierRL = US - UR;
end

%% ============================================================
