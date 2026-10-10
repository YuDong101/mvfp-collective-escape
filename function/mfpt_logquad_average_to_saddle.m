function [tauAvg, logTauAvg, localDx, nGrid] = mfpt_logquad_average_to_saddle( ...
    xLeft, xRight, par, F, domainGridN, absorbSide)
tauAvg = nan;
logTauAvg = nan;
localDx = nan;
nGrid = 0;

if ~(isfinite(xLeft) && isfinite(xRight) && xRight > xLeft)
    return;
end
if ~(isfinite(domainGridN) && domainGridN >= 5 && isfinite(par.D) && par.D > 0)
    return;
end

nGrid = max(5, round(domainGridN));
x = linspace(xLeft, xRight, nGrid).';
localDx = x(2) - x(1);
U = effective_potential_from_F(x, par, F);
logWeight = -U ./ par.D;
logZ = log_trapz_exp(x, logWeight);
if ~isfinite(logZ)
    return;
end

if strcmpi(absorbSide, 'right')
    logInner = log_cumtrapz_exp(x, logWeight);
    logOuterIntegrand = U ./ par.D + logInner;
    logTauProfile = log_cumtrapz_exp_reverse(x, logOuterIntegrand) - log(par.D);
elseif strcmpi(absorbSide, 'left')
    logInner = log_cumtrapz_exp_reverse(x, logWeight);
    logOuterIntegrand = U ./ par.D + logInner;
    logTauProfile = log_cumtrapz_exp(x, logOuterIntegrand) - log(par.D);
else
    return;
end

logNumerator = log_trapz_exp(x, logTauProfile + logWeight);
logTauAvg = logNumerator - logZ;
tauAvg = exp_if_representable(logTauAvg);
end

%% ============================================================
