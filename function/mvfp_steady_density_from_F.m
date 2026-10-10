function [rho, m1, m2, Pss, threshold, readout] = mvfp_steady_density_from_F(F, par, x)
readout = empty_readout();
U = effective_potential_from_F(x, par, F);
logWeight = -U ./ par.D;
shift = max(logWeight);
w = exp(logWeight - shift);
Z = trapz(x, w);
if ~(isfinite(Z) && Z > 0)
    rho = [];
    m1 = nan;
    m2 = nan;
    Pss = nan;
    threshold = nan;
    return;
end
rho = w ./ Z;
m1 = trapz(x, x .* rho);
m2 = trapz(x, x.^2 .* rho);
[ok, ~, xs, ~] = extrema_from_F(F, par);
if ok
    threshold = xs;
    readout.thresholdSaddle = xs;
    readout.PssSaddle = integrate_density_left_of_threshold(x, rho, xs);
else
    threshold = 0;
end
readout.PssFixed = integrate_density_left_of_threshold(x, rho, 0);
readout.legacyThreshold = threshold;
readout.hasSaddle = logical(ok);
if ok
    readout.legacyDefinition = 'saddle-threshold';
    Pss = readout.PssSaddle;
else
    readout.legacyDefinition = 'fixed-threshold-fallback';
    Pss = readout.PssFixed;
end
readout.Pss = Pss;
end

function readout = empty_readout()
readout = struct('Pss', nan, 'PssFixed', nan, 'PssSaddle', nan, ...
    'thresholdFixed', 0, 'thresholdSaddle', nan, ...
    'legacyThreshold', nan, 'hasSaddle', false, ...
    'fixedDefinition', 'integral-x-less-than-or-equal-to-zero', ...
    'saddleDefinition', 'integral-x-less-than-or-equal-to-saddle', ...
    'legacyDefinition', 'unavailable');
end

%% ============================================================
