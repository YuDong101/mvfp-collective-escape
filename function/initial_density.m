function rho = initial_density(x, par, P0, initStd, cfg)
if nargin < 5 || ~isstruct(cfg)
    cfg = struct();
end

% An explicitly prepared density is shared across intervention parameters.
% The new nodal control-volume path uses these same trapezoidal weights.
if isfield(cfg, 'rho0') && ~isempty(cfg.rho0)
    if ~isfield(cfg, 'rho0X') || ~isnumeric(cfg.rho0X) || ...
            ~isreal(cfg.rho0X) || numel(cfg.rho0X) ~= numel(x) || ...
            any(~isfinite(cfg.rho0X(:))) || any(abs(cfg.rho0X(:)-x(:)) > 1e-12)
        error('MVFP:InitialDensity:GridMismatch', ...
            'Explicit rho0 requires the same finite rho0X grid.');
    end
    rho = cfg.rho0(:);
    if ~isnumeric(rho) || ~isreal(rho) || numel(rho) ~= numel(x) || ...
            any(~isfinite(rho)) || any(rho < 0)
        error('MVFP:InitialDensity:InvalidDensity', ...
            'rho0 must be a finite nonnegative density vector on x.');
    end
    mass0 = trapz(x, rho);
    if ~isfinite(mass0) || abs(mass0-1) > 1e-10
        error('MVFP:InitialDensity:InvalidMass', ...
            'Explicit rho0 must already have unit trapezoidal mass.');
    end
    return;
end

info0 = frozenP_initial_info(P0, par);
[xL, xR] = initial_representatives(info0, par);
variant = initial_density_variant(cfg);

switch lower(variant.type)
    case {'two-packet-gaussian', 'baseline', 'width'}
        rho = two_packet_density(x, P0, xL, xR, initStd);

    case {'shifted-packets', 'inward-shifted', 'outward-shifted'}
        [xLShift, xRShift] = shifted_packet_centers(P0, par, xL, xR, variant);
        rho = two_packet_density(x, P0, xLShift, xRShift, initStd);

    case {'shoulder-near-saddle', 'saddle-shoulder'}
        rho = shoulder_density(x, par, P0, xL, xR, initStd, variant);

    otherwise
        error('Unknown initial density variant type: %s', variant.type);
end
rho = normalize_density(x, rho);
end

%% ============================================================
function variant = initial_density_variant(cfg)
variant = struct('type', 'two-packet-gaussian');
if isfield(cfg, 'initialDensityVariant') && isstruct(cfg.initialDensityVariant)
    variant = merge_struct_fields(variant, cfg.initialDensityVariant);
end
if ~isfield(variant, 'type') || isempty(variant.type)
    variant.type = 'two-packet-gaussian';
end
end

%% ============================================================
function rho = two_packet_density(x, P0, xL, xR, initStd)
rho = P0 .* gaussian_pdf(x, xL, initStd) ...
    + (1 - P0) .* gaussian_pdf(x, xR, initStd);
end

%% ============================================================
function [xLShift, xRShift] = shifted_packet_centers(P0, par, xL, xR, variant)
shiftFraction = variant_value(variant, 'shiftFraction', 0.35);
direction = variant_value(variant, 'direction', 'inward');
xs = initial_saddle(P0, par);
if ~isfinite(xs)
    xs = 0;
end

if strcmpi(direction, 'outward') || strcmpi(variant.type, 'outward-shifted')
    signFactor = -1;
else
    signFactor = 1;
end
xLShift = xL + signFactor .* shiftFraction .* (xs - xL);
xRShift = xR + signFactor .* shiftFraction .* (xs - xR);
end

%% ============================================================
function rho = shoulder_density(x, par, P0, xL, xR, initStd, variant)
xs = initial_saddle(P0, par);
if ~isfinite(xs)
    xs = 0;
end
shoulderWeight = variant_value(variant, 'shoulderWeight', 0.2);
shoulderWeight = min(max(shoulderWeight, 0), 0.9);
shoulderFraction = variant_value(variant, 'shoulderFraction', 0.75);
shoulderStdScale = variant_value(variant, 'shoulderStdScale', 1.5);
shoulderStd = max(initStd .* shoulderStdScale, eps);

xLShoulder = xL + shoulderFraction .* (xs - xL);
xRShoulder = xR + shoulderFraction .* (xs - xR);

rhoLeft = (1 - shoulderWeight) .* gaussian_pdf(x, xL, initStd) ...
    + shoulderWeight .* gaussian_pdf(x, xLShoulder, shoulderStd);
rhoRight = (1 - shoulderWeight) .* gaussian_pdf(x, xR, initStd) ...
    + shoulderWeight .* gaussian_pdf(x, xRShoulder, shoulderStd);
rho = P0 .* rhoLeft + (1 - P0) .* rhoRight;
end

%% ============================================================
function xs = initial_saddle(P0, par)
F0 = initial_field_guess(P0, par);
[ok, ~, xs, ~] = extrema_from_F(F0, par);
if ~ok
    xs = nan;
end
end

%% ============================================================
function value = variant_value(variant, fieldName, defaultValue)
if isfield(variant, fieldName) && ~isempty(variant.(fieldName))
    value = variant.(fieldName);
else
    value = defaultValue;
end
end

%% ============================================================
