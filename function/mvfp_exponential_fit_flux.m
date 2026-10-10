function J = mvfp_exponential_fit_flux(x, rho, par, F)
%MVFP_EXPONENTIAL_FIT_FLUX Exponentially fitted face flux for a fixed field.
if ~isscalar(par.D) || ~isfinite(par.D) || par.D <= 0
    error('MVFP:Flux:NonpositiveDiffusion', ...
        'The exponential-fit flux requires a finite scalar D > 0.');
end

x = x(:);
rho = rho(:);
dx = x(2) - x(1);
U = effective_potential_from_F(x,par,F);
z = diff(U) ./ par.D;
J = (par.D ./ dx) .* (bernoulli_stable(z) .* rho(1:end-1) ...
    - bernoulli_stable(-z) .* rho(2:end));
end

function B = bernoulli_stable(z)
B = zeros(size(z));
small = abs(z) < 1e-4;
zs = z(small);
B(small) = 1 - zs./2 + zs.^2./12 - zs.^4./720;

positive = z >= 50;
zp = z(positive);
ez = exp(-zp);
B(positive) = zp .* ez ./ (1 - ez);

negative = z <= -50;
zn = -z(negative);
B(negative) = zn ./ (1 - exp(-zn));

middle = ~(small | positive | negative);
B(middle) = z(middle) ./ expm1(z(middle));
end
