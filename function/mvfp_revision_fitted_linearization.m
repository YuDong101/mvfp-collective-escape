function lin = mvfp_revision_fitted_linearization(grid,rho,par)
%MVFP_REVISION_FITTED_LINEARIZATION Full fitted-flux Jacobian A+b*c'.
% Includes the normalized moment denominator and both no-flux end faces.
% This opt-in helper deliberately rejects the legacy discretization.
if ~isfield(grid,'volumeWeights') || ~isfield(grid,'fluxScheme') || ...
        ~strcmp(grid.fluxScheme,'exponential-fit') || ...
        norm(grid.volumeWeights-grid.w,inf)>10*eps || par.sigma2~=0
    error('MVFP:RevisionSpectrum:Convention','Requires pairwise fitted flux and nodal volumes.');
end
x=grid.x; w=grid.w; n=numel(x); rho=rho(:); m=w'*rho;
if m<=0 || any(~isfinite(rho)), error('MVFP:RevisionSpectrum:Density','Invalid density.'); end
M=(w.*x)'*rho/m; F=par.sigma1*M+par.bias;
z=diff(effective_potential_from_F(x,par,F))/par.D;
[bp,dbp]=bernoulli(z); [bm,dbm]=bernoulli(-z);
face=sparse([(1:n-1)';(1:n-1)'],[(1:n-1)';(2:n)'], ...
    par.D/grid.dx*[bp;-bm],n-1,n);
div=sparse([(1:n-1)';(2:n)'],[(1:n-1)';(1:n-1)'], ...
    [-1./w(1:n-1);1./w(2:n)],n,n-1);
A=div*face;
% dz/dF=-dx/D: derivative of both Bernoulli factors is included.
b=div*(-dbp.*rho(1:n-1)-dbm.*rho(2:n));
c=par.sigma1*w.*(x-M)/m;
lin=struct('A',A,'b',b,'c',c,'F',F,'mass',m,'mean',M, ...
    'apply',@(v) A*v+b*(c'*v));
end

function [B,dB]=bernoulli(z)
B=zeros(size(z)); dB=B;
s=abs(z)<1e-3; t=z(s);
B(s)=1-t/2+t.^2/12-t.^4/720+t.^6/30240;
dB(s)=-1/2+t/6-t.^3/180+t.^5/5040;
p=z>=50; t=z(p); q=exp(-t);
B(p)=t.*q./(1-q); dB(p)=q./(1-q)-t.*q./(1-q).^2;
v=z<=-50; t=-z(v); q=exp(-t);
B(v)=t./(1-q); dB(v)=-1./(1-q)+t.*q./(1-q).^2;
k=~(s|p|v); t=z(k); q=expm1(t);
B(k)=t./q; dB(k)=(q-t.*exp(t))./q.^2;
end
