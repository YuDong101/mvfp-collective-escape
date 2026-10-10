function [t,rho]=mvfp_bordered_integrate(g,p,times,rho0,cfg,Ffixed)
% Exact index-one reformulation of fitted/nodal pairwise mean-field dynamics.
% F is an algebraic variable, leaving a sparse bordered Newton system.
% The spatial flux, initial state, boundary conditions and tolerances match
% the original eliminated-mean-field ODE. No clipping is applied to raw rho.
n=numel(g.x);
if ~isempty(Ffixed)
    q=p;q.bias=Ffixed-p.sigma1*(g.wx.'*rho0)/(g.w.'*rho0);
    lin=mvfp_revision_fitted_linearization(g,rho0,q);
    op=odeset('RelTol',cfg.odeRelTol,'AbsTol',cfg.odeAbsTol,'NonNegative',1:n,'Jacobian',lin.A);
    [t,rho]=ode15s(@(~,r)mvfp_rhs_fixed_feedback(g,r,p,Ffixed),times,rho0,op);
    return
end
F0=p.sigma1*(g.wx.'*rho0)/(g.w.'*rho0)+p.bias;
mass=spdiags([ones(n,1);0],0,n+1,n+1);
op=odeset('RelTol',cfg.odeRelTol,'AbsTol',cfg.odeAbsTol, ...
    'Mass',mass,'MassSingular','yes','MStateDependence','none','Jacobian',@(~,y) jac(g,y,p));
[t,y]=ode15s(@(~,y) rhs(g,y,p),times,[rho0;F0],op);
rho=y(:,1:n);
constraint=y(:,end)-p.sigma1*(rho*g.wx)./(rho*g.w)-p.bias;
assert(max(abs(constraint))<max(1e-8,100*cfg.odeAbsTol),'Algebraic feedback constraint failed');
assert(min(rho(:))>=-max(1e-10,10*cfg.odeAbsTol),'Negative raw density exceeds declared tolerance');
end

function dy=rhs(g,y,p)
r=y(1:end-1);F=y(end);M=(g.wx.'*r)/(g.w.'*r);
dy=[mvfp_rhs_fixed_feedback(g,r,p,F);p.sigma1*M+p.bias-F];
end

function J=jac(g,y,p)
r=y(1:end-1);F=y(end);q=p;q.bias=F-p.sigma1*(g.wx.'*r)/(g.w.'*r);
lin=mvfp_revision_fitted_linearization(g,r,q);
J=[lin.A,sparse(lin.b);sparse(lin.c'),-1];
end
