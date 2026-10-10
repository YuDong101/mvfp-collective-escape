function out = compute_mvfp_revision_local_spectrum(par,cfg)
%COMPUTE_MVFP_REVISION_LOCAL_SPECTRUM Opt-in exploratory B spectrum preflight.
% Positive secular root uses sparse tridiagonal solves, never a dense fine
% grid eigenproblem. A full weighted zero-mass spectrum is computed only
% when Nx <= denseCheckMaxNx. Fine grids verify a positive mode, not a full
% numerical count of positive directions; chi is reported separately.
if nargin<2, cfg=struct; end
cfg=defaults(cfg,struct('Nx',1001,'L',1,'rootBracket',[0.03 0.06], ...
    'denseCheckMaxNx',201,'runGrowth',false,'growthAmplitudes',[1e-6 1e-5], ...
    'growthRelTol',1e-9,'growthAbsTol',1e-12,'growthMaxStep',0.1));
x=linspace(-cfg.L,cfg.L,cfg.Nx)'; grid=mvfp_grid_cache_from_x(x);
grid.volumeWeights=grid.w; grid.fluxScheme='exponential-fit'; w=grid.w;
F=fzero(@(f) residual(f,par,x),cfg.rootBracket,optimset('TolX',1e-14));
[rho,M,M2]=mvfp_steady_density_from_F(F,par,x);
lin=mvfp_revision_fitted_linearization(grid,rho,par); n=numel(x); I=speye(n);
secular=@(s) 1-lin.c'*((s*I-lin.A)\lin.b);
lo=1e-6; hi=1;
while secular(hi)<0 && hi<1e5, hi=hi*2; end
if secular(lo)*secular(hi)>=0
    error('MVFP:RevisionSpectrum:NoBracket','Positive secular root not bracketed.');
end
lambda=fzero(secular,[lo hi],optimset('TolX',1e-12));
e=(lambda*I-lin.A)\lin.b; e=e/sum(w.*abs(e));
if (w.*x)'*e<0, e=-e; end
q=(lambda*I-lin.A')\lin.c; q=q/(q'*e); psi=q./w;
R=mvfp_rhs(grid,rho,par); Je=lin.apply(e);
leftResidual=lin.A'*q+lin.c*(lin.b'*q)-lambda*q;
% Smooth relative direction is positive for all finite difference samples.
v=rho.*(sin(3*x)+0.3*cos(5*x)); v=v-rho*(w'*v)/(w'*rho);
v=v/max(abs(v./max(rho,realmin)));
epsList=[1e-3 1e-4 1e-5]; fd=zeros(3,1); fdMass=fd;
for k=1:3
    h=epsList(k); d=(mvfp_rhs(grid,rho+h*v,par)-mvfp_rhs(grid,rho-h*v,par))/(2*h);
    fd(k)=sum(w.*abs(d-lin.apply(v)))/max(sum(w.*abs(lin.apply(v))),realmin);
    fdMass(k)=abs(w'*d);
end
out=struct('status','exploratory-preflight','cfg',cfg,'par',par,'grid',grid, ...
    'F',F,'rhoSteady',rho,'M1',M,'chi',par.sigma1*(M2-M*M)/par.D, ...
    'lambda',lambda,'rightMode',e,'leftMode',psi,'steadyRhs',R, ...
    'steadyResidualL1',sum(w.*abs(R)),'selfConsistencyResidual',residual(F,par,x), ...
    'massSteady',w'*rho,'massMode',w'*e,'modePairing',q'*e, ...
    'rightResidualL1',sum(w.*abs(Je-lambda*e)), ...
    'leftRelativeResidual',norm(leftResidual)/max(norm(lambda*q),realmin), ...
    'massJacobianA',norm(w'*lin.A,inf),'massJacobianB',abs(w'*lin.b), ...
    'residualSeed',abs(q'*R)/lambda,'fdAmplitudes',epsList,'fdRelativeErrors',fd, ...
    'fdMass',fdMass,'fdMinimumDensity',min(rho-1e-3*abs(v)), ...
    'numericalPositiveCount',NaN,'positiveCountVerified',false, ...
    'positiveModeVerified',false,'growth',[]);
out.positiveModeVerified=lambda>0 && out.rightResidualL1<1e-7 && out.leftRelativeResidual<1e-6;
if n<=cfg.denseCheckMaxNx
    % Euclidean orthonormal basis after sqrt(W) transform: w'*eta=0.
    sw=sqrt(w); Q=null(sw');
    J=full(lin.A)+lin.b*lin.c';
    J0=Q'*(sw.*J./sw')*Q;
    ev=eig(J0); [~,order]=sort(real(ev),'descend'); ev=ev(order);
    out.zeroMassEigenvalues=ev;
    out.numericalPositiveCount=sum(real(ev)>1e-7);
    out.positiveCountVerified=true;
    out.denseLeadingDifference=abs(real(ev(1))-lambda);
else
    out.zeroMassEigenvalues=[]; out.denseLeadingDifference=NaN;
end
if cfg.runGrowth
    times=linspace(0,min(5,1/lambda),41)';
    for k=1:numel(cfg.growthAmplitudes)
        a=cfg.growthAmplitudes(k);
        % Smooth exponential tilt preserves tiny positive tails. Limiting
        % |e/rho| guards roundoff in resolvent tails; no rho is clipped.
        ratio=e./max(rho,realmin); ratio=max(-1e4,min(1e4,ratio));
        rho0=rho.*exp(a*ratio); rho0=rho0/(w'*rho0);
        jac=@(~,r) fullJac(grid,r,par);
        op=odeset('RelTol',cfg.growthRelTol,'AbsTol',cfg.growthAbsTol, ...
            'MaxStep',cfg.growthMaxStep,'Jacobian',jac);
        tic; [t,y]=ode15s(@(~,r) mvfp_rhs(grid,r,par),times,rho0,op); elapsed=toc;
        amp=(y-rho')*q; predicted=amp(1)*exp(lambda*t);
        fit=polyfit(t,log(abs(amp)),1);
        g=struct('requestedAmplitude',a,'t',t,'rhoRaw',y', ...
            'projection',amp,'linearPrediction',predicted, ...
            'relativePredictionError',max(abs(amp-predicted)./abs(predicted)), ...
            'fittedGrowthRate',fit(1),'massRaw',y*w,'minRawDensity',min(y,[],2), ...
            'initialDirectionRelativeError',sum(w.*abs((rho0-rho)/amp(1)-e)), ...
            'residualSeedFraction',out.residualSeed/abs(amp(1)),'wallSeconds',elapsed);
        out.growth=[out.growth;g]; %#ok<AGROW>
    end
end
end

function r=residual(F,par,x)
[~,m]=mvfp_steady_density_from_F(F,par,x); r=par.sigma1*m+par.bias-F;
end
function J=fullJac(grid,rho,par)
lin=mvfp_revision_fitted_linearization(grid,rho,par); J=full(lin.A)+lin.b*lin.c';
end
function cfg=defaults(cfg,d)
names=fieldnames(d); for k=1:numel(names), if ~isfield(cfg,names{k}), cfg.(names{k})=d.(names{k}); end, end
end
