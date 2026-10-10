function out = compute_mvfp_revision_residence(cfg)
%COMPUTE_MVFP_REVISION_RESIDENCE Case-B residence prediction feasibility.
% The unstable rate is loaded from an independent spectrum artifact and is
% never fitted to the residence trajectories.
if nargin<1 || ~isstruct(cfg), error('MVFP:RevisionResidence:ConfigRequired','A configuration struct is required.'); end
cfg=apply_defaults(cfg);
s=load_spectrum(cfg.spectrumFile);
[x,w,rhoU,eU,psi,lambda,par]=validate_spectrum(s,cfg);
validate_uncertainties(cfg,lambda);
artifacts=validate_artifacts(cfg,par,x);

center=mean(cfg.PbInterval); multipliers=[1 2 4 8 16 32];
P0=[center-cfg.deltaMin*multipliers center+cfg.deltaMin*multipliers];
side=[repmat("lower",1,6) repmat("upper",1,6)];
rank=[1:6 1:6]; isCalibration=ismember(rank,[2 4 6]);
design=struct('P0',num2cell(P0),'side',cellstr(side), ...
    'distanceRank',num2cell(rank),'distance',num2cell(abs(P0-center)), ...
    'isCalibration',num2cell(isCalibration),'isHoldout',num2cell(~isCalibration));
out=struct('status','prepared-no-physical-runs','createdAt',char(datetime('now')), ...
    'config',cfg,'par',par,'spectrumFile',cfg.spectrumFile, ...
    'spectrumHash',mvfp_streaming_sha256(cfg.spectrumFile), ...
    'lambda',lambda,'lambdaUncertainty',cfg.lambdaUncertainty, ...
    'PbInterval',cfg.PbInterval,'Pb',center,'basinUncertainty',cfg.basinUncertainty, ...
    'deltaMin',cfg.deltaMin,'offsetMultipliers',multipliers,'design',design, ...
    'readoutDefinition','PFixed(t)=integral over x<=0', ...
    'distanceDefinition','weighted L1 distance to the named unstable steady density');
out.sourceHashes=source_hashes();
out.artifacts=artifacts;
out.matlabVersion=version;
out.budgetLimitation='Checks occur before/after each case; the shared ODE solver has no hard wall-time interrupt.';
out.configHash=mvfp_value_sha256(strip_runtime(cfg));
if cfg.prepareOnly, return; end
if any(P0<0 | P0>1), error('MVFP:RevisionResidence:P0Range','Requested offsets leave [0,1].'); end

cases=repmat(empty_case(),1,0); started=tic;
if ~isempty(cfg.outputFile) && isfile(cfg.outputFile)
    previous=load(cfg.outputFile,'out');
    if isfield(previous,'out')
        old=previous.out;
        fields={'configHash','spectrumHash','artifacts','sourceHashes','matlabVersion'};
        compatible=all(isfield(old,fields));
        for j=1:numel(fields)
            compatible=compatible && isequaln(old.(fields{j}),out.(fields{j}));
        end
        if ~compatible
            error('MVFP:RevisionResidence:CheckpointMismatch','Use a new outputFile for changed physics, processing, or provenance.');
        end
        if isfield(old,'cases'), cases=old.cases; end
    end
end
out.cases=cases; out.status='clean-partial'; checkpoint(out,cfg);
for k=numel(cases)+1:numel(design)
    if toc(started)>=cfg.stageBudgetSeconds, out.status='clean-partial-budget'; checkpoint(out,cfg); return; end
    caseStarted=tic;
    solver=cfg.solver; solver.xGridN=numel(x); solver.xMax=max(abs(x));
    solver.rho0=[]; solver.rho0X=[]; solver.storeRawDensity=true;
    gridCheck=mvfp_grid_cache(solver);
    assert_grid(x,gridCheck.x,'solver');
    times=unique([0:cfg.traceStep:cfg.T cfg.T]);
    try
        trace=solve_mvfp_density_trace(par,design(k).P0,times,solver);
    catch err
        out.status='clean-partial-solver-failure'; out.failureIdentifier=err.identifier;
        out.failureMessage=err.message; checkpoint(out,cfg); rethrow(err);
    end
    assert_grid(x,trace.x,'trace');
    cases(k)=summarize_case(trace,design(k),x,w,rhoU,eU,psi,lambda,par,cfg,s);
    out.cases=cases; out.completedCases=k; checkpoint(out,cfg);
    if toc(caseStarted)>cfg.maxCaseSeconds, out.status='clean-partial-case-budget'; checkpoint(out,cfg); return; end
end
out.cases=cases;
out.prediction=calibrate_prediction(cases,lambda,cfg);
out.status=out.prediction.status;
if ~isempty(cfg.outputFile)
    folder=fileparts(cfg.outputFile); if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
    save(cfg.outputFile,'out','-v7.3');
end
end

function cfg=apply_defaults(cfg)
d=struct('T',200,'traceStep',.1,'entryRadius',.05,'exitRadii',[.08 .10 .12], ...
    'baseExitRadius',.10,'prepareOnly',false,'outputFile','', ...
    'densitySampleTimes',[0 25 50 100 150 200], ...
    'stageBudgetSeconds',900,'maxCaseSeconds',120, ...
    'quality',struct('minLocalSamples',5,'minLocalDuration',.5,'growthRmsMax',.2,'stableFractionMax',.2, ...
    'rawMassErrorMax',1e-6,'boundaryMassMax',1e-6,'negativeDensityMax',1e-10, ...
    'residualSeedFractionMax',.01,'heldoutRelativeErrorMax',.1,'predictionRelativeWidthMax',.2));
f=fieldnames(d); for k=1:numel(f), if ~isfield(cfg,f{k}), cfg.(f{k})=d.(f{k}); end, end
qfields=fieldnames(d.quality);
for k=1:numel(qfields)
    f=qfields{k};
    if ~isfield(cfg.quality,f) || ~isscalar(cfg.quality.(f)) || ~isfinite(cfg.quality.(f)) || cfg.quality.(f)<=0
        error('MVFP:RevisionResidence:InvalidQualityGate','All quality thresholds must be predeclared finite positive scalars.');
    end
end
if ~isscalar(cfg.T) || ~isfinite(cfg.T) || cfg.T<=0 || ~isscalar(cfg.traceStep) || ~isfinite(cfg.traceStep) || cfg.traceStep<=0 || ...
    ~isscalar(cfg.stageBudgetSeconds) || ~isfinite(cfg.stageBudgetSeconds) || cfg.stageBudgetSeconds<=0 || ...
    ~isscalar(cfg.maxCaseSeconds) || ~isfinite(cfg.maxCaseSeconds) || cfg.maxCaseSeconds<=0
    error('MVFP:RevisionResidence:InvalidBudget','Time/grid budgets must be positive finite scalars.');
end
required={'spectrumFile','PbInterval','basinUncertainty','lambdaUncertainty','PbGridN','basinFile','uncertaintyFile'};
for k=1:numel(required)
    if ~isfield(cfg,required{k}) || isempty(cfg.(required{k}))
        error('MVFP:RevisionResidence:MissingInput','cfg.%s must be supplied empirically.',required{k});
    end
end
if ~isfield(cfg,'deltaMin') || isempty(cfg.deltaMin), cfg.deltaMin=10*cfg.basinUncertainty; end
if ~isfield(cfg,'solver')
    cfg.solver=struct('xMax',1,'xGridN',[],'initStd',.02,'thresholdMode','fixed', ...
        'odeRelTol',1e-8,'odeAbsTol',1e-12,'storeRawDensity',true, ...
        'fluxScheme','exponential-fit','massDiscretization','nodal-control-volume');
end
end

function s=load_spectrum(file)
if ~isfile(file), error('MVFP:RevisionResidence:MissingSpectrum','Spectrum file does not exist.'); end
q=load(file);
if ~isfield(q,'s') || ~isstruct(q.s), error('MVFP:RevisionResidence:InvalidSpectrum','MAT file must contain struct variable s.'); end
s=q.s;
end

function [x,w,rho,e,psi,lambda,par]=validate_spectrum(s,cfg)
needed={'grid','rhoSteady','rightMode','leftMode','lambda','par'};
if ~all(isfield(s,needed)), error('MVFP:RevisionResidence:InvalidSpectrum','Spectrum lacks named grid/density/mode fields.'); end
x=s.grid.x(:); rho=s.rhoSteady(:); e=s.rightMode(:); psi=s.leftMode(:); lambda=s.lambda; par=s.par;
if ~isscalar(lambda) || ~isfinite(lambda) || ~isreal(lambda), error('MVFP:RevisionResidence:InvalidSpectrum','Invalid lambda.'); end
if isfield(s.grid,'volumeWeights'), w=s.grid.volumeWeights(:); elseif isfield(s.grid,'w'), w=s.grid.w(:); else, error('MVFP:RevisionResidence:InvalidSpectrum','Spectrum grid lacks weights.'); end
if any([numel(w) numel(rho) numel(e) numel(psi)]~=numel(x)) || any(~isfinite([x;w;rho;e;psi])) || ~isreal([x;w;rho;e;psi]) || lambda<=0
    error('MVFP:RevisionResidence:InvalidSpectrum','Spectrum arrays or positive lambda are invalid.');
end
if any(w<=0) || any(diff(x)<=0) || abs(w'*rho-1)>1e-8 || abs(w'*e)>1e-8 || abs((w.*psi)'*e-1)>1e-6
    error('MVFP:RevisionResidence:InvalidSpectrum','Density mass, zero-mass mode or left/right pairing invalid.');
end
if ~isfield(s.grid,'fluxScheme') || ~strcmp(s.grid.fluxScheme,'exponential-fit') || ~isfield(s.grid,'volumeWeights')
    error('MVFP:RevisionResidence:DiscretizationMismatch','Spectrum requires exponential-fit and nodal-control-volume.');
end
if ~isfield(cfg.solver,'fluxScheme') || ~strcmp(cfg.solver.fluxScheme,'exponential-fit') || ...
    ~isfield(cfg.solver,'massDiscretization') || ~strcmp(cfg.solver.massDiscretization,'nodal-control-volume') || ~strcmp(cfg.solver.thresholdMode,'fixed')
    error('MVFP:RevisionResidence:DiscretizationMismatch','Solver must match spectrum and fixed x=0 readout.');
end
check=cfg.solver; check.xGridN=numel(x); check.xMax=max(abs(x)); gg=mvfp_grid_cache(check);
if abs(cfg.solver.xMax-check.xMax)>1e-12
    error('MVFP:RevisionResidence:SpectrumGridMismatch','Solver domain differs from spectrum.');
end
assert_grid(x,gg.x,'solver'); assert_grid(w,gg.volumeWeights,'weights');
if isfield(cfg,'expectedGrid'), assert_grid(x,cfg.expectedGrid,'expected'); end
if isfield(cfg,'PbGrid'), assert_grid(x,cfg.PbGrid,'basin'); end
if ~isscalar(cfg.PbGridN) || cfg.PbGridN~=numel(x), error('MVFP:RevisionResidence:SpectrumGridMismatch','Pb and spectrum grid sizes differ.'); end
end

function validate_uncertainties(cfg,lambda)
if ~isnumeric(cfg.PbInterval) || numel(cfg.PbInterval)~=2 || any(~isfinite(cfg.PbInterval)) || cfg.PbInterval(2)<cfg.PbInterval(1)
    error('MVFP:RevisionResidence:InvalidPbInterval','PbInterval must be a finite ordered pair.');
end
if ~isscalar(cfg.basinUncertainty) || ~isfinite(cfg.basinUncertainty) || cfg.basinUncertainty<=0 || ...
        ~isscalar(cfg.lambdaUncertainty) || ~isfinite(cfg.lambdaUncertainty) || cfg.lambdaUncertainty<0
    error('MVFP:RevisionResidence:InvalidUncertainty','Empirical basin and lambda uncertainties are required.');
end
if cfg.basinUncertainty < diff(cfg.PbInterval)/2
    error('MVFP:RevisionResidence:UncertaintyGate','basinUncertainty must cover interval halfwidth.');
end
if ~isscalar(cfg.deltaMin) || ~isfinite(cfg.deltaMin) || cfg.deltaMin<=0 || cfg.deltaMin < 10*cfg.basinUncertainty
    error('MVFP:RevisionResidence:UncertaintyGate','deltaMin must be at least ten times basinUncertainty.');
end
if cfg.lambdaUncertainty>=lambda, error('MVFP:RevisionResidence:InvalidUncertainty','lambda uncertainty must be smaller than lambda.'); end
end

function c=summarize_case(trace,design,x,w,rhoU,eU,psi,lambda,par,cfg,spectrum)
rho=trace.rhoRaw; eta=rho-rhoU; l1=sum(w.*abs(eta),1).';
a=(w.*psi)'*eta; stable=eta-eU*a; stableNorm=sum(w.*abs(stable),1).';
M1=((w.*x)'*rho)./(w'*rho); F=par.sigma1*M1+par.bias;
raw=trace.rhoRaw; Praw=trapz(x(x<=0),raw(x<=0,:),1).';
events=mvfp_revision_residence_events(trace.t,l1,cfg.exitRadii,cfg.entryRadius);
sampleTimes=cfg.densitySampleTimes(cfg.densitySampleTimes<=cfg.T);
[~,sampleIndex]=arrayfun(@(z)min(abs(trace.t-z)),sampleTimes);
seed=Inf; if isfield(spectrum,'residualSeed') && isscalar(spectrum.residualSeed), seed=abs(spectrum.residualSeed); end
windowCfg=cfg.quality; windowCfg.entryRadius=cfg.entryRadius;
window=mvfp_revision_outgoing_modal_window(trace.t,l1,a,stableNorm./max(l1,realmin),seed,windowCfg);
local=false(size(l1)); local(window.indices)=true;
growthResidual=nan(size(a(:)));
if nnz(local)>=2
    av=a(:); tv=trace.t(:); da=gradient(av(local),tv(local));
    growthResidual(local)=(da-lambda*av(local))./max(abs(lambda*av(local)),realmin);
end
initialRho=rho(:,1); initialActualP0=trapz(x(x<=0),initialRho(x<=0));
c=empty_case(); c.P0Requested=design.P0; c.side=design.side; c.distanceRank=design.distanceRank;
c.distance=design.distance; c.isCalibration=design.isCalibration; c.isHoldout=design.isHoldout;
c.t=trace.t(:); c.distanceL1=l1; c.projection=a(:); c.stableComponentL1=stableNorm;
c.F=F(:); c.PFixed=trace.P(:); c.Praw=Praw; c.mass=trace.mass(:); c.rawMass=trace.rawMass(:);
c.rawMinDensity=trace.rawMinDensity(:); c.rawBoundaryMass=trace.rawBoundaryMass(:);
c.initialRho=initialRho; c.initialActualP0=initialActualP0; c.initialM1=M1(1); c.initialF=F(1);
c.t50=mvfp_revision_first_arrival(trace.t,trace.P,.5);
c.densitySampleTimes=trace.t(sampleIndex); c.densitySamples=rho(:,sampleIndex);
c.rawDensitySamples=raw(:,sampleIndex); c.events=events; c.localProjectionGrowthResidual=growthResidual;
c.outgoingModalWindow=window;
c.localProjectionGrowthRms=sqrt(mean(growthResidual(local).^2,'omitnan'));
q=cfg.quality; c.qualityThresholds=q;
c.localStableFraction=max(stableNorm(local)./max(l1(local),realmin),[],'omitnan');
c.residualSeedFraction=Inf;
if isfield(spectrum,'residualSeed') && isscalar(spectrum.residualSeed) && any(local), c.residualSeedFraction=abs(spectrum.residualSeed)/min(abs(a(local))); end
c.qualityPassed=nnz(local)>=q.minLocalSamples && all(isfinite(growthResidual(local))) && ...
    c.localProjectionGrowthRms<=q.growthRmsMax && ~isempty(c.localStableFraction) && c.localStableFraction<=q.stableFractionMax && ...
    all(isfinite([c.rawMass;c.rawMinDensity;c.rawBoundaryMass])) && ...
    max(abs(c.rawMass-1))<=q.rawMassErrorMax && max(c.rawBoundaryMass)<=q.boundaryMassMax && ...
    min(c.rawMinDensity)>=-q.negativeDensityMax && isfinite(c.residualSeedFraction) && ...
    c.residualSeedFraction<=q.residualSeedFractionMax;
c.validationStatus='exploratory-local-gates-failed';
if c.qualityPassed, c.validationStatus='local-gates-passed-review-required'; end
end

function p=calibrate_prediction(cases,lambda,cfg)
[~,base]=min(abs(cfg.exitRadii-cfg.baseExitRadius));
rows=struct([]); sideNames={"lower","upper"}; allEnough=true;
for sidx=1:2
    ids=find(strcmp({cases.side},sideNames{sidx})); usable=false(size(ids)); tres=nan(size(ids)); eu=nan(size(ids));
    for j=1:numel(ids)
        e=cases(ids(j)).events; usable(j)=e.status(base)=="exited" && cases(ids(j)).qualityPassed; tres(j)=e.tRes(base); eu(j)=e.eventTimeUncertainty(base);
    end
    cal=usable & [cases(ids).isCalibration]; hold=usable & [cases(ids).isHoldout];
    enough=nnz(cal)==3 && nnz(hold)>=2; allEnough=allEnough && enough;
    if any(cal) %#ok<ALIGN>
        d=[cases(ids(cal)).distance]; Cvals=tres(cal)+log(d)/lambda;
        C=mean(Cvals);
    else, C=nan; end
    for j=1:numel(ids)
        d=cases(ids(j)).distance; pred=C-log(d)/lambda;
        bounds=[NaN NaN]; predUnc=NaN;
        if any(cal)
            bounds=mvfp_revision_residence_prediction_bounds([cases(ids(cal)).P0Requested], ...
                cases(ids(j)).P0Requested,tres(cal),eu(cal), ...
                mean(cfg.PbInterval)+[-1 1]*cfg.basinUncertainty,lambda+[-1 1]*cfg.lambdaUncertainty);
            predUnc=max(abs(bounds-pred));
        end
        err=tres(j)-pred; rel=abs(err)/max(abs(tres(j)),realmin);
        errorBounds=[tres(j)-eu(j)-bounds(2),tres(j)+eu(j)-bounds(1)];
        totalUnc=predUnc+eu(j);
        rows=[rows struct('side',sideNames{sidx},'distanceRank',cases(ids(j)).distanceRank, ...
            'P0',cases(ids(j)).P0Requested,'distance',d,'role',ternary(cases(ids(j)).isCalibration,"calibration","holdout"), ...
            'usable',usable(j),'eventStatus',cases(ids(j)).events.status(base), ...
            'measuredResidence',tres(j),'eventTimeUncertainty',eu(j), ...
            'predictedResidence',pred,'predictionUncertainty',predUnc, ...
            'error',err,'relativeError',rel,'predictionLower',bounds(1),'predictionUpper',bounds(2), ...
            'errorLower',errorBounds(1),'errorUpper',errorBounds(2),'totalErrorUncertainty',totalUnc, ...
            'heldoutPass',~cases(ids(j)).isCalibration && usable(j) && enough && ...
            rel<=cfg.quality.heldoutRelativeErrorMax && totalUnc/max(abs(pred),realmin)<=cfg.quality.predictionRelativeWidthMax)]; %#ok<AGROW>
    end
end
p=struct('fixedSlope',-1/lambda,'baseExitRadius',cfg.exitRadii(base), ...
    'table',struct2table(rows),'uncertaintyMethod', ...
    'conservative shared-Pb/shared-lambda prediction interval plus calibration event brackets; total half-width gate also includes measured event uncertainty');
if allEnough, p.status='feasibility-computed-review-heldout-errors'; else, p.status='insufficient-prespecified-usable-points'; end
end

function y=ternary(q,a,b), if q, y=a; else, y=b; end, end
function assert_grid(a,b,label)
b=b(:); if numel(a)~=numel(b) || any(abs(a-b)>1e-12*max(1,max(abs(a)))) %#ok<ALIGN>
    error('MVFP:RevisionResidence:SpectrumGridMismatch','%s grid differs from the named spectrum grid.',label);
end
end
function h=source_hashes()
listing=dir(fullfile(fileparts(mfilename('fullpath')),'*.m')); names=erase({listing.name},'.m');
h=struct('file',{},'sha256',{}); root=fileparts(fileparts(mfilename('fullpath')));
for k=1:numel(names), f=fullfile(root,'function',[names{k} '.m']); if isfile(f), h(end+1)=struct('file',['function/' names{k} '.m'],'sha256',mvfp_streaming_sha256(f)); end, end %#ok<AGROW>
end
function c=strip_runtime(c), c=rmfield(c,intersect(fieldnames(c),{'outputFile','prepareOnly'})); end
function c=empty_case()
c=struct('P0Requested',[],'side','','distanceRank',[],'distance',[], ...
    'isCalibration',false,'isHoldout',false,'t',[],'distanceL1',[], ...
    'projection',[],'stableComponentL1',[],'F',[],'PFixed',[],'Praw',[], ...
    'mass',[],'rawMass',[],'rawMinDensity',[],'rawBoundaryMass',[], ...
    'initialRho',[],'initialActualP0',[],'initialM1',[],'initialF',[],'t50',[], ...
    'densitySampleTimes',[],'densitySamples',[],'rawDensitySamples',[], ...
    'events',[],'localProjectionGrowthResidual',[],'localProjectionGrowthRms',[], ...
    'qualityThresholds',[],'localStableFraction',[],'residualSeedFraction',[], ...
    'qualityPassed',false,'validationStatus','exploratory-unvalidated','outgoingModalWindow',struct());
end

function checkpoint(out,cfg)
if isempty(cfg.outputFile), error('MVFP:RevisionResidence:CheckpointRequired','Physical execution requires outputFile.'); end
folder=fileparts(cfg.outputFile); if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
save(cfg.outputFile,'out','-v7.3');
end

function a=validate_artifacts(cfg,par,x)
if ~isfile(cfg.basinFile) || ~isfile(cfg.uncertaintyFile)
    error('MVFP:RevisionResidence:MissingArtifact','Named basin and uncertainty artifacts are required.');
end
b=load(cfg.basinFile); u=load(cfg.uncertaintyFile);
a=struct('basinFile',cfg.basinFile,'basinHash',mvfp_streaming_sha256(cfg.basinFile), ...
    'uncertaintyFile',cfg.uncertaintyFile,'uncertaintyHash',mvfp_streaming_sha256(cfg.uncertaintyFile));
if ~isfield(b,'result') || ~all(isfield(b.result,{'par','config','brackets'})) || ~isequaln(b.result.par,par)
    error('MVFP:RevisionResidence:ArtifactMismatch','Basin parameters must match spectrum.');
end
bs=b.result.config.solver; cs=cfg.solver;
fields={'xMax','initStd','thresholdMode','fluxScheme','massDiscretization','initialDensityVariant','rho0','rho0X'};
for k=1:numel(fields)
    f=fields{k}; bv=[]; cv=[];
    if isfield(bs,f), bv=bs.(f); end
    if isfield(cs,f), cv=cs.(f); end
    if ~isequaln(bv,cv), error('MVFP:RevisionResidence:ArtifactMismatch','Initial protocol/solver mismatch: %s.',f); end
end
br=b.result.brackets; hit=find([br.Nx]==numel(x));
if numel(hit)~=1 || ~strcmp(br(hit).status,'bracketed') || ~isequal([br(hit).lowerP0 br(hit).upperP0],cfg.PbInterval)
    error('MVFP:RevisionResidence:ArtifactMismatch','Named basin grid and interval must match.');
end
if ~isfield(u,'uncertainty') || ~all(isfield(u.uncertainty,{'spectrumHash','basinHash','PbInterval','basinUncertainty','lambdaUncertainty'}))
    error('MVFP:RevisionResidence:ArtifactMismatch','Uncertainty artifact must bind empirical values to both input hashes.');
end
v=u.uncertainty;
if ~strcmp(v.spectrumHash,mvfp_streaming_sha256(cfg.spectrumFile)) || ~strcmp(v.basinHash,a.basinHash) || ...
        ~isequal(v.PbInterval,cfg.PbInterval) || v.basinUncertainty~=cfg.basinUncertainty || v.lambdaUncertainty~=cfg.lambdaUncertainty
    error('MVFP:RevisionResidence:ArtifactMismatch','Uncertainty provenance or values mismatch.');
end
end
