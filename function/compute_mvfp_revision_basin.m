function result = compute_mvfp_revision_basin(par,cfg)
%COMPUTE_MVFP_REVISION_BASIN Basin-B and multi-time-window analysis.
%   This is an isolated revision-v2 research driver. It records numerical
%   branch brackets and unresolved cases; it does not infer an exact basin
%   root or an infinite-time limit.

if nargin < 2, cfg=struct(); end
if isfield(cfg,'mode') && strcmp(cfg.mode,'classify')
    result=classify_tail(cfg.classifierInput,cfg.classifierOptions);
    return
end
cfg=defaults(cfg);
validate_inputs(par,cfg);
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
sources=source_hashes(cfg.projectRoot);
configHash=mvfp_value_sha256(struct('par',par,'cfg',strip_runtime(cfg)));
started=tic;
result=struct('status','exploratory-basin-brackets-not-infinite-time-proof', ...
    'createdAt',char(datetime('now')),'par',par,'config',cfg, ...
    'configHash',configHash,'sourceHashes',sources, ...
    'readoutDefinition','PFixed(t)=integral over x<=0', ...
    'basinDefinition','intersection of this two-packet initial family with numerically classified stable branches');

allCases=struct([]); allBrackets=struct([]); stopped=false;
for nx=cfg.gridN
    if toc(started) >= cfg.stageBudgetSeconds, stopped=true; break; end
    solver=cfg.solver; solver.xGridN=nx;
    refs=make_references(par,solver,cfg);
    if nnz([refs.stable]) < 2
        error('MVFP:RevisionBasin:StableBranches','Grid %d has fewer than two stable references.',nx);
    end
    points=unique([cfg.coarseP0(:); cfg.seedP0]);
    gridCases=struct([]);
    for k=1:numel(points)
        if toc(started) >= cfg.stageBudgetSeconds, stopped=true; break; end
        c=run_case(par,solver,refs,points(k),cfg.baseMaxT,cfg,sources,configHash);
        gridCases=append_struct(gridCases,c); allCases=append_struct(allCases,c);
        write_progress(cfg,allCases,allBrackets);
        if c.wallSeconds > cfg.maxCaseSeconds, stopped=true; break; end
    end
    if stopped, break; end
    bracket=locate_bracket(gridCases);
    while strcmp(bracket.status,'bracketed') && bracket.width > target_for_grid(nx,cfg)
        if toc(started) >= cfg.stageBudgetSeconds, stopped=true; break; end
        mid=mean([bracket.lowerP0 bracket.upperP0]);
        [known,idx]=case_at(gridCases,mid);
        if known, c=gridCases(idx); else
            c=run_case(par,solver,refs,mid,cfg.baseMaxT,cfg,sources,configHash);
            gridCases=append_struct(gridCases,c); allCases=append_struct(allCases,c);
        end
        if strcmp(c.classification,'unresolved')
            for extendedT=cfg.extensionT
                if extendedT <= c.maxT, continue; end
                c=run_case(par,solver,refs,mid,extendedT,cfg,sources,configHash);
                gridCases=replace_case(gridCases,c); allCases=append_struct(allCases,c);
                if ~strcmp(c.classification,'unresolved'), break; end
                if c.wallSeconds > cfg.maxCaseSeconds, break; end
            end
        end
        previous=bracket; bracket=locate_bracket(gridCases);
        if strcmp(c.classification,'unresolved') || isequaln(previous,bracket), break; end
        write_progress(cfg,allCases,allBrackets);
        if c.wallSeconds > cfg.maxCaseSeconds, stopped=true; break; end
    end
    bracket.Nx=nx; bracket.dx=solver.xMax*2/(nx-1);
    bracket.targetWidth=target_for_grid(nx,cfg);
    bracket.isTargetMet=strcmp(bracket.status,'bracketed') && bracket.width<=bracket.targetWidth ...
        && bracket.unresolvedInside==0;
    allBrackets=append_struct(allBrackets,bracket);
    result.(['grid' num2str(nx)])=struct('references',refs,'cases',gridCases,'bracket',bracket);
    write_progress(cfg,allCases,allBrackets);
    if stopped, break; end
end
result.cases=allCases; result.brackets=allBrackets;
result.stoppedAtBudget=stopped; result.wallSeconds=toc(started);
result.unresolvedSpans=unresolved_spans(allCases);
result.status=final_status(result,cfg);
save(fullfile(cfg.outputDir,'basin_B_results.mat'),'result','-v7.3');
write_outputs(cfg,result);
end

function cfg=defaults(cfg)
root=fileparts(fileparts(mfilename('fullpath')));
d=struct('projectRoot',root,'outputDir',fullfile(root,'outputs','pairwise_paper_revision_v2','second_batch','basin_B'), ...
    'gridN',[2001 4001],'coarseP0',(0:.01:.12).','seedP0',.0704, ...
    'readoutTimes',[50 100 200 400],'traceStep',5,'baseMaxT',400, ...
    'extensionT',[800 1600],'generalBracketWidth',1e-4,'fineBracketWidth',1e-3, ...
    'stageBudgetSeconds',900,'maxCaseSeconds',120,'steadyFGrid',linspace(-3,3,2401), ...
    'reuseCheckpointDirs',{{}}, ...
    'steadyRootTol',1e-11,'classifier',struct('tailWindow',20,'l1Tol',.01, ...
    'tailChangeTol',.001,'rhsL1Tol',.001,'fieldTol',.01,'tieMargin',.002), ...
    'solver',struct('xMax',1,'xGridN',2001,'initStd',.02,'thresholdMode','fixed', ...
    'odeRelTol',1e-8,'odeAbsTol',1e-12,'massDiscretization','nodal-control-volume', ...
    'fluxScheme','exponential-fit','storeRawDensity',true));
cfg=merge_struct_fields(d,cfg);
cfg.classifier=merge_struct_fields(d.classifier,cfg.classifier);
cfg.solver=merge_struct_fields(d.solver,cfg.solver);
end

function validate_inputs(par,cfg)
needed={'u','v','bias','D','sigma1','sigma2'};
assert(all(isfield(par,needed)),'Missing model parameter.');
assert(par.u==4 && par.v==20 && par.bias==-.5 && par.D==.01 && par.sigma1==2 && par.sigma2==0, ...
    'MVFP:RevisionBasin:WrongRegime','This driver is restricted to revision case B.');
assert(cfg.solver.xMax==1 && cfg.solver.initStd==.02, ...
    'MVFP:RevisionBasin:WrongProtocol','Case B requires xMax=1 and initStd=.02.');
assert(strcmp(cfg.solver.massDiscretization,'nodal-control-volume') && ...
    strcmp(cfg.solver.fluxScheme,'exponential-fit'),'MVFP:RevisionBasin:WrongSolver', ...
    'Case B requires nodal control volumes and exponential-fit flux.');
end

function refs=make_references(par,solver,cfg)
opts=struct('mvfp',solver,'steadyFGrid',cfg.steadyFGrid,'steadyRootTol',cfg.steadyRootTol);
b=compute_mvfp_steady_branches_for_param(par,opts,nan);
g=mvfp_grid_cache(solver); refs=struct([]);
for k=1:numel(b.F)
    [rho,m1,m2,~,~,readout]=mvfp_steady_density_from_F(b.F(k),par,g.x);
    rhs=mvfp_rhs(g,rho,par);
    r=struct('index',k,'F',b.F(k),'m1',m1,'m2',m2,'PFixed',readout.PssFixed, ...
        'stable',logical(b.stable(k)>0),'chi',b.chi(k),'rho',rho, ...
        'weightedRhsL1',g.w.'*abs(rhs),'selfConsistencyResidual',b.residualAtRoot(k));
    refs=append_struct(refs,r);
end
end

function c=run_case(par,solver,refs,P0,maxT,cfg,sources,configHash)
times=unique([0:cfg.traceStep:maxT cfg.readoutTimes(cfg.readoutTimes<=maxT) maxT]);
signature=struct('par',par,'solver',solver,'P0',P0,'maxT',maxT, ...
    'times',times,'classifier',cfg.classifier,'sourceHashes',sources,'configHash',configHash);
caseHash=mvfp_value_sha256(signature);
folder=fullfile(cfg.outputDir,'checkpoints'); if ~exist(folder,'dir'), mkdir(folder); end
path=fullfile(folder,[caseHash '.mat']);
if isfile(path)
    s=load(path,'c'); c=s.c;
    if ~isequaln(c.signature,signature), error('MVFP:RevisionBasin:SignatureMismatch','Checkpoint signature mismatch.'); end
    c=normalize_case_fields(c,sources);
    return
end
[compatible,old]=compatible_checkpoint([{folder} cfg.reuseCheckpointDirs(:).'],signature);
if compatible
    c=reclassify_checkpoint(old,par,solver,refs,cfg,signature,caseHash,sources);
    partial=[path '.partial.mat']; save(partial,'c','-v7.3'); movefile(partial,path,'f');
    return
end
fprintf('B basin P0=%.8f Nx=%d T=%g\n',P0,solver.xGridN,maxT);
watch=tic; trace=solve_mvfp_density_trace(par,P0,times,solver); wall=toc(watch);
g=mvfp_grid_cache(solver);
M1=g.wx.'*trace.rho; M2=g.wx2.'*trace.rho;
F=par.sigma1*M1+par.sigma2*M2.*M1+par.bias;
input=struct('x',g.x,'w',g.w,'t',trace.t,'rho',trace.rho,'F',F, ...
    'endpointRhs',mvfp_rhs(g,trace.rho(:,end),par),'references',refs, ...
    'diagnostics',struct('maxMassError',max(abs(trace.rawMass-1)), ...
    'maxBoundaryMass',max(trace.rawBoundaryMass),'minRawDensity',min(trace.rawMinDensity)));
classification=classify_tail(input,cfg.classifier);
[~,readIdx]=min(abs(trace.t(:)-cfg.readoutTimes),[],1);
compact=struct('t',trace.t,'PFixed',trace.P,'M1',M1.','M2',M2.','F',F.', ...
    'rawMass',trace.rawMass,'rawBoundaryMass',trace.rawBoundaryMass, ...
    'rawMinDensity',trace.rawMinDensity,'readoutTimes',trace.t(readIdx), ...
    'readoutP',trace.P(readIdx),'rhoInitial',trace.rho(:,1), ...
    'rhoTail',trace.rho(:,trace.t>=maxT-cfg.classifier.tailWindow),'tailTimes',trace.t(trace.t>=maxT-cfg.classifier.tailWindow));
% Full-time density is intentionally omitted; compact moments and tail checkpoints suffice.
c=struct('P0',P0,'Nx',solver.xGridN,'dx',g.dx,'maxT',maxT, ...
    'classification',classification.label,'branchIndex',classification.branchIndex, ...
    'metrics',classification,'trace',compact,'wallSeconds',wall,'caseHash',caseHash, ...
    'signature',signature,'maxMassError',max(abs(trace.rawMass-1)), ...
    'maxBoundaryMass',max(trace.rawBoundaryMass),'minRawDensity',min(trace.rawMinDensity), ...
    'originalPhysicalCaseHash','','physicalSourceHashes',sources.numerical, ...
    'classifierSourceHashes',sources.classifier,'reusedPhysicalSolve',false);
partial=[path '.partial.mat']; save(partial,'c','-v7.3'); movefile(partial,path,'f');
end

function out=classify_tail(in,opt)
stable=find([in.references.stable]); n=numel(stable);
tail=find(in.t>=in.t(end)-opt.tailWindow);
if numel(tail)<2, tail=max(1,numel(in.t)-1):numel(in.t); end
l1=nan(n,numel(tail));
for j=1:n
    ref=in.references(stable(j)).rho(:);
    l1(j,:)=in.w.'*abs(in.rho(:,tail)-ref);
end
[best,bestJ]=min(l1(:,end)); sorted=sort(l1(:,end));
if numel(sorted)>1, separation=sorted(2)-sorted(1); else, separation=inf; end
tailChange=in.w.'*abs(in.rho(:,tail(end))-in.rho(:,tail(1)));
rhsL1=in.w.'*abs(in.endpointRhs(:));
fieldDistance=abs(in.F(end)-in.references(stable(bestJ)).F);
sustained=all(l1(bestJ,:)<=opt.l1Tol);
passed=sustained && tailChange<=opt.tailChangeTol && rhsL1<=opt.rhsL1Tol && ...
    fieldDistance<=opt.fieldTol && separation>=opt.tieMargin;
diagnosticsPassed=true;
if isfield(in,'diagnostics')
    diagnosticsPassed=in.diagnostics.maxMassError<=.001 && ...
        in.diagnostics.maxBoundaryMass<=.001 && in.diagnostics.minRawDensity>=-1e-12;
end
passed=passed && diagnosticsPassed;
label='unresolved'; branchIndex=nan;
if passed
    branchIndex=stable(bestJ);
    ordered=sort([in.references(stable).F]);
    if in.references(branchIndex).F==ordered(1), label='left-stable';
    elseif in.references(branchIndex).F==ordered(end), label='right-stable';
    else, label=sprintf('stable-%d',branchIndex); end
end
out=struct('label',label,'branchIndex',branchIndex,'bestReferenceIndex',stable(bestJ), ...
    'endpointL1',best,'maxTailL1',max(l1(bestJ,:)),'tailChangeL1',tailChange, ...
    'endpointRhsWeightedL1',rhsL1,'fieldDistance',fieldDistance, ...
    'referenceSeparation',separation,'sustained',sustained, ...
    'diagnosticsPassed',diagnosticsPassed,'passed',passed);
end

function b=locate_bracket(cases)
[~,ord]=sort([cases.P0]); cases=cases(ord);
resolved=~strcmp({cases.classification},'unresolved'); labels={cases.classification};
pairs=[];
ri=find(resolved);
for k=1:numel(ri)-1
    i=ri(k); j=ri(k+1);
    if ~strcmp(labels{i},labels{j}), pairs(end+1,:)=[i j]; end %#ok<AGROW>
end
if isempty(pairs)
    b=struct('status','not-bracketed','lowerP0',nan,'upperP0',nan,'width',nan, ...
        'lowerClass','none','upperClass','none','unresolvedInside',sum(~resolved), ...
        'transitionIntervals',struct([])); return
end
intervals=struct([]);
for k=1:size(pairs,1)
    ij=pairs(k,:); intervals=append_struct(intervals,struct( ...
        'lowerP0',cases(ij(1)).P0,'upperP0',cases(ij(2)).P0, ...
        'lowerClass',cases(ij(1)).classification,'upperClass',cases(ij(2)).classification, ...
        'unresolvedInside',sum(~resolved(ij(1):ij(2)))));
end
if size(pairs,1)>1
    b=struct('status','ambiguous-multiple-transitions','lowerP0',nan,'upperP0',nan, ...
        'width',nan,'lowerClass','multiple','upperClass','multiple', ...
        'unresolvedInside',sum(~resolved),'transitionIntervals',intervals); return
end
ij=pairs(1,:);
b=struct('status','bracketed','lowerP0',cases(ij(1)).P0,'upperP0',cases(ij(2)).P0, ...
    'width',cases(ij(2)).P0-cases(ij(1)).P0,'lowerClass',cases(ij(1)).classification, ...
    'upperClass',cases(ij(2)).classification,'unresolvedInside', ...
    sum(~resolved(ij(1):ij(2))),'transitionIntervals',intervals);
end

function spans=unresolved_spans(cases)
spans=struct([]);
for nx=unique([cases.Nx])
    q=cases([cases.Nx]==nx); [~,o]=sort([q.P0]); q=q(o); mask=strcmp({q.classification},'unresolved');
    starts=find(mask & [true ~mask(1:end-1)]); stops=find(mask & [~mask(2:end) true]);
    for k=1:numel(starts)
        spans=append_struct(spans,struct('Nx',nx,'lowerP0',q(starts(k)).P0, ...
            'upperP0',q(stops(k)).P0,'nCases',stops(k)-starts(k)+1));
    end
end
end

function write_progress(cfg,cases,brackets)
s=struct('cases',cases,'brackets',brackets,'updatedAt',char(datetime('now')));
save(fullfile(cfg.outputDir,'basin_B_progress.mat'),'s','-v7.3');
if ~isempty(cases), writetable(case_table(cases),fullfile(cfg.outputDir,'basin_B_progress.csv')); end
end

function write_outputs(cfg,result)
writetable(case_table(result.cases),fullfile(cfg.outputDir,'basin_B_summary.csv'));
if ~isempty(result.brackets), writetable(bracket_table(result.brackets),fullfile(cfg.outputDir,'basin_B_brackets.csv')); end
fid=fopen(fullfile(cfg.outputDir,'basin_B_report.md'),'w','n','UTF-8'); clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Revision v2 — case B basin and multiple time windows\n\n');
fprintf(fid,'Status: `%s`. Total wall time: %.1f s.\n\n',result.status,result.wallSeconds);
fprintf(fid,'Classification requires sustained weighted L1 proximity to a grid-specific stable Boltzmann density, tail density change, endpoint RHS, and feedback-field agreement. `P(T)` is reported separately and is not used alone for branch assignment.\n\n');
for k=1:numel(result.brackets)
    b=result.brackets(k); fprintf(fid,'- Nx=%d: %s; Pb interval [%.8g, %.8g], width %.3g; unresolved samples inside=%d; target met=%d.\n', ...
        b.Nx,b.status,b.lowerP0,b.upperP0,b.width,b.unresolvedInside,b.isTargetMet);
end
fprintf(fid,'\nThese finite-window classifications do not prove an infinite-time basin boundary. Unresolved states are retained in the MAT/CSV outputs. Source and configuration SHA-256 values are stored in the MAT file and each checkpoint signature.\n');
end

function t=bracket_table(brackets)
rows=struct([]);
for k=1:numel(brackets)
    b=brackets(k); rows=append_struct(rows,struct('Nx',b.Nx,'dx',b.dx, ...
        'status',string(b.status),'lowerP0',b.lowerP0,'upperP0',b.upperP0, ...
        'width',b.width,'lowerClass',string(b.lowerClass),'upperClass',string(b.upperClass), ...
        'unresolvedInside',b.unresolvedInside,'targetWidth',b.targetWidth,'isTargetMet',b.isTargetMet, ...
        'nTransitionIntervals',numel(b.transitionIntervals)));
end
t=struct2table(rows);
end

function [found,c]=compatible_checkpoint(folders,sig)
found=false; c=[];
for d=1:numel(folders)
folder=folders{d}; if ~isfolder(folder), continue; end
files=dir(fullfile(folder,'*.mat'));
for k=1:numel(files)
    try
        q=load(fullfile(files(k).folder,files(k).name),'c');
        if ~isfield(q,'c') || ~isfield(q.c,'signature'), continue; end
        old=q.c.signature;
        fields={'par','solver','P0','maxT','times'};
        same=true;
        for j=1:numel(fields), same=same && isfield(old,fields{j}) && isequaln(old.(fields{j}),sig.(fields{j})); end
        if same, found=true; c=q.c; return; end
    catch
    end
end
end
end

function c=reclassify_checkpoint(c,par,solver,refs,cfg,signature,caseHash,currentSources)
g=mvfp_grid_cache(solver); tailRho=c.trace.rhoTail; tailT=c.trace.tailTimes(:).';
tailF=c.trace.F(c.trace.t>=c.maxT-cfg.classifier.tailWindow).';
input=struct('x',g.x,'w',g.w,'t',tailT,'rho',tailRho,'F',tailF, ...
    'endpointRhs',mvfp_rhs(g,tailRho(:,end),par),'references',refs, ...
    'diagnostics',struct('maxMassError',c.maxMassError, ...
    'maxBoundaryMass',c.maxBoundaryMass,'minRawDensity',c.minRawDensity));
classification=classify_tail(input,cfg.classifier);
oldHash=c.caseHash;
if isfield(c,'physicalSourceHashes'), oldSources=c.physicalSourceHashes;
elseif isfield(c.signature,'sourceHashes') && isstruct(c.signature.sourceHashes) && isfield(c.signature.sourceHashes,'numerical')
    oldSources=c.signature.sourceHashes.numerical;
else, oldSources=c.signature.sourceHashes; end
c.classification=classification.label; c.branchIndex=classification.branchIndex;
c.metrics=classification; c.signature=signature; c.caseHash=caseHash;
c.originalPhysicalCaseHash=oldHash; c.physicalSourceHashes=oldSources;
c.classifierSourceHashes=currentSources.classifier; c.reusedPhysicalSolve=true;
c=normalize_case_fields(c,currentSources);
end

function c=normalize_case_fields(c,currentSources)
if ~isfield(c,'originalPhysicalCaseHash'), c.originalPhysicalCaseHash=''; end
if ~isfield(c,'physicalSourceHashes'), c.physicalSourceHashes=currentSources.numerical; end
if ~isfield(c,'classifierSourceHashes'), c.classifierSourceHashes=currentSources.classifier; end
if ~isfield(c,'reusedPhysicalSolve'), c.reusedPhysicalSolve=false; end
end

function t=case_table(cases)
rows=struct([]);
for k=1:numel(cases)
    c=cases(k); p=nan(1,4); p(1:numel(c.trace.readoutP))=c.trace.readoutP;
    r=struct('Nx',c.Nx,'dx',c.dx,'P0',c.P0,'maxT',c.maxT,'classification',string(c.classification), ...
        'branchIndex',c.branchIndex,'P50',p(1),'P100',p(2),'P200',p(3),'P400',p(4), ...
        'endpointL1',c.metrics.endpointL1,'maxTailL1',c.metrics.maxTailL1, ...
        'tailChangeL1',c.metrics.tailChangeL1,'endpointRhsWeightedL1',c.metrics.endpointRhsWeightedL1, ...
        'fieldDistance',c.metrics.fieldDistance,'maxMassError',c.maxMassError, ...
        'maxBoundaryMass',c.maxBoundaryMass,'minRawDensity',c.minRawDensity, ...
        'wallSeconds',c.wallSeconds,'caseHash',string(c.caseHash));
    rows=append_struct(rows,r);
end
t=struct2table(rows);
end

function s=source_hashes(root)
s=struct();
s.numerical=hash_closure(root,{'solve_mvfp_density_trace'});
s.classifier=hash_closure(root,{'compute_mvfp_revision_basin', ...
    'compute_mvfp_steady_branches_for_param','mvfp_steady_density_from_F','mvfp_rhs'});
end

function hashes=hash_closure(root,seeds)
folder=fullfile(root,'function'); files=dir(fullfile(folder,'*.m'));
available=erase({files.name},'.m'); queue=seeds(:).'; seen={};
while ~isempty(queue)
    name=queue{1}; queue(1)=[];
    if ismember(name,seen) || ~ismember(name,available), continue; end
    seen{end+1}=name; %#ok<AGROW>
    text=fileread(fullfile(folder,[name '.m']));
    calls=regexp(text,'(?<![.''"])\<([A-Za-z]\w*)\>\s*\(','tokens');
    if ~isempty(calls), queue=[queue cellfun(@(z)z{1},calls,'UniformOutput',false)]; end %#ok<AGROW>
end
seen=sort(seen); hashes=struct('file',{},'sha256',{});
for k=1:numel(seen)
    rel=['function/' seen{k} '.m']; hashes(k)=struct('file',rel, ...
        'sha256',mvfp_streaming_sha256(fullfile(root,rel))); %#ok<AGROW>
end
end
function c=strip_runtime(c), c=rmfield(c,intersect(fieldnames(c),{'outputDir','stageBudgetSeconds','maxCaseSeconds'})); end
function x=append_struct(x,v), if isempty(x), x=v; else, x(end+1)=v; end, end
function x=replace_case(x,v)
idx=find([x.Nx]==v.Nx & abs([x.P0]-v.P0)<1e-12,1,'last');
if isempty(idx), x=append_struct(x,v); else, x(idx)=v; end
end
function [yes,idx]=case_at(c,p), idx=find(abs([c.P0]-p)<1e-12,1); yes=~isempty(idx); end
function v=target_for_grid(nx,cfg), if nx==max(cfg.gridN), v=cfg.fineBracketWidth; else, v=cfg.generalBracketWidth; end, end
function s=final_status(r,cfg)
if r.stoppedAtBudget, s='budget-stopped-pending';
elseif isempty(r.brackets) || any(~[r.brackets.isTargetMet]), s='numerically-unresolved-or-bracket-pending';
else, s='finite-window-branch-brackets-computed-general-comparison-pending-review'; end
end
