function result = run_mvfp_paper_revision_stage(cfg,stage)
% Explicit stages only; research gates are never auto-passed.
allowed={'preflight','controls','convergence_pilot'};
if ~ismember(stage,allowed)
    error('MVFP:Revision:StageNotImplemented', ...
        'Implemented stages: %s. Requested stage has not been implemented.',strjoin(allowed,', '));
end
if cfg.parA.sigma2~=0 || cfg.parB.sigma2~=0
    error('MVFP:Revision:PairwiseRequired','Revision v2 first batch requires sigma2=0.');
end
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
provenance=source_provenance(cfg.projectRoot);
result=struct('stage',stage,'status','exploratory-not-convergence-certified', ...
    'config',cfg,'sourceHashes',provenance,'createdAt',char(datetime('now')), ...
    'readoutDefinition','PFixed=integral[-L,0] rho dx; P0 is packet weight');
switch stage
    case 'preflight'
        [result,summary]=preflight(cfg,result);
    case 'controls'
        [result,summary]=controls(cfg,result);
    case 'convergence_pilot'
        [result,summary]=convergence_pilot(cfg,result);
end
save(fullfile(cfg.outputDir,[stage '.mat']),'result','-v7.3');
writetable(summary,fullfile(cfg.outputDir,[stage '_summary.csv']));
write_report(cfg,result,summary);
fprintf('Saved %s: %s\n',stage,cfg.outputDir);
end

function sources=source_provenance(root)
files=dir(fullfile(root,'function','*.m'));
sources=struct('file',{},'sha256',{});
for k=1:numel(files)
    sources(k)=struct('file',['function/' files(k).name], ...
        'sha256',mvfp_streaming_sha256(fullfile(files(k).folder,files(k).name)));
end
sources(end+1)=struct('file','run_recomputation.m', ...
    'sha256',mvfp_streaming_sha256(fullfile(root,'run_recomputation.m')));
end

function [result,summary]=preflight(cfg,result)
opts=struct('mvfp',cfg.mvfp,'steadyFGrid',cfg.steadyFGrid,'steadyRootTol',cfg.steadyRootTol);
rows=struct([]);
for label={'A','B'}
    name=label{1}; par=cfg.(['par' name]);
    branches=compute_mvfp_steady_branches_for_param(par,opts,nan);
    result.(['branches' name])=branches;
    for k=1:numel(branches.F)
        chi=par.sigma1*(branches.m2(k)-branches.m1(k)^2)/par.D;
        row=struct('regime',string(name),'F',branches.F(k),'M1',branches.m1(k), ...
            'chi',chi,'legacyStable',branches.stable(k),'residual',branches.residualAtRoot(k));
        if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
    end
end
summary=struct2table(rows);
gridCfg=cfg.mvfp; gridCfg.xGridN=101;
g=mvfp_grid_cache(gridCfg);
rho=normalize_density(g.x,exp(-g.x.^2).*(1+.2*g.x));
new=mvfp_rhs(g,rho,cfg.parB);
old=mvfp_rhs(g.x,rho,cfg.parB);
fixed=mvfp_rhs_fixed_feedback(g,rho,cfg.parB,.04);
result.massProbe=struct('newTrapzDerivative',g.w.'*new, ...
    'fixedTrapzDerivative',g.w.'*fixed,'oldRectangleDerivative',g.dx*sum(old), ...
    'oldTrapzDerivative',g.w.'*old,'scope','synthetic operator test, not a production error estimate');
assert(abs(result.massProbe.newTrapzDerivative)<1e-10);
assert(abs(result.massProbe.fixedTrapzDerivative)<1e-10);
result.status='protocol-preflight-complete-dynamics-unvalidated';
end

function [result,summary]=controls(cfg,result)
rows=struct([]); cases=cell(1,numel(cfg.controlP0));
for k=1:numel(cfg.controlP0)
    P0=cfg.controlP0(k);
    opts=struct('mvfp',cfg.mvfp,'T',cfg.T,'feedbackControlP0List',P0, ...
        'feedbackControlTimes',unique([0:cfg.outputTimeStep:cfg.T,cfg.T]), ...
        'feedbackControlsCommonRho0',true,'readoutThreshold',.5, ...
        'dynamicFeedbackReadoutLevels',.2:.2:.8);
    signature=struct('par',cfg.parA,'opts',opts,'sourceHashes',result.sourceHashes);
    [entry,path,found]=read_case(cfg,'controls',signature);
    if ~found
        fprintf('Control P0=%.7g, Nx=%d\n',P0,cfg.mvfp.xGridN);
        watch=tic;
        diagnostic=compute_mvfp_feedback_control_diagnostics(cfg.parA,opts);
        entry=struct('diagnostic',diagnostic,'wallSeconds',toc(watch),'signature',signature);
    end
    [entry.diagnostic,eventFieldsChanged]=apply_revision_control_events(entry.diagnostic);
    if ~found || eventFieldsChanged
        save_case(path,entry);
    end
    c=entry.diagnostic.cases(1);
    assert(c.initialL1.fullVsFrozen<1e-12 && c.initialL1.fullVsNoCoupling<1e-12);
    assert(c.initialDriftL1FullVsFrozen<1e-12);
    cases{k}=entry;
    for control={'full','frozenFeedback','noCoupling'}
        s=c.(control{1}); event=mvfp_revision_first_arrival(s.t,s.PFixed,.5);
        row=struct('P0',P0,'control',string(control{1}),'F0',s.F(1), ...
            'oldGuessForCoupled',initial_field_guess(P0,cfg.parA), ...
            'PFixed0',s.PFixed(1),'PFixedT',s.PFixed(end),'t50',event.time, ...
            'eventStatus',string(event.status),'maxMassError',max(s.massError), ...
            'maxBoundaryMass',max(s.boundaryMass),'rawMinDensity',min(s.rawMinDensity), ...
            'initialL1FullFrozen',c.initialL1.fullVsFrozen, ...
            'initialL1FullUncoupled',c.initialL1.fullVsNoCoupling, ...
            'wallSecondsForThree',entry.wallSeconds);
        if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
    end
end
result.cases=cases;
summary=struct2table(rows);
result.status='controls-computed-spatial-and-event-convergence-pending';
plot_controls(cfg,result);
end

function [result,summary]=convergence_pilot(cfg,result)
rows=struct([]); cases={}; clockStart=tic; lastSeconds=0; stopped=false;
for nx=cfg.pilotGridN
    for label={'A','B'}
        name=label{1}; par=cfg.(['par' name]); points=cfg.(['pilotP0' name]);
        for P0=points
            if toc(clockStart)>cfg.pilotMaxStageSeconds || lastSeconds>cfg.pilotMaxLastCaseSeconds
                stopped=true; break;
            end
            solver=cfg.mvfp; solver.xGridN=nx;
            times=unique([0:cfg.outputTimeStep:cfg.T,cfg.T]);
            signature=struct('par',par,'solver',solver,'P0',P0,'times',times, ...
                'sourceHashes',result.sourceHashes);
            [entry,path,found]=read_case(cfg,'convergence_pilot',signature);
            if ~found
                fprintf('Timing %s P0=%.7g Nx=%d T=%g\n',name,P0,nx,cfg.T);
                watch=tic; trace=solve_mvfp_density_trace(par,P0,times,solver);
                entry=struct('trace',trace,'wallSeconds',toc(watch),'signature',signature);
                save_case(path,entry);
            end
            lastSeconds=entry.wallSeconds;
            trace=entry.trace; event=mvfp_revision_first_arrival(trace.t,trace.P,.5);
            row=struct('regime',string(name),'P0',P0,'Nx',nx,'dx',2*solver.xMax/(nx-1), ...
                'T',cfg.T,'PFixedT',trace.P(end),'t50',event.time,'eventStatus',string(event.status), ...
                'maxMassError',max(abs(trace.rawMass-1)), ...
                'maxBoundaryMass',max(trace.rawBoundaryMass),'rawMinDensity',min(trace.rawMinDensity), ...
                'wallSeconds',entry.wallSeconds);
            if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
            cases{end+1}=entry; %#ok<AGROW>
            writetable(struct2table(rows),fullfile(cfg.outputDir,'convergence_pilot_progress.csv'));
        end
        if stopped, break; end
    end
    if stopped, break; end
end
result.cases=cases; result.stoppedAtBudget=stopped;
result.status='representative-timing-only-no-production-grid-certified';
summary=struct2table(rows);
end

function [diagnostic,changed]=apply_revision_control_events(diagnostic)
changed=false;
for k=1:numel(diagnostic.cases)
    for control={'full','frozenFeedback','noCoupling'}
        name=control{1}; s=diagnostic.cases(k).(name);
        levels=s.readoutLevels(:).';
        events=repmat(mvfp_revision_first_arrival(s.t,s.PFixed,levels(1)),1,numel(levels));
        for j=1:numel(levels)
            events(j)=mvfp_revision_first_arrival(s.t,s.PFixed,levels(j));
        end
        transitionEvent=mvfp_revision_first_arrival(s.t,s.PFixed,.5);
        arrivalTimes=[events.time];
        needsUpdate=~isfield(s,'firstArrivalEvents') || ...
            ~isequaln(s.firstArrivalEvents,events) || ...
            ~isfield(s,'transitionEvent') || ~isequaln(s.transitionEvent,transitionEvent) || ...
            ~isequaln(s.arrivalTimes,arrivalTimes) || ~isequaln(s.transitionTime,transitionEvent.time);
        if needsUpdate
            s.firstArrivalEvents=events;
            s.transitionEvent=transitionEvent;
            s.arrivalTimes=arrivalTimes;
            s.transitionTime=transitionEvent.time;
            diagnostic.cases(k).(name)=s;
            changed=true;
        end
    end
end
end

function [entry,path,found]=read_case(cfg,stage,signature)
folder=fullfile(cfg.outputDir,'cases',stage);
if ~exist(folder,'dir'), mkdir(folder); end
path=fullfile(folder,[mvfp_value_sha256(signature) '.mat']);
found=isfile(path); entry=[];
if found
    loaded=load(path,'entry'); entry=loaded.entry;
    assert(isequaln(entry.signature,signature),'Case provenance differs.');
end
end
function save_case(path,entry)
partial=[path '.partial.mat'];
save(partial,'entry','-v7.3');
movefile(partial,path,'f');
end

function write_report(cfg,result,summary)
path=fullfile(cfg.outputDir,[result.stage '_report.md']);
fid=fopen(path,'w','n','UTF-8'); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Revision v2: %s\n\nStatus: `%s`.\n\n',result.stage,result.status);
fprintf(fid,'Generated: %s. Summary rows: %d.\n\n',result.createdAt,height(summary));
fprintf(fid,'Source: `%s.mat`; table: `%s_summary.csv`.\n\n',result.stage,result.stage);
fprintf(fid,'Solver: %s, Nx=%d, RelTol=%g, AbsTol=%g.\n\n', ...
    cfg.mvfp.massDiscretization,cfg.mvfp.xGridN,cfg.mvfp.odeRelTol,cfg.mvfp.odeAbsTol);
fprintf(fid,'These first-batch results do not certify spatial convergence, a new basin boundary, a spectral growth rate, or a publication claim.\n');
if strcmp(result.stage,'convergence_pilot')
    fprintf(fid,'\nBudget stop: %d. Actual grids are listed in the CSV. Threshold bracketing, domain and tolerance checks are still required.\n',result.stoppedAtBudget);
end
end

function plot_controls(cfg,result)
fig=figure('Visible','off','Color','w','Position',[50,50,1100,740]);
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
n=numel(result.cases);
cols=ceil(sqrt(n)); rows=ceil(n/cols);
tiledlayout(rows,cols,'TileSpacing','compact');
for k=1:numel(result.cases)
    nexttile; c=result.cases{k}.diagnostic.cases(1);
    plot(c.full.t,c.full.PFixed,'LineWidth',1.5); hold on;
    plot(c.frozenFeedback.t,c.frozenFeedback.PFixed,'--','LineWidth',1.5);
    plot(c.noCoupling.t,c.noCoupling.PFixed,':','LineWidth',1.5);
    ylim([0,1]); xlabel('t'); ylabel('P_{fixed}(t)');
    title(sprintf('Same initial density: P_0=%.4g',c.P0));
    legend('full','frozen F(0)','uncoupled','Location','best');
end
sgtitle(sprintf('Exploratory controls: D=.005, sigma_1=.5, Nx=%d; convergence pending',cfg.mvfp.xGridN));
exportgraphics(fig,fullfile(cfg.outputDir,'controls_first_batch.png'),'Resolution',180);
exportgraphics(fig,fullfile(cfg.outputDir,'controls_first_batch.pdf'),'ContentType','vector');
end
