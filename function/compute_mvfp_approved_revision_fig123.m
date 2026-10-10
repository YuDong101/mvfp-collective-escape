function result=compute_mvfp_approved_revision_fig123(cfg)
%COMPUTE_MVFP_APPROVED_REVISION_FIG123 Isolated approved D1/D2/D3 data work.
% Each figure is independently checkpointed and closed figures are preserved.
out=cfg.outputDir; if ~isfolder(out), mkdir(out); end
root=cfg.projectRoot;
base=mvfp_paper_revision_config(root); par=base.parA;
solver=base.mvfp; solver.xGridN=2001; solver.fluxScheme='exponential-fit';
sourceRoot=fullfile(root,'outputs','pairwise_paper_revision_v2');
coarsePath=fullfile(sourceRoot,'fitted_flux_controls','controls.mat');
finePath=fullfile(sourceRoot,'fourth_batch','figures123','fine_controls_nx4001','controls.mat');
twPath=fullfile(sourceRoot,'second_batch','time_windows_A','nx2001','time_windows.mat');
twFinePath=fullfile(sourceRoot,'second_batch','time_windows_A','nx4001','time_windows.mat');
if cfg.smoke
    solver.xGridN=201;
    tr=solve_mvfp_density_trace(par,.08,0:.1:1,solver);
    s=pack_trace(tr,par);
    assert(max(abs(s.rawMass-1))<1e-7 && all(isfinite(s.PFixed)));
    g=mvfp_grid_cache(solver); rho0=initial_density(g.x,par,.08,solver.initStd,solver);
    assert(abs(trapz(g.x,rho0)-1)<1e-12);
    cf=solver; cf.rho0=rho0; cf.rho0X=g.x;
    frozen=solve_mvfp_density_trace_fixed_feedback(par,.08,0:.1:1,cf,s.F(1));
    assert(max(abs(frozen.rhoRaw(:,1)-tr.rhoRaw(:,1)))<1e-12);
    ev=mvfp_revision_first_arrival([0 1],[.1 .2],.5);
    assert(~isfinite(ev.time));
    result=struct('smokePassed',true,'maxRawMassError',max(abs(s.rawMass-1)), ...
        'initialDriftDifference',max(abs(mvfp_rhs(g,rho0,par)-mvfp_rhs_fixed_feedback(g,rho0,par,s.F(1)))));
    assert(result.initialDriftDifference<1e-9);
    write_json(fullfile(out,'smoke_validation.json'),result); return
end
cs=load(coarsePath,'result'); coarse=cs.result;
fs=load(finePath,'result'); fine=fs.result;
ts=load(twPath,'result'); tw=ts.result;
tfs=load(twFinePath,'result'); twFine=tfs.result;
assert(isequaln(coarse.config.parA,par) && isequaln(tw.config.par,par));
assert(isequaln(coarse.config.mvfp,solver),'Coarse cache solver does not match.');
assert(tw.config.mvfp.xGridN==2001 && strcmp(tw.config.mvfp.fluxScheme,'exponential-fit'));
verify_sources(root,tw.sourceHashes);
verify_sources(root,twFine.sourceHashes);
verify_control_sources(root,coarse.sourceHashes);
verify_control_sources(root,fine.sourceHashes);
sources={coarsePath,finePath,twPath,twFinePath};
result=struct('createdAt',char(datetime('now')),'par',par,'config',cfg, ...
    'sources',{sources},'status','candidate data; specific plots require user approval');
for k=1:numel(sources), result.sourceSHA256{k}=mvfp_streaming_sha256(sources{k}); end
thresholds=[struct2table(tw.thresholds);struct2table(twFine.thresholds)];
thr=thresholds(thresholds.T==100,:);
result.thresholdEnvelopeT100=[min(thr.P0Lower),max(thr.P0Upper)];
writetable(thresholds,fullfile(out,'fig3_threshold_brackets.csv'));
writetable(struct2table(tw.windowStatus),fullfile(out,'fig3_window_status.csv'));
if cfg.runFigure1
    fprintf('Computing Figure 1 actual-density frozen diagnostics.\n');
    P0=unique([linspace(0,1,501),.06,.0647,.07,.08,.09]);
    g=mvfp_grid_cache(solver); n=numel(P0); vals=nan(n,11);
    for k=1:n
        rho=initial_density(g.x,par,P0(k),solver.initStd,solver);
        F=par.sigma1*trapz(g.x,g.x.*rho)+par.bias;
        ex=effective_potential_extrema_from_F(F,par);
        [~,~,lLR,lRL,~,~,~,~,~,~,~,~,~]=compute_mvfp_mfpt_point_from_F(F,par,g.x,-1,1,20001);
        bLR=nan; bRL=nan;
        if isfinite(ex.xs)
            u=effective_potential_from_F([ex.xL ex.xs ex.xR],par,F);
            bLR=u(2)-u(1); bRL=u(2)-u(3);
        end
        vals(k,:)=[P0(k),F,ex.xL,ex.xs,ex.xR,ex.xSingle,bLR,bRL,lLR/log(10),lRL/log(10),trapz(g.x,rho)];
    end
    names={'P0','F0','xLeftMinimum','xSaddle','xRightMinimum','xSingleMinimum', ...
        'barrierLeftToRight','barrierRightToLeft','log10TauLeftToSaddle','log10TauRightToSaddle','initialMass'};
    writetable(array2table(vals,'VariableNames',names),fullfile(out,'fig1_frozen_summary.csv'));
    x=linspace(-1,1,2001).'; p=[0 .25 .5 .75 1]; U=zeros(numel(x),5); F0=zeros(1,5);
    for k=1:5
        rho=initial_density(g.x,par,p(k),solver.initStd,solver);
        F0(k)=par.sigma1*trapz(g.x,g.x.*rho)+par.bias;
        U(:,k)=effective_potential_from_F(x,par,F0(k)); U(:,k)=U(:,k)-min(U(:,k));
    end
    write_json(fullfile(out,'fig1_potentials.json'),struct('x',x,'P0',p,'F0',F0,'UminusMin',U));
    f1=struct('mfptDefinition','Boltzmann basin averaged MFPT, reflecting outer edge at -1 or +1 and absorbing saddle; direct log quadrature', ...
        'NxInitialDensity',2001,'mfptDomainGridN',20001,'dynamicThresholdEnvelopeT100',result.thresholdEnvelopeT100, ...
        'actualInitialField',true,'oldGuessUsedOnlyToPreparePacketCenters',true);
    write_json(fullfile(out,'fig1_metadata.json'),f1);
end
if cfg.runFigure2
    wanted=[0 .04 .06 .0647 .07 .08 .10 .16]; cases=cell(size(wanted)); missing=[];
    for k=1:numel(wanted)
        c=find_control(coarse,wanted(k));
        if ~isempty(c), cases{k}=from_control(c.full); else, missing(end+1)=k; end %#ok<AGROW>
    end
    fprintf('Figure 2: reusing %d exact coarse controls; %d new full traces.\n',numel(wanted)-numel(missing),numel(missing));
    jobs=cell(size(missing));
    for j=1:numel(missing), jobs{j}=make_job(wanted(missing(j)),100,.1,solver,par,out,'fig2'); end
    computed=run_jobs(jobs,cfg);
    for j=1:numel(missing), cases{missing(j)}=computed{j}; end
    tables=cell(1,numel(cases)); events=struct([]); diagnostic=struct([]);
    for k=1:numel(cases)
        s=cases{k}; tables{k}=trace_table(wanted(k),2001,'full',s);
        e=event_rows(wanted(k),s); events=[events,e]; %#ok<AGROW>
        diagnostic=[diagnostic,diag_row(wanted(k),'full',s)]; %#ok<AGROW>
    end
    writetable(vertcat(tables{:}),fullfile(out,'fig2_dynamic_full_traces.csv'));
    writetable(struct2table(events),fullfile(out,'fig2_event_timing.csv'));
    controls={}; validation=struct([]);
    for p=[.06 .08]
        c=find_control(coarse,p); f=find_control(fine,p);
        assert(c.initialL1.fullVsFrozen<1e-12 && c.initialL1.fullVsNoCoupling<1e-12);
        for key={'full','frozenFeedback','noCoupling'}
            a=from_control(c.(key{1})); b=from_control(f.(key{1}));
            controls{end+1}=trace_table(p,2001,key{1},a); %#ok<AGROW>
            ea=mvfp_revision_first_arrival(a.t,a.PFixed,.5); eb=mvfp_revision_first_arrival(b.t,b.PFixed,.5);
            row=struct('P0',p,'control',string(key{1}),'PFixedT_Nx2001',a.PFixed(end), ...
                'PFixedT_Nx4001',b.PFixed(end),'absEndpointDifference',abs(a.PFixed(end)-b.PFixed(end)), ...
                't50_Nx2001',ea.time,'t50_Nx4001',eb.time, ...
                'maxMassErrorNx4001',max(abs(b.rawMass-1)),'maxBoundaryMassNx4001',max(b.boundaryMass));
            validation=[validation,row]; %#ok<AGROW>
        end
    end
    writetable(vertcat(controls{:}),fullfile(out,'fig2_same_density_controls.csv'));
    writetable(struct2table(validation),fullfile(out,'fig2_fine_grid_validation.csv'));
    writetable(struct2table(diagnostic),fullfile(out,'fig2_numerical_diagnostics.csv'));
    s=cases{find(abs(wanted-.08)<1e-12,1)}; x=linspace(-1,1,801).';
    U=effective_potential_from_F(x,par,s.F); U=U-min(U,[],1);
    write_json(fullfile(out,'fig2_dynamic_potential.json'),struct('P0',.08,'Nx',2001,'x',x,'t',s.t, ...
        'UminusMin',U,'xLeftMinimum',s.xL,'xSaddle',s.xs,'xRightMinimum',s.xR,'xSingleMinimum',s.xSingle));
end
if cfg.runFigure3
    tight=tw.config.mvfp; pExisting=[tw.entries.P0];
    % Merge numerically identical preparation weights (binary .01 grids can
    % differ by a few ulps); 12 decimals preserve all search-bracket digits.
    wanted=unique(round([0:.01:.3,pExisting],12)); times=0:.5:100;
    P=nan(numel(wanted),numel(times)); t50=nan(size(wanted)); diagnostics=struct([]); source=strings(size(wanted));
    missing=[];
    for k=1:numel(wanted)
        idx=find(abs(pExisting-wanted(k))<1e-12,1);
        if isempty(idx), missing(end+1)=k; continue; end %#ok<AGROW>
        e=tw.entries(idx); P(k,:)=interp1(e.t,e.PFixed,times);
        ev=mvfp_revision_first_arrival(times,P(k,:),.5); t50(k)=ev.time;
        diagnostics=[diagnostics,struct('P0',wanted(k),'maxMassError',e.maxMassError, ...
            'maxBoundaryMass',e.maxBoundaryMass,'minRawDensity',e.minRawDensity)]; %#ok<AGROW>
        source(k)="existing Nx2001 time-window cache";
    end
    fprintf('Figure 3: reusing %d readout traces; %d new traces; density .07 to200.\n',numel(wanted)-numel(missing),numel(missing));
    jobs=cell(size(missing));
    for j=1:numel(missing), jobs{j}=make_job(wanted(missing(j)),100,.1,tight,par,out,'fig3'); end
    computed=run_jobs(jobs,cfg);
    for j=1:numel(missing)
        k=missing(j); s=computed{j}; P(k,:)=interp1(s.t,s.PFixed,times);
        ev=mvfp_revision_first_arrival(s.t,s.PFixed,.5); t50(k)=ev.time;
        diagnostics=[diagnostics,struct('P0',wanted(k),'maxMassError',max(abs(s.rawMass-1)), ...
            'maxBoundaryMass',max(s.boundaryMass),'minRawDensity',min(s.rawMinDensity))]; %#ok<AGROW>
        source(k)="new Nx2001 tight-tolerance .1-time sample cache";
    end
    write_json(fullfile(out,'fig3_readout_heatmap.json'),struct('P0',wanted,'t',times,'PFixed',P, ...
        'Nx',2001,'solver',tight,'source',source));
    status=repmat("reached",size(wanted)); status(~isfinite(t50))="right-censored";
    plotted=t50; plotted(~isfinite(plotted))=100;
    writetable(table(wanted(:),t50(:),plotted(:),status(:),source(:), ...
        'VariableNames',{'P0','t50Actual','plottedTime','statusByT100','source'}),fullfile(out,'fig3_t50_censoring_T100.csv'));
    reps={}; initials={}; g=mvfp_grid_cache(tight);
    for p=[.06 .07 .08]
        e=tw.entries(find(abs(pExisting-p)<1e-12,1)); ix=e.t<=200;
        reps{end+1}=table(repmat(p,sum(ix),1),e.t(ix).',e.PFixed(ix).',e.M1(ix).', ...
            'VariableNames',{'P0','t','PFixed','M1'}); %#ok<AGROW>
        rho0=initial_density(g.x,par,p,tight.initStd,tight);
        initials{end+1}=table(repmat(p,numel(g.x),1),g.x,rho0,'VariableNames',{'P0','x','rho0'}); %#ok<AGROW>
    end
    writetable(vertcat(reps{:}),fullfile(out,'fig3_representative_readout_traces.csv'));
    writetable(vertcat(initials{:}),fullfile(out,'fig3_initial_densities.csv'));
    job=make_job(.07,200,.1,tight,par,out,'fig3_density'); job.keepDensity=true;
    s=run_one(job);
    ii=1:2:numel(s.x); jj=1:5:numel(s.t);
    write_json(fullfile(out,'fig3_density_evolution.json'),struct('P0',.07,'Nx',2001, ...
        'x',s.x(ii),'t',s.t(jj),'rhoRaw',s.rhoRaw(ii,jj),'rawMass',s.rawMass(jj), ...
        'boundaryMass',s.boundaryMass(jj),'displaySampling','every second spatial node and every fifth .1 time sample; no density renormalization'));
    e=tw.entries(find(abs(pExisting-.07)<1e-12,1)); idx=e.t<=200;
    delta=max(abs(interp1(s.t,s.PFixed,e.t(idx))-e.PFixed(idx)));
    assert(delta<1e-4,'New density trace and existing representative trace disagree.');
    write_json(fullfile(out,'fig3_density_validation.json'),struct('maxReadoutDifferenceVsTimeWindowCache',delta, ...
        'maxRawMassError',max(abs(s.rawMass-1)),'maxBoundaryMass',max(s.boundaryMass)));
    writetable(struct2table(diagnostics),fullfile(out,'fig3_numerical_diagnostics.csv'));
end
result.solverFigure2=solver; result.solverFigure3=tw.config.mvfp;
write_json(fullfile(out,'manifest.json'),result);
fprintf('Approved Figure 1--3 data exports complete: %s\n',out);
end

function c=find_control(r,p)
c=[];
for k=1:numel(r.cases)
    a=r.cases{k}.diagnostic.cases(1);
    if abs(a.P0-p)<1e-12, c=a; return; end
end
end

function verify_sources(root,sources)
% Existing processing changes do not change the stored integration.
for k=1:numel(sources)
    rel=strrep(sources(k).file,'\','/');
    if contains(rel,'compute_mvfp_revision_time_windows') || contains(rel,'mvfp_revision_refine_crossings'), continue; end
    assert(strcmp(sources(k).sha256,mvfp_streaming_sha256(fullfile(root,rel))), ...
        'Integration source changed since cache: %s',rel);
end
end

function verify_control_sources(root,sources)
names={'solve_mvfp_density_trace.m','solve_mvfp_density_trace_fixed_feedback.m','mvfp_rhs.m', ...
    'mvfp_rhs_fixed_feedback.m','mvfp_exponential_fit_flux.m','initial_density.m', ...
    'mvfp_grid_cache.m','mvfp_grid_cache_from_x.m','normalize_density.m'};
for j=1:numel(names)
    rel=['function/' names{j}]; idx=find(strcmp({sources.file},rel),1);
    assert(~isempty(idx) && strcmp(sources(idx).sha256,mvfp_streaming_sha256(fullfile(root,rel))), ...
        'Control integration source changed since cache: %s',rel);
end
end

function j=make_job(p,T,dt,solver,par,out,tag)
j=struct('P0',p,'T',T,'dt',dt,'solver',solver,'par',par,'keepDensity',false);
folder=fullfile(out,'cases'); if ~isfolder(folder),mkdir(folder);end
j.path=fullfile(folder,sprintf('%s_P%0.8f_Nx%d.mat',tag,p,solver.xGridN));
end

function results=run_jobs(jobs,cfg)
results=cell(size(jobs));
if isempty(jobs),return;end
if cfg.useParallel && license('test','Distrib_Computing_Toolbox')
    pool=gcp('nocreate');
    if isempty(pool), parpool('Processes',cfg.numWorkers); end
    parfor k=1:numel(jobs), results{k}=run_one(jobs{k}); end
else
    for k=1:numel(jobs), results{k}=run_one(jobs{k}); end
end
end

function s=run_one(job)
signature=rmfield(job,'path');
if isfile(job.path)
    old=load(job.path,'s','signature');
    assert(abs(old.signature.P0-signature.P0)<1e-12 && ...
        isequaln(rmfield(old.signature,'P0'),rmfield(signature,'P0')),'Local cache config differs.');
    s=old.s; return
end
fprintf('Integrating P0=%g Nx=%d T=%g reltol=%g\n',job.P0,job.solver.xGridN,job.T,job.solver.odeRelTol);
w=tic; tr=solve_mvfp_density_trace(job.par,job.P0,0:job.dt:job.T,job.solver);
s=pack_trace(tr,job.par); s.wallSeconds=toc(w);
assert(max(abs(s.rawMass-1))<1e-6 && max(s.boundaryMass)<1e-6 && min(s.rawMinDensity)>-1e-10, ...
    'Numerical diagnostics failed P0=%g',job.P0);
if job.keepDensity, s.rhoRaw=tr.rhoRaw; end
save(job.path,'s','signature','-v7.3');
fprintf('Completed P0=%g in %.2fs; P(T)=%.9g\n',job.P0,s.wallSeconds,s.PFixed(end));
end

function s=pack_trace(tr,par)
r=mvfp_trace_moment_readouts(tr,par);
land=mvfp_trace_dynamic_landscape(tr,r,par,.2:.2:.8);
s=struct('x',tr.x,'t',tr.t(:).','PFixed',r.PFixed(:).','F',r.F(:).', ...
    'xL',land.xL(:).','xs',land.xs(:).','xR',land.xR(:).','xSingle',land.xSingle(:).', ...
    'barrierRL',land.barrierRL(:).','barrierLR',land.barrierLR(:).', ...
    'rawMass',tr.rawMass(:).','boundaryMass',tr.rawBoundaryMass(:).','rawMinDensity',tr.rawMinDensity(:).');
end

function s=from_control(c)
s=struct('x',c.x,'t',c.t(:).','PFixed',c.PFixed(:).','F',c.F(:).', ...
    'xL',c.xL(:).','xs',c.xs(:).','xR',c.xR(:).','xSingle',c.xSingle(:).', ...
    'barrierRL',c.barrierRL(:).','barrierLR',c.barrierLR(:).', ...
    'rawMass',trapz(c.x,c.rhoRaw,1),'boundaryMass',c.boundaryMass(:).','rawMinDensity',c.rawMinDensity(:).');
end

function t=trace_table(p,nx,key,s)
n=numel(s.t);
t=table(repmat(p,n,1),repmat(nx,n,1),repmat(string(key),n,1),s.t(:),s.PFixed(:),s.F(:), ...
    s.barrierRL(:),s.barrierLR(:),s.rawMass(:),abs(s.rawMass(:)-1),s.boundaryMass(:),s.rawMinDensity(:), ...
    s.xL(:),s.xs(:),s.xR(:),s.xSingle(:), ...
    'VariableNames',{'P0','Nx','control','t','PFixed','F','barrierRightToLeft','barrierLeftToRight', ...
    'rawMass','rawMassError','boundaryMass','rawMinDensity','xLeftMinimum','xSaddle','xRightMinimum','xSingleMinimum'});
end

function rows=event_rows(p,s)
last=find(isfinite(s.barrierRL),1,'last'); bracket=[nan nan];
if ~isempty(last)
    bracket(1)=s.t(last); if last<numel(s.t), bracket(2)=s.t(last+1); end
end
rows=struct([]);
for q=[.2 .4 .5 .6 .8]
    e=mvfp_revision_first_arrival(s.t,s.PFixed,q);
    rows=[rows,struct('P0',p,'q',q,'time',e.time,'status',string(e.status), ...
        'lastBistableLower',bracket(1),'lastBistableUpper',bracket(2))]; %#ok<AGROW>
end
end

function d=diag_row(p,key,s)
d=struct('P0',p,'control',string(key),'PFixedT',s.PFixed(end), ...
    'maxMassError',max(abs(s.rawMass-1)),'maxBoundaryMass',max(s.boundaryMass),'minRawDensity',min(s.rawMinDensity));
end

function write_json(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>0); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value,'ConvertInfAndNaN',true));
end
