function result = compute_mvfp_revision_time_windows(cfg)
% Exploratory multi-time readouts and all sampled crossings, with checkpoints.
required={'projectRoot','outputDir','par','mvfp','times','P0Values'};
for k=1:numel(required)
    assert(isfield(cfg,required{k}),'Missing time-window configuration: %s',required{k});
end
defaults=struct('q',.5,'outputTimeStep',.5,'bracketWidth',1.25e-4, ...
    'maxMassError',1e-6,'maxBoundaryMass',1e-6,'minRawDensity',-1e-10, ...
    'maxCaseSeconds',180,'maxStageSeconds',1800,'maxRefinements',20, ...
    'reusePhysicsFromPriorProcessing',false);
cfg=merge_struct_fields(defaults,cfg);
if ~isfield(cfg.mvfp,'thresholdMode') || ~strcmp(cfg.mvfp.thresholdMode,'fixed') || cfg.q~=.5
    error('MVFP:Revision:InvalidReadoutProtocol', ...
        'Time-window t50 metadata requires thresholdMode=''fixed'' and q=0.5.');
end
assert(cfg.par.sigma2==0 && strcmp(cfg.mvfp.fluxScheme,'exponential-fit') && ...
    strcmp(cfg.mvfp.massDiscretization,'nodal-control-volume'), ...
    'Revision time-window runs require the verified pairwise fitted/nodal protocol.');
assert(all(isfinite(cfg.times)) && all(cfg.times>0) && cfg.outputTimeStep>0 && ...
    cfg.bracketWidth>0 && all(cfg.P0Values>=0 & cfg.P0Values<=1), 'Invalid scan configuration.');
cfg.times=unique(cfg.times(:).'); cfg.P0Values=unique(cfg.P0Values(:).');
assert(~isempty(cfg.times) && numel(cfg.P0Values)>=2,'At least two initial weights are required.');
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
folder=fullfile(cfg.outputDir,'cases'); if ~exist(folder,'dir'), mkdir(folder); end
sourceHashes=hash_integration_sources(cfg.projectRoot);
processingSourceHashes=hash_processing_sources(cfg.projectRoot);
times=unique([0:cfg.outputTimeStep:max(cfg.times),cfg.times]);
solver=cfg.mvfp; solver.storeRawDensity=true;
entries=struct([]); clockStart=tic; lastSeconds=0;
result=struct('config',cfg,'sourceHashes',sourceHashes, ...
    'processingSourceHashes',processingSourceHashes,'createdAt',char(datetime('now')), ...
    'status','exploratory-general-convergence-not-yet-certified', ...
    'readoutDefinition','fixed x=0 occupancy; P0 is packet preparation weight', ...
    'samplingLimitation','All crossings resolved on sampled grid; hidden crossings between samples are not excluded.');
for P0=cfg.P0Values
    evaluate(P0);
end
coarse=sort_entries(entries);
[~,~,~,allBrackets]=mvfp_revision_refine_crossings([coarse.P0], ...
    vertcat(coarse.PAtTimes).',cfg.q,[coarse.reliable],@evaluate_readouts, ...
    cfg.bracketWidth,cfg.maxRefinements);
coarse=sort_entries(entries);
thresholds=struct([]); windows=struct([]);
for j=1:numel(cfg.times)
    y=arrayfun(@(e)e.PAtTimes(j),coarse); reliable=[coarse.reliable];
    brackets=allBrackets{j};
    label='sampled-crossings';
    if isempty(brackets)
        if ~all(reliable), label='unresolved-samples';
        elseif all(y<cfg.q), label='all-sampled-below';
        elseif all(y>cfg.q), label='all-sampled-above';
        end
    end
    win=struct('T',cfg.times(j),'q',cfg.q,'nSampledCrossings',numel(brackets),'status',string(label));
    if isempty(windows), windows=win; else, windows(end+1)=win; end %#ok<AGROW>
    for k=1:numel(brackets)
        b=brackets(k);
        lo=evaluate(b.lower); hi=evaluate(b.upper);
        row=struct('T',cfg.times(j),'q',cfg.q,'crossingIndex',k, ...
            'P0Lower',b.lower,'P0Upper',b.upper,'midpoint',.5*(b.lower+b.upper), ...
            'width',b.upper-b.lower,'PLower',b.yLower,'PUpper',b.yUpper, ...
            't50Lower',lo.firstArrival.time,'t50Upper',hi.firstArrival.time, ...
            'sampledRecrossLower',lo.sampledRecrossAfterFirstHit, ...
            'sampledRecrossUpper',hi.sampledRecrossAfterFirstHit, ...
            'status',string(b.status),'Nx',solver.xGridN,'xMax',solver.xMax);
        if isempty(thresholds), thresholds=row; else, thresholds(end+1)=row; end %#ok<AGROW>
    end
end
result.entries=sort_entries(entries); result.thresholds=thresholds;
result.windowStatus=windows; result.wallSeconds=toc(clockStart);
result.totalIntegrationSeconds=sum([entries.wallSeconds]);
save(fullfile(cfg.outputDir,'time_windows.mat'),'result','-v7.3');
writetable(struct2table(thresholds),fullfile(cfg.outputDir,'thresholds_summary.csv'));
writetable(struct2table(windows),fullfile(cfg.outputDir,'window_status.csv'));
write_samples();
fprintf('Saved time windows: %d unique traces, %s\n',numel(entries),cfg.outputDir);

    function e=evaluate(P0)
        if ~isempty(entries)
            found=find([entries.P0]==P0,1);
            if ~isempty(found), e=entries(found); return; end
        end
        signature=struct('par',cfg.par,'solver',solver,'P0',P0,'times',times,'sourceHashes',sourceHashes, ...
            'readoutTimes',cfg.times,'q',cfg.q,'diagnosticLimits',[cfg.maxMassError cfg.maxBoundaryMass cfg.minRawDensity]);
        path=fullfile(folder,[mvfp_value_sha256(signature) '.mat']);
        if toc(clockStart)>cfg.maxStageSeconds || lastSeconds>cfg.maxCaseSeconds
            error('MVFP:Revision:TimeBudget','Stopped between cases; checkpoints are retained.');
        end
        reusedPath='';
        if isfile(path)
            reusedPath=path;
        elseif cfg.reusePhysicsFromPriorProcessing
            reusedPath=find_compatible_prior_case(folder,signature,cfg.projectRoot);
        end
        if ~isempty(reusedPath)
            cached=load(reusedPath,'entry'); e=cached.entry;
            if strcmp(reusedPath,path)
                assert(isequaln(e.signature,signature),'Cached trace protocol differs.');
            else
                assert(signatures_match_physics(e.signature,signature), ...
                    'Cached trace numerical integration protocol differs.');
            end
            e.firstArrival=mvfp_revision_first_arrival(e.t,e.PFixed,cfg.q);
            hit=find(e.PFixed>=cfg.q,1); e.sampledRecrossAfterFirstHit=false;
            if ~isempty(hit), e.sampledRecrossAfterFirstHit=any(e.PFixed(hit:end)<cfg.q); end
        else
            fprintf('Time windows Nx=%d P0=%.9g Tmax=%g\n',solver.xGridN,P0,max(times));
            watch=tic; trace=solve_mvfp_density_trace(cfg.par,P0,times,solver); elapsed=toc(watch);
            g=mvfp_grid_cache(solver); [~,idx]=ismember(cfg.times,trace.t);
            first=mvfp_revision_first_arrival(trace.t,trace.P,cfg.q);
            hit=find(trace.P>=cfg.q,1); recross=false;
            if ~isempty(hit), recross=any(trace.P(hit:end)<cfg.q); end
            rawP=trapz(trace.x(trace.x<=0),trace.rhoRaw(trace.x<=0,:),1);
            e=struct('P0',P0,'signature',signature,'t',trace.t,'PFixed',trace.P(:).', ...
                'PRaw',rawP,'PAtTimes',trace.P(idx).','rhoAtTimes',trace.rhoRaw(:,idx), ...
                'x',trace.x,'rho0',trace.rhoRaw(:,1),'M1',g.wx.'*trace.rhoRaw, ...
                'rawMass',trace.rawMass,'firstArrival',first,'sampledRecrossAfterFirstHit',recross, ...
                'maxMassError',max(abs(trace.rawMass-1)), ...
                'maxBoundaryMass',max(trace.rawBoundaryMass),'minRawDensity',min(trace.rawMinDensity), ...
                'maxReadoutNormalizationChange',max(abs(rawP-trace.P(:).')), ...
                'wallSeconds',elapsed,'reliable',false);
            e.reliable=e.maxMassError<=cfg.maxMassError && e.maxBoundaryMass<=cfg.maxBoundaryMass ...
                && e.minRawDensity>=cfg.minRawDensity && all(isfinite(e.PFixed));
            entry=e; partial=[path '.partial.mat']; %#ok<NASGU>
            save(partial,'entry','-v7.3'); movefile(partial,path,'f');
        end
        e.executionSourceHashes=e.signature.sourceHashes;
        e.processingSourceHashes=processingSourceHashes;
        e.reusedFromPriorProcessing=~isempty(reusedPath) && ~strcmp(reusedPath,path);
        lastSeconds=e.wallSeconds;
        if isempty(entries), entries=e; else, entries(end+1)=e; end
        write_samples();
    end

    function [values,isReliable]=evaluate_readouts(P0)
        e=evaluate(P0); values=e.PAtTimes(:); isReliable=e.reliable;
    end

    function write_samples()
        sampleRows=struct([]);
        ordered=sort_entries(entries);
        for ii=1:numel(ordered)
            e=ordered(ii);
            for jj=1:numel(cfg.times)
                r=struct('P0',e.P0,'T',cfg.times(jj),'PFixed',e.PAtTimes(jj), ...
                    't50',e.firstArrival.time,'eventStatus',string(e.firstArrival.status), ...
                    'eventObservationEnd',max(times),'sampledRecross',e.sampledRecrossAfterFirstHit, ...
                    'reliable',e.reliable,'maxMassError',e.maxMassError,'maxBoundaryMass',e.maxBoundaryMass, ...
                    'minRawDensity',e.minRawDensity,'wallSeconds',e.wallSeconds);
                if isempty(sampleRows), sampleRows=r; else, sampleRows(end+1)=r; end %#ok<AGROW>
            end
        end
        writetable(struct2table(sampleRows),fullfile(cfg.outputDir,'readouts_summary.csv'));
    end
end

function ordered=sort_entries(entries)
if isempty(entries), ordered=entries; return; end
[~,idx]=sort([entries.P0]); ordered=entries(idx);
end

function sources=hash_integration_sources(root)
entry=fullfile(root,'function','solve_mvfp_density_trace.m');
files=matlab.codetools.requiredFilesAndProducts(entry);
files=files(startsWith(files,fullfile(root,'function')));
sources=struct('file',{},'sha256',{});
for k=1:numel(files)
    [~,name,ext]=fileparts(files{k});
    sources(end+1)=struct('file',['function/' name ext], ...
        'sha256',mvfp_streaming_sha256(files{k})); %#ok<AGROW>
end
[~,order]=sort({sources.file}); sources=sources(order);
end

function sources=hash_processing_sources(root)
names={'compute_mvfp_revision_time_windows.m','mvfp_revision_refine_crossings.m', ...
    'mvfp_revision_crossing_intervals.m','mvfp_revision_first_arrival.m'};
sources=struct('file',{},'sha256',{});
for k=1:numel(names)
    path=fullfile(root,'function',names{k});
    sources(end+1)=struct('file',['function/' names{k}], ...
        'sha256',mvfp_streaming_sha256(path)); %#ok<AGROW>
end
end

function path=find_compatible_prior_case(folder,signature,root)
path=''; files=dir(fullfile(folder,'*.mat'));
for k=1:numel(files)
    candidate=load(fullfile(files(k).folder,files(k).name),'entry');
    if ~isfield(candidate,'entry') || ~isfield(candidate.entry,'signature') || ...
            ~signatures_match_physics(candidate.entry.signature,signature)
        continue;
    end
    if legacy_hash_difference_is_processing_only(candidate.entry.signature.sourceHashes,root)
        path=fullfile(files(k).folder,files(k).name); return;
    end
end
end

function tf=signatures_match_physics(a,b)
if ~isstruct(a) || ~isstruct(b) || ~isfield(a,'sourceHashes') || ~isfield(b,'sourceHashes')
    tf=false; return;
end
a=rmfield(a,'sourceHashes'); b=rmfield(b,'sourceHashes'); tf=isequaln(a,b);
end

function tf=legacy_hash_difference_is_processing_only(oldHashes,root)
allowed={'function/compute_mvfp_revision_time_windows.m','function/mvfp_revision_refine_crossings.m'};
tf=isstruct(oldHashes) && all(isfield(oldHashes,{'file','sha256'}));
for k=1:numel(oldHashes)
    rel=strrep(oldHashes(k).file,'\','/');
    path=fullfile(root,strrep(rel,'/',filesep));
    if ~isfile(path), tf=tf && ismember(rel,allowed); continue; end
    if ~strcmp(oldHashes(k).sha256,mvfp_streaming_sha256(path)) && ~ismember(rel,allowed)
        tf=false; return;
    end
end
end
