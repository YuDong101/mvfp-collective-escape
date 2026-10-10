function run_auxiliary(root)
addpath(fullfile(root,'function'));maxNumCompThreads(1);set(0,'DefaultFigureVisible','off');
base=mvfp_paper_revision_config(root);par=base.parA;par.D=.01;
opts=struct('P0',0,'T',100,'useParallel',false,'steadyFGrid',linspace(-8,8,1601), ...
    'steadyRootTol',1e-8,'xSigmaBifurcationScan',unique([10.^linspace(-2,1,96),linspace(2,6,81)]), ...
    'xSigmaBifurcationFocusScan',linspace(1,4,121), ...
    'mvfp',struct('xMax',1,'xGridN',1001,'initStd',.02,'thresholdMode','fixed','odeRelTol',1e-4,'odeAbsTol',1e-8), ...
    'freeEnergyBasinDensityTol',.04,'freeEnergyBasinFeedbackTol',.005,'freeEnergyBasinUnstableDistanceTol',.1, ...
    'freeEnergyBasinD',.01,'freeEnergyBasinSigma1',2);
folder=fullfile(root,'data','fig45');if ~isfolder(folder),mkdir(folder);end
file=fullfile(folder,'branch_results.mat');
if isfile(file),q=load(file,'xb');xb=q.xb;else
    xb=compute_mvfp_x_bifurcation_vs_sigma(par,opts);save(file,'xb','opts','par','-v7.3');
end
fprintf('Fresh stationary branches computed: %d roots\n',numel(xb.F));
opts.freeEnergyBasinTimes=0:.5:100;opts.freeEnergyBasinDissipationTimes=unique([0:.5:100,0:.01:10]);
opts.freeEnergyBasinP0List=[0 .05 .06 .065 .068 .0695 .07 .0703 .0705 .0707 .071 .0713 .0715 .072 .075 .08 .1];
file=fullfile(folder,'basin_trajectories.mat');
if isfile(file),q=load(file,'diag');diag=q.diag;else
    diag=compute_mvfp_free_energy_basin_diagnostics(par,opts);save(file,'diag','opts','par','-v7.3');
end
compact=diag;compact.cases=rmfield(compact.cases,'trace');compact.dissipationCase=rmfield(compact.dissipationCase,'trace');
subset=struct('source','fresh numerical solves','source_sha256',mvfp_streaming_sha256(file), ...
    'meta',struct('createdAt',char(datetime('now')),'runID','cp_fresh_recompute_2026_10_10'), ...
    'par',par,'xBifurcationVsSigma',xb,'freeEnergyBasinDiagnostics',compact);
write_json(fullfile(folder,'source_subset.json'),subset);
clear diag compact subset;
opts.freeEnergyBasinP0List=[0 .071 .1];
opts.freeEnergyBasinTimes=unique([0:.1:100,0:.005:8,35:.005:46]);
opts.freeEnergyBasinDissipationTimes=opts.freeEnergyBasinTimes;
file=fullfile(folder,'dense_trajectories.mat');
if isfile(file),q=load(file,'diag');diag=q.diag;else
    diag=compute_mvfp_free_energy_basin_diagnostics(par,opts);save(file,'diag','opts','-v7.3');
end
compact=struct('par',diag.par,'mvfp',opts.mvfp,'branches',diag.branches, ...
    'cases',rmfield(diag.cases,'trace'),'xGridN',numel(diag.x),'matlabVersion',version);
write_json(fullfile(folder,'dense_trajectories.json'),compact);clear diag compact;
fprintf('Fresh Figure 5 trajectories computed\n');
write_json(fullfile(root,'logs','auxiliary_complete.json'),struct('status','computed-awaiting-audit','finishedAt',char(datetime('now'))));
end

function write_json(path,value)
fid=fopen(path,'w');assert(fid>=0);fwrite(fid,jsonencode(value));fclose(fid);
end
