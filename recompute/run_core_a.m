function run_core_a(root)
addpath(fullfile(root,'function'));addpath(fullfile(root,'recompute'));
set(0,'DefaultFigureVisible','off');maxNumCompThreads(1);
base=mvfp_paper_revision_config(root);base.mvfp.fluxScheme='exponential-fit';base.mvfp.xGridN=2001;
data=fullfile(root,'outputs','pairwise_paper_revision_v2');
cfg=base;cfg.outputTimeStep=.1;cfg.outputDir=fullfile(data,'fitted_flux_controls');
run_mvfp_paper_revision_stage(cfg,'controls');
cfg.mvfp.xGridN=4001;cfg.controlP0=[.06 .08];
cfg.outputDir=fullfile(data,'fourth_batch','figures123','fine_controls_nx4001');
run_mvfp_paper_revision_stage(cfg,'controls');
cfg=struct('projectRoot',root,'outputDir',fullfile(data,'second_batch','time_windows_A','nx2001'), ...
    'par',base.parA,'mvfp',base.mvfp,'times',[50 100 200 400],'P0Values',0:.01:.12, ...
    'outputTimeStep',.5,'bracketWidth',1.25e-4,'maxStageSeconds',86400,'maxCaseSeconds',7200);
cfg.mvfp.odeRelTol=1e-8;cfg.mvfp.odeAbsTol=1e-12;
coarse=compute_mvfp_revision_time_windows(cfg);
points=[0 .04 .12];
for k=1:numel(coarse.thresholds)
    points=[points coarse.thresholds(k).P0Lower-.001 coarse.thresholds(k).P0Upper+.001]; %#ok<AGROW>
end
cfg.P0Values=unique(max(0,min(.12,points)));cfg.mvfp.xGridN=4001;
cfg.outputDir=fullfile(data,'second_batch','time_windows_A','nx4001');
fine=compute_mvfp_revision_time_windows(cfg);
fig=struct('projectRoot',root,'outputDir',fullfile(root,'data','fig123'), ...
    'runFigure1',true,'runFigure2',true,'runFigure3',true,'smoke',false,'useParallel',false,'numWorkers',1);
result=compute_mvfp_approved_revision_fig123(fig);
save(fullfile(root,'data','fig123','fresh_result.mat'),'result','-v7.3');
fid=fopen(fullfile(root,'logs','core_a_complete.json'),'w');
fwrite(fid,jsonencode(struct('status','computed-awaiting-audit','finishedAt',char(datetime('now')))));fclose(fid);
end
