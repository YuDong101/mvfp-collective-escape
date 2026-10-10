function finalize_local_data(root,figureRoot)
% Run after core A, core B and auxiliary computations have completed.
if nargin<2,figureRoot=root;end
addpath(fullfile(root,'function'));set(0,'DefaultFigureVisible','off');maxNumCompThreads(1);
for name={'core_a_complete.json','core_b_complete.json'}
    assert(isfile(fullfile(root,'logs',name{1})),'Required stage incomplete: %s',name{1});
end
assert(isfile(fullfile(figureRoot,'logs','auxiliary_complete.json')));
validate_mvfp_approved_revision_fig123(root,fullfile(root,'data','fig123'));
data=fullfile(root,'outputs','pairwise_paper_revision_v2');
folder=fullfile(figureRoot,'data','fig45');
q=load(fullfile(folder,'basin_trajectories.mat'),'diag','opts','par');
b=load(fullfile(data,'second_batch','basin_B_local','basin_B_results.mat'),'result');
br=b.result.brackets(1);interval=[br.lowerP0 br.upperP0];
extraOpts=q.opts;extraOpts.freeEnergyBasinP0List=interval;
extra=compute_mvfp_free_energy_basin_diagnostics(q.par,extraOpts);
save(fullfile(folder,'fresh_basin_endpoint_trajectories.mat'),'extra','extraOpts','-v7.3');
diag=q.diag;diag.cases=[diag.cases extra.cases];[diag.P0,order]=sort([diag.cases.P0]);diag.cases=diag.cases(order);
diag.selectedStableIndex=[diag.cases.selectedStableIndex];diag.convergedMask=[diag.cases.convergedToStableBranch];
diag.PbInterval=interval;diag.cases=rmfield(diag.cases,'trace');diag.dissipationCase=rmfield(diag.dissipationCase,'trace');
subset=jsondecode(fileread(fullfile(folder,'source_subset.json')));subset.freeEnergyBasinDiagnostics=diag;
xb=subset.xBifurcationVsSigma;par=subset.par;
combined=fullfile(folder,'combined_source.mat');save(combined,'diag','xb','par','-v7.3');
subset.source='combined_source.mat';subset.source_sha256=mvfp_streaming_sha256(combined);
fid=fopen(fullfile(folder,'source_subset.json'),'w');fwrite(fid,jsonencode(subset));fclose(fid);
res=fullfile(data,'third_batch','residence_feasibility');
copyfile(fullfile(res,'residence_prediction.csv'),fullfile(folder,'residence_prediction.csv'));
copyfile(fullfile(res,'residence_prediction.mat'),fullfile(folder,'residence_prediction.mat'));
fid=fopen(fullfile(figureRoot,'logs','local_validation_complete.json'),'w');
fwrite(fid,jsonencode(struct('status','local-results-validated','finishedAt',char(datetime('now')))));fclose(fid);
end
