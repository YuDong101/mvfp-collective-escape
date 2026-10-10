% Figures 1--5 and supporting local theory. Use a short path on Windows.
% Enabled stages compute/resume their own outputs; disabled stages are kept.
runBranchesAndFreeEnergy=true;
runCaseA=true;
runCaseB=true;
finalizeLocalData=true;

root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'function'),fullfile(root,'recompute'));
assert(strcmp(which('solve_mvfp_density_trace'),fullfile(root,'function','solve_mvfp_density_trace.m')));
for folder={'logs','data','provenance'}
    if ~isfolder(fullfile(root,folder{1})),mkdir(fullfile(root,folder{1}));end
end
inputs=jsondecode(fileread(fullfile(root,'provenance','input_config.json')));
base=mvfp_paper_revision_config(root);
assert(isequaln(inputs.parA,base.parA) && isequaln(inputs.parB,base.parB));
fid=fopen(fullfile(root,'provenance','rerun_started.json'),'w');assert(fid>=0);
fwrite(fid,jsonencode(struct('startedAt',char(datetime('now')),'matlab',version, ...
    'scope','Figures 1--5 and local theory','layout','shared-function-library')));fclose(fid);
if runBranchesAndFreeEnergy,run_auxiliary(root);end
if runCaseA,run_core_a(root);end
if runCaseB,run_core_b(root);end
if finalizeLocalData,finalize_local_data(root);end
fprintf('Numerical stages complete. Continue with the Python commands in README.md.\n');
