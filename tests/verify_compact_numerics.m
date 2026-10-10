function report=verify_compact_numerics(root,outputDir,checkDependencies)
% Actual new calculations; no stored study results are loaded.
% The optional third argument is used only for comparison with the old snapshot.
if nargin<3,checkDependencies=true;end
root=strrep(root,'/',filesep);outputDir=strrep(outputDir,'/',filesep);
savedPath=path;cleanup=onCleanup(@()path(savedPath)); %#ok<NASGU>
restoredefaultpath;addpath(fullfile(root,'function'));
maxNumCompThreads(1);set(0,'DefaultFigureVisible','off');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(strcmp(which('solve_mvfp_density_trace'),fullfile(root,'function','solve_mvfp_density_trace.m')));
base=mvfp_paper_revision_config(root);cfg=base.mvfp;cfg.xGridN=201;
cfg.fluxScheme='upwind';times=0:.1:3;
upwind=solve_mvfp_density_trace(base.parA,.08,times,cfg);
assert(max(abs(upwind.rawMass-1))<1e-5 && all(upwind.threshold==0));
cfg.fluxScheme='exponential-fit';cfg.odeRelTol=1e-8;cfg.odeAbsTol=1e-12;
fitted=solve_mvfp_density_trace(base.parB,.073,times,cfg);
grid=mvfp_grid_cache(cfg);rho0=initial_density(grid.x,base.parB,.073,cfg.initStd,cfg);
field=base.parB.sigma1*trapz(grid.x,grid.x.*rho0)+base.parB.bias;
frozen=solve_mvfp_density_trace_fixed_feedback(base.parB,.073,times,cfg,field);
assert(max(abs(fitted.rawMass-1))<1e-6 && min(fitted.rawMinDensity)>-1e-10);
assert(max(abs(frozen.rawMass-1))<1e-6 && min(frozen.rawMinDensity)>-1e-10);
particle=struct('N',100,'nRealizations',3,'dt',.02,'xClip',8);
[pm,ps,pv]=langevin_final_left_fraction(base.parA,.08,2,particle,1011);
assert(all(isfinite(pv)) && abs(pm-mean(pv))<eps && abs(ps-std(pv))<eps);

% Exercise source bookkeeping, checkpoint creation and matched-density controls.
controlCfg=base;controlCfg.mvfp=cfg;controlCfg.T=3;
controlCfg.controlP0=.08;controlCfg.outputTimeStep=.1;
controlCfg.outputDir=fullfile(outputDir,'controls');
controls=run_mvfp_paper_revision_stage(controlCfg,'controls');
c=controls.cases{1}.diagnostic.cases(1);
assert(c.initialL1.fullVsFrozen<1e-12 && c.initialL1.fullVsNoCoupling<1e-12);
controlValues=[c.full.PFixed;c.frozenFeedback.PFixed;c.noCoupling.PFixed];

% Exercise the density export route without requiring a pre-existing cache.
export=compute_mvfp_approved_revision_fig123(struct('projectRoot',root, ...
    'outputDir',fullfile(outputDir,'export'),'smoke',true));
assert(export.smokePassed);
windowCfg=struct('projectRoot',root,'outputDir',fullfile(outputDir,'windows'), ...
    'par',base.parA,'mvfp',cfg,'times',[1 2],'P0Values',[0 .04], ...
    'outputTimeStep',.1,'maxStageSeconds',60,'maxCaseSeconds',30);
windows=compute_mvfp_revision_time_windows(windowCfg);
assert(numel(windows.windowStatus)==2);

% Stationary branches, free energy and the unstable eigenmode use different
% parts of the shared library and catch missing non-solver dependencies.
opts=struct('P0',0,'T',3,'useParallel',false,'steadyFGrid',linspace(-8,8,1601), ...
    'steadyRootTol',1e-8,'xSigmaBifurcationScan',[.5 2], ...
    'xSigmaBifurcationFocusScan',[1 2], ...
    'freeEnergyBasinP0List',[0 .071 .1],'freeEnergyBasinTimes',0:.1:3, ...
    'freeEnergyBasinDissipationTimes',0:.1:3,'freeEnergyBasinD',.01, ...
    'freeEnergyBasinSigma1',2,'mvfp',base.mvfp);
opts.mvfp.xGridN=201;
par=base.parA;par.D=.01;
branches=compute_mvfp_x_bifurcation_vs_sigma(par,opts);
energy=compute_mvfp_free_energy_basin_diagnostics(par,opts);
assert(~isempty(branches.F) && numel(energy.cases)==3);
spectrum=compute_mvfp_revision_local_spectrum(base.parB,struct('Nx',201,'runGrowth',true));
assert(spectrum.positiveModeVerified && spectrum.numericalPositiveCount==1);

report=struct('status','passed','no_precomputed_results_loaded',true, ...
    'gridN',201,'upwindMassError',max(abs(upwind.rawMass-1)), ...
    'fittedMassError',max(abs(fitted.rawMass-1)),'frozenMassError',max(abs(frozen.rawMass-1)), ...
    'particleRealizations',3,'matchedInitialDensity',true,'densityExport',export.smokePassed, ...
    'timeWindows',numel(windows.windowStatus),'stationaryRoots',numel(branches.F), ...
    'freeEnergyTrajectories',numel(energy.cases),'lambdaNx201',spectrum.lambda,'matlab',version);
if checkDependencies
    addpath(fullfile(root,'recompute'));
    entry={fullfile(root,'run_recomputation.m'),fullfile(root,'recompute','run_auxiliary.m'), ...
        fullfile(root,'recompute','run_core_a.m'),fullfile(root,'recompute','run_core_b.m'), ...
        fullfile(root,'recompute','finalize_local_data.m'), ...
        fullfile(root,'function','langevin_trace_realization.m'), ...
        fullfile(root,'tests','verify_compact_numerics.m')};
    [files,products]=matlab.codetools.requiredFilesAndProducts(entry);
    for k=1:numel(files)
        assert(startsWith(files{k},[root filesep],'IgnoreCase',true) || ...
            startsWith(files{k},[matlabroot filesep],'IgnoreCase',true), ...
            'External project dependency: %s',files{k});
    end
    report.dependencyFiles=files;report.requiredProducts=products;
end
numeric=struct('upwindP',upwind.P,'upwindRho',upwind.rho, ...
    'fittedP',fitted.P,'fittedRho',fitted.rho,'frozenP',frozen.P,'frozenRho',frozen.rho, ...
    'particleValues',pv,'controlValues',controlValues,'branchF',branches.F, ...
    'lambda',spectrum.lambda,'rightMode',spectrum.rightMode);
save(fullfile(outputDir,'quick_numerics.mat'),'numeric','report','-v7.3');
fid=fopen(fullfile(outputDir,'quick_report.json'),'w');assert(fid>=0);
fwrite(fid,jsonencode(report));fclose(fid);disp(report);
end
