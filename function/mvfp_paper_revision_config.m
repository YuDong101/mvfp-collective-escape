function cfg = mvfp_paper_revision_config(projectRoot)
% 独立实验设置；不改变 mvfp_prl_analysis 的用户参数。
if nargin<1, projectRoot=fileparts(fileparts(mfilename('fullpath'))); end
cfg.projectRoot=projectRoot;
cfg.outputDir=fullfile(projectRoot,'outputs','pairwise_paper_revision_v2');
cfg.parA=struct('u',4,'v',20,'bias',-.5,'D',.005,'sigma1',.5,'sigma2',0);
cfg.parB=struct('u',4,'v',20,'bias',-.5,'D',.01,'sigma1',2,'sigma2',0);
cfg.mvfp=struct('xMax',1,'xGridN',1001,'initStd',.02, ...
    'thresholdMode','fixed','odeRelTol',1e-6,'odeAbsTol',1e-10, ...
    'massDiscretization','nodal-control-volume','storeRawDensity',true, ...
    'frozenFieldMode','actual-density');
cfg.T=100;
cfg.controlP0=[.060,.0647,.070,.080];
cfg.outputTimeStep=.5; % 首轮输出；尚未认证事件插值误差。
cfg.steadyFGrid=linspace(-3,3,1201);
cfg.steadyRootTol=1e-10;
cfg.pilotGridN=[1001,2001,4001];
cfg.pilotP0A=[.05,.07,.09];
cfg.pilotP0B=.0704; % 历史边界附近试点，不是新版 Pb。
cfg.pilotMaxStageSeconds=900; % 案例间检查预算；不截断正在积分的轨迹。
cfg.pilotMaxLastCaseSeconds=180; % 显著变慢时先评估算法，再放大网格。
cfg.stageConfig=struct(); % 新阶段必须显式提供配置；不继承首批粗网格。
cfg.dispatchOnly=false; % true 只检查新阶段入口参数，不计算、不写文件、不认证科学门槛。
end
