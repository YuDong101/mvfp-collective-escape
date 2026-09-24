function results = run_minimal_reproduction
%RUN_MINIMAL_REPRODUCTION Minimal, standalone MVFP reproduction.
%
%   RESULTS = RUN_MINIMAL_REPRODUCTION computes two representative figures
%   for the collective-escape mechanism:
%     1) density evolution, self-consistent potential, feedback, and readout;
%     2) final left fraction P(T) and first threshold-arrival time versus P0.
%
% The script uses only base MATLAB. It does not depend on the parent project's
% function/ directory. Edit the parameter block below to explore other cases.
%
% Example from PowerShell:
%   matlab -batch "run('run_minimal_reproduction.m')"
%
% Outputs are written beside this file in figures/ and results/.

%% Editable experiment configuration
cfg = struct();
cfg.u = 4.0;                         % Bare quadratic coefficient.
cfg.v = 20.0;                        % Bare quartic coefficient.
cfg.bias = -0.5;                     % Constant feedback bias.
cfg.D = 0.005;                       % Diffusion/noise strength.
cfg.sigma1 = 0.5;                    % Pairwise feedback strength.
cfg.xMax = 1.0;                      % Spatial interval [-xMax,xMax].
cfg.xGridN = 241;                    % Grid points; increase for a finer run.
cfg.initStd = 0.035;                 % Initial packet width.
cfg.T = 60;                          % Finite observation time.
cfg.snapshotTimes = [0 5 15 30 60];
cfg.snapshotP0 = 0.08;               % Representative near-transition initial P0.
cfg.readoutP0 = [0 0.04 0.08 0.16];  % Curves shown in the mechanism panel.
cfg.P0Scan = linspace(0,1,21);        % Small scan for the second figure.
cfg.readoutThreshold = 0.5;           % Left-fraction threshold for arrival time.
cfg.odeRelTol = 2e-5;
cfg.odeAbsTol = 1e-8;
cfg.seed = 11;                        % Kept for reproducible future extensions.

packageRoot = fileparts(mfilename('fullpath'));
figureDir = fullfile(packageRoot,'figures');
resultDir = fullfile(packageRoot,'results');
if ~isfolder(figureDir), mkdir(figureDir); end
if ~isfolder(resultDir), mkdir(resultDir); end
rng(cfg.seed,'twister');

x = linspace(-cfg.xMax,cfg.xMax,cfg.xGridN).';
tSnapshot = cfg.snapshotTimes(:).';
trace = simulate_case(x,cfg,cfg.snapshotP0,tSnapshot);

% Figure 1: dynamic feedback mechanism.
fig1 = figure('Color','w','Name','Minimal MVFP core','Visible','off');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
colors = lines(numel(tSnapshot));
for k = 1:numel(tSnapshot)
    plot(x,trace.rho(:,k),'LineWidth',1.5,'Color',colors(k,:));
end
xlabel('x'); ylabel('\rho(x,t)'); title(sprintf('Density, P_0 = %.3g',cfg.snapshotP0));
legend(compose('t = %.3g',tSnapshot),'Location','best'); grid on; box on;
nexttile; hold on;
for k = 1:numel(tSnapshot)
    U = effective_potential(x,cfg,trace.F(k));
    plot(x,U-min(U),'LineWidth',1.5,'Color',colors(k,:));
end
xlabel('x'); ylabel('U_{eff}(x,t)-min U'); title('Self-consistent potential');
legend(compose('t = %.3g',tSnapshot),'Location','best'); grid on; box on;
nexttile; hold on;
plot(trace.t,trace.F,'k-','LineWidth',1.8);
plot(trace.t,trace.M1,'--','LineWidth',1.5);
yline(0,'Color',[.5 .5 .5]); xlabel('t'); ylabel('value');
title('Feedback field and first moment'); legend('F(t)','M_1(t)','Location','best'); grid on; box on;
nexttile; hold on;
for j = 1:numel(cfg.readoutP0)
    tr = simulate_case(x,cfg,cfg.readoutP0(j),tSnapshot);
    plot(tr.t,tr.P,'LineWidth',1.5,'DisplayName',sprintf('P_0 = %.3g',cfg.readoutP0(j)));
end
xlabel('t'); ylabel('P_{fixed}(t)'); ylim([0 1]);
title('Finite-time left readout'); legend('Location','best'); grid on; box on;
sgtitle('Minimal MVFP collective-escape mechanism');
export_figure(fig1,fullfile(figureDir,'minimal_core'));
close(fig1);

% Figure 2: a compact P0 scan and threshold-arrival diagnostic.
P0 = cfg.P0Scan(:);
PT = nan(size(P0)); arrivalTime = nan(size(P0));
scanTimes = linspace(0,cfg.T,31);
for j = 1:numel(P0)
    tr = simulate_case(x,cfg,P0(j),scanTimes);
    PT(j) = tr.P(end);
    arrivalTime(j) = first_crossing(tr.t,tr.P,cfg.readoutThreshold);
end
fig2 = figure('Color','w','Name','Minimal MVFP readout and MFPT','Visible','off');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(P0,PT,'o-','LineWidth',1.6,'MarkerFaceColor',[.2 .4 .8]);
xlabel('initial left fraction P_0'); ylabel(sprintf('P_{fixed}(T), T = %.3g',cfg.T));
title('Finite-time readout curve'); ylim([0 1]); grid on; box on;
nexttile; semilogy(P0,arrivalTime,'o-','LineWidth',1.6,'MarkerFaceColor',[.8 .3 .2]);
xlabel('initial left fraction P_0'); ylabel('first time P_{fixed}(t) reaches 0.5');
title('Threshold-arrival diagnostic'); grid on; box on;
sgtitle('Minimal P_0 response and switching-time summary');
export_figure(fig2,fullfile(figureDir,'minimal_readout_mfpt'));
close(fig2);

results = struct('cfg',cfg,'x',x,'snapshot',trace,'P0',P0,'PT',PT, ...
    'arrivalTime',arrivalTime,'threshold',cfg.readoutThreshold);
save(fullfile(resultDir,'minimal_summary.mat'),'results','-v7');
writematrix([P0 PT arrivalTime],fullfile(resultDir,'minimal_readout_summary.csv'));
fprintf('Minimal MVFP reproduction finished.\n');
fprintf('Figures: %s\n',figureDir);
fprintf('Results: %s\n',resultDir);
end

function trace = simulate_case(x,cfg,P0,times)
times = times(:).';
rho0 = P0*gaussian_packet(x,-0.78,cfg.initStd) + ...
    (1-P0)*gaussian_packet(x,0.78,cfg.initStd);
rho0 = normalize_density(x,rho0);
if numel(times) == 1 && times == 0
    tOut = 0; Y = rho0.';
else
    odeOpts = odeset('RelTol',cfg.odeRelTol,'AbsTol',cfg.odeAbsTol, ...
        'NonNegative',1:numel(x));
    rhs = @(~,rho) mvfp_rhs_minimal(x,rho,cfg);
    [tOut,Y] = ode15s(rhs,times,rho0,odeOpts);
end
rho = zeros(numel(x),numel(times));
P = zeros(1,numel(times)); M1 = zeros(1,numel(times)); M2 = zeros(1,numel(times)); F = zeros(1,numel(times));
for k = 1:numel(times)
    [~,idx] = min(abs(tOut-times(k)));
    r = normalize_density(x,max(Y(idx,:).',0));
    rho(:,k) = r;
    P(k) = min(max(trapz(x,r.*(x<=0)),0),1);
    M1(k) = trapz(x,x.*r);
    M2(k) = trapz(x,x.^2.*r);
    F(k) = cfg.sigma1*M1(k) + cfg.bias;
end
trace = struct('x',x,'t',times,'rho',rho,'P',P,'M1',M1,'M2',M2,'F',F,'P0',P0);
end

function drho = mvfp_rhs_minimal(x,rho,cfg)
rho = max(rho(:),0);
dx = x(2)-x(1);
mass = trapz(x,rho);
if mass <= 0 || ~isfinite(mass), mass = 1; end
M1 = trapz(x,x.*rho)/mass;
F = cfg.sigma1*M1 + cfg.bias;
b = (cfg.u-cfg.sigma1)*x - cfg.v*x.^3 + F;
bFace = 0.5*(b(1:end-1)+b(2:end));
rhoL = rho(1:end-1); rhoR = rho(2:end);
rhoUp = rhoL;
mask = bFace < 0; rhoUp(mask) = rhoR(mask);
J = bFace.*rhoUp - cfg.D*(rhoR-rhoL)/dx;
Jfull = [0;J;0];
dr = -(Jfull(2:end)-Jfull(1:end-1))/dx;
dr = dr - mean(dr); % maintain the discrete mass to solver tolerance.
drho = dr;
end

function rho = gaussian_packet(x,mu,sigma)
rho = exp(-0.5*((x-mu)/sigma).^2)/(sqrt(2*pi)*sigma);
end

function rho = normalize_density(x,rho)
rho = max(rho(:),0);
mass = trapz(x,rho);
if mass <= 0 || ~isfinite(mass), error('MinimalMVFP:InvalidDensity','Density mass is invalid.'); end
rho = rho/mass;
end

function U = effective_potential(x,cfg,F)
U = -0.5*(cfg.u-cfg.sigma1)*x.^2 + 0.25*cfg.v*x.^4 - F*x;
end

function crossing = first_crossing(t,y,threshold)
idx = find(y >= threshold,1,'first');
if isempty(idx), crossing = NaN; return; end
if idx == 1, crossing = t(1); return; end
dy = y(idx)-y(idx-1);
if dy <= 0, crossing = t(idx); else
    crossing = t(idx-1) + (threshold-y(idx-1))*(t(idx)-t(idx-1))/dy;
end
end

function export_figure(fig,fileStem)
exportgraphics(fig,[fileStem '.png'],'Resolution',180);
exportgraphics(fig,[fileStem '.pdf'],'ContentType','vector');
end
