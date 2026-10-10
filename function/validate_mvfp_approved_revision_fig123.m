function validation=validate_mvfp_approved_revision_fig123(root,out)
%VALIDATE_MVFP_APPROVED_REVISION_FIG123 Data/export and focused quadrature checks.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,out=fullfile(root,'outputs','pairwise_figure_revision_2026_09_13','data','fig123');end
cfg=mvfp_paper_revision_config(root); p=cfg.parA; sol=cfg.mvfp;
sol.xGridN=2001; sol.fluxScheme='exponential-fit'; g=mvfp_grid_cache(sol);
f1=readtable(fullfile(out,'fig1_frozen_summary.csv'));
assert(all(abs(f1.initialMass-1)<1e-12));
assert(min(f1.P0)==0 && max(f1.P0)==1 && any(isfinite(f1.xSingleMinimum)));
both=isfinite(f1.xSaddle);
assert(all(f1.xLeftMinimum(both)<f1.xSaddle(both)) && all(f1.xSaddle(both)<f1.xRightMinimum(both)));
assert(all(f1.barrierLeftToRight(both)>=0) && all(f1.barrierRightToLeft(both)>=0));
cross=find((f1.log10TauRightToSaddle(1:end-1)-2).*(f1.log10TauRightToSaddle(2:end)-2)<=0);
assert(numel(cross)==1,'Expected one frozen tau_RL=T crossing on sampled range.');
bracket=[f1.P0(cross),f1.P0(cross+1)];
frozenCrossing=fzero(@(q)logtau(q,20001,p,sol,g)-2,bracket,optimset('TolX',1e-10));
rows=struct([]);
for q=[.06 .08 frozenCrossing]
    [l1,f2001]=logtau(q,20001,p,sol,g); [l2,~]=logtau(q,40001,p,sol,g);
    sf=sol; sf.xGridN=4001; gf=mvfp_grid_cache(sf); [~,f4001]=logtau(q,20001,p,sf,gf);
    rows=[rows,struct('P0',q,'log10TauRL_N20001',l1,'log10TauRL_N40001',l2, ...
        'absLog10Difference',abs(l1-l2),'F0_Nx2001',f2001,'F0_Nx4001',f4001, ...
        'absF0Difference',abs(f2001-f4001))]; %#ok<AGROW>
end
assert(max([rows.absLog10Difference])<1e-5 && max([rows.absF0Difference])<1e-10);
writetable(struct2table(rows),fullfile(out,'fig1_quadrature_validation.csv'));
f2=readtable(fullfile(out,'fig2_dynamic_full_traces.csv'));
assert(numel(unique(f2.P0))==8 && all(f2.Nx==2001));
assert(all(f2.PFixed>=0 & f2.PFixed<=1) && max(f2.rawMassError)<1e-6);
e=readtable(fullfile(out,'fig2_event_timing.csv'));
e08=e(abs(e.P0-.08)<1e-12 & e.q==.5,:); e0647=e(abs(e.P0-.0647)<1e-12 & e.q==.5,:);
assert(isfinite(e08.time) && isnan(e0647.time));
h=jsondecode(fileread(fullfile(out,'fig3_readout_heatmap.json')));
assert(min(h.P0)==0 && max(h.P0)==.3 && all(diff(h.P0)>1e-12));
assert(size(h.PFixed,1)==numel(h.P0) && size(h.PFixed,2)==numel(h.t));
assert(all(isfinite(h.PFixed(:))));
d=jsondecode(fileread(fullfile(out,'fig3_density_evolution.json')));
assert(max(d.t)==200 && size(d.rhoRaw,1)==numel(d.x) && size(d.rhoRaw,2)==numel(d.t));
c=readtable(fullfile(out,'fig3_t50_censoring_T100.csv'));
censored=~isfinite(c.t50Actual); assert(all(c.plottedTime(censored)==100));
assert(all(c.t50Actual(~censored)<=100));
v=readtable(fullfile(out,'fig2_fine_grid_validation.csv'));
assert(max(v.absEndpointDifference)<.01);
valid=isfinite(v.t50_Nx2001) & isfinite(v.t50_Nx4001);
assert(max(abs(v.t50_Nx2001(valid)-v.t50_Nx4001(valid)))<1);
% Reuse both already computed grids to compare all requested arrival levels
% and the fold sampling intervals, not only the half-occupancy event.
cache=fullfile(root,'outputs','pairwise_paper_revision_v2');
co=load(fullfile(cache,'fitted_flux_controls','controls.mat'),'result');
fi=load(fullfile(cache,'fourth_batch','figures123','fine_controls_nx4001','controls.mat'),'result');
evrows=struct([]);
for q0=[.06 .08]
    a=get_case(co.result,q0); b=get_case(fi.result,q0);
    ba=fold_bracket(a.full); bb=fold_bracket(b.full);
    for level=[.2 .4 .5 .6 .8]
        ea=mvfp_revision_first_arrival(a.full.t,a.full.PFixed,level);
        eb=mvfp_revision_first_arrival(b.full.t,b.full.PFixed,level);
        delta=nan;
        if isfinite(ea.time)&&isfinite(eb.time),delta=abs(ea.time-eb.time);assert(delta<1);end
        assert(isfinite(ea.time)==isfinite(eb.time),'Grid event status mismatch.');
        evrows=[evrows,struct('P0',q0,'q',level,'timeNx2001',ea.time,'timeNx4001',eb.time, ...
            'absTimeDifference',delta,'statusNx2001',string(ea.status),'statusNx4001',string(eb.status), ...
            'foldLowerNx2001',ba(1),'foldUpperNx2001',ba(2), ...
            'foldLowerNx4001',bb(1),'foldUpperNx4001',bb(2))]; %#ok<AGROW>
    end
end
writetable(struct2table(evrows),fullfile(out,'fig2_fine_grid_event_validation.csv'));
validation=struct('passed',true,'createdAt',char(datetime('now')), ...
    'frozenTauRLEquals100P0',frozenCrossing,'frozenTauCrossingSearchBracket',bracket, ...
    'maxQuadratureLog10Difference',max([rows.absLog10Difference]), ...
    'maxInitialFieldGridDifference',max([rows.absF0Difference]), ...
    'P008_t50',e08.time,'P00647_rightCensored',isnan(e0647.time), ...
    'figure3UniqueP0Count',numel(h.P0),'maxFigure2MassError',max(f2.rawMassError), ...
    'limitations','Local grid and quadrature checks; not global convergence certification. Event times interpolate saved .1 or .5 output samples.');
write_json(fullfile(out,'validation.json'),validation);
meta=jsondecode(fileread(fullfile(out,'fig1_metadata.json')));
meta.frozenTauRLEquals100P0=frozenCrossing; meta.quadratureValidationFile='fig1_quadrature_validation.csv';
write_json(fullfile(out,'fig1_metadata.json'),meta);
disp(validation);
end

function c=get_case(r,p)
for k=1:numel(r.cases)
    c=r.cases{k}.diagnostic.cases(1); if abs(c.P0-p)<1e-12,return;end
end
error('Missing control case.');
end

function b=fold_bracket(s)
i=find(isfinite(s.barrierRL),1,'last');b=[nan nan];
if ~isempty(i),b(1)=s.t(i);if i<numel(s.t),b(2)=s.t(i+1);end;end
end

function [val,F]=logtau(q,n,p,s,g)
rho=initial_density(g.x,p,q,s.initStd,s); F=p.sigma1*trapz(g.x,g.x.*rho)+p.bias;
[~,~,~,lr]=compute_mvfp_mfpt_point_from_F(F,p,g.x,-1,1,n); val=lr/log(10);
end

function write_json(path,value)
f=fopen(path,'w','n','UTF-8'); assert(f>0); c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s',jsonencode(value,'ConvertInfAndNaN',true));
end
