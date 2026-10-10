function out=mvfp_revision_outgoing_modal_window(t,d,a,stableFraction,residualSeed,cfg)
% Select by geometry and mode resolution, independently of growth-fit errors.
t=t(:); d=d(:); a=a(:); stableFraction=stableFraction(:);
assert(numel(t)==numel(d) && numel(a)==numel(t) && numel(stableFraction)==numel(t));
out=struct('indices',[],'times',[],'reason','no-entry','minimumIndex',NaN);
first=find(d<cfg.entryRadius,1); if isempty(first), return; end
last=find((1:numel(t))'>first & d>=cfg.entryRadius,1);
if isempty(last), last=numel(t)+1; end
minimum=[];
for j=max(2,first):min(last-2,numel(t)-1)
    if d(j)<=d(j-1) && d(j+1)>=d(j) && (d(j)<d(j-1) || d(j+1)>d(j))
        minimum=j; break;
    end
end
out.reason='no-resolved-local-minimum'; if isempty(minimum), return; end
out.minimumIndex=minimum;
good=isfinite(a) & abs(a)>0 & isfinite(stableFraction) & ...
    stableFraction<=cfg.stableFractionMax & ...
    abs(residualSeed)./max(abs(a),realmin)<=cfg.residualSeedFractionMax;
good(1:minimum-1)=false; good(last:end)=false;
out.reason='no-resolved-outgoing-modal-window';
starts=find(good & [true;~good(1:end-1)]);
ends=find(good & [~good(2:end);true]);
for k=1:numel(starts)
    ids=starts(k):ends(k);
    if numel(ids)>=cfg.minLocalSamples && t(ids(end))-t(ids(1))>=cfg.minLocalDuration
        out.indices=ids; out.times=t(ids); out.reason='selected-by-geometry-and-resolution'; return;
    end
end
end
