function event = mvfp_revision_first_arrival(t,p,q)
% Fixed-threshold occupancy event; a downward recrossing is not a first hit.
t=t(:); p=p(:);
if numel(t)~=numel(p) || isempty(t) || any(~isfinite(t)) || ...
        any(diff(t)<=0) || any(~isfinite(p)) || ~isscalar(q) || ~isfinite(q)
    error('MVFP:Revision:InvalidEventInput','Use matching finite samples on increasing times.');
end
idx=find(p>=q,1,'first');
event=struct('level',q,'time',nan,'status','right-censored', ...
    'bracket',[t(end),inf],'observationEnd',t(end), ...
    'method','first-upward-hit-linear-output-interpolation');
if isempty(idx), return; end
event.status='reached';
if idx==1
    event.time=t(1); event.bracket=[t(1),t(1)];
else
    event.bracket=[t(idx-1),t(idx)];
    event.time=t(idx-1)+(q-p(idx-1))*(t(idx)-t(idx-1))/(p(idx)-p(idx-1));
end
end
