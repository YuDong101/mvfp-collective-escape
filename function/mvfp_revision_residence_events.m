function event = mvfp_revision_residence_events(t,distance,exitRadii,entryRadius)
%MVFP_REVISION_RESIDENCE_EVENTS First entry and subsequent outward exits.
t=t(:); distance=distance(:); exitRadii=exitRadii(:).';
if numel(t)~=numel(distance) || numel(t)<1 || any(diff(t)<=0) || ...
        any(~isfinite(t)) || any(~isfinite(distance))
    error('MVFP:RevisionResidence:InvalidEventTrace', ...
        'Time and distance must be finite vectors of equal length with increasing time.');
end
if ~isscalar(entryRadius) || ~isfinite(entryRadius) || entryRadius<=0 || isempty(exitRadii) || any(~isfinite(exitRadii)) || any(exitRadii<=entryRadius)
    error('MVFP:RevisionResidence:InvalidEventRadii','Exit radii must exceed the entry radius.');
end

inIndex=find(distance<entryRadius,1);
event=struct('entryRadius',entryRadius,'exitRadii',exitRadii, ...
    'tIn',nan,'tInBracket',[nan nan],'tOut',nan(size(exitRadii)), ...
    'tOutBracket',nan(numel(exitRadii),2),'tRes',nan(size(exitRadii)), ...
    'status',strings(size(exitRadii)),'outward',false(size(exitRadii)), ...
    'censorTime',t(end),'eventTimeUncertainty',nan(size(exitRadii)));
if isempty(inIndex)
    event.status(:)="never-entered";
    return
end
if inIndex==1
    event.tIn=t(1); event.tInBracket=[t(1) t(1)];
else
    event.tIn=linear_crossing(t(inIndex-1:inIndex),distance(inIndex-1:inIndex),entryRadius);
    event.tInBracket=t(inIndex-1:inIndex).';
end
for k=1:numel(exitRadii)
    r=exitRadii(k);
    j=find((1:numel(t)-1)'>=inIndex & distance(1:end-1)<r & distance(2:end)>=r,1);
    if isempty(j)
        event.status(k)="entered-censored";
        continue
    end
    event.tOut(k)=linear_crossing(t(j:j+1),distance(j:j+1),r);
    event.tOutBracket(k,:)=t(j:j+1).';
    event.tRes(k)=event.tOut(k)-event.tIn;
    event.status(k)="exited"; event.outward(k)=distance(j+1)>distance(j);
    event.eventTimeUncertainty(k)=max(abs(event.tInBracket-event.tIn)) + ...
        max(abs(event.tOutBracket(k,:)-event.tOut(k)));
end
end

function tc=linear_crossing(t,y,level)
if y(2)==y(1), tc=t(2); else, tc=t(1)+(level-y(1))*(t(2)-t(1))/(y(2)-y(1)); end
end
