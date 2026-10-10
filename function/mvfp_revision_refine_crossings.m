function [x,y,reliable,byReadout] = mvfp_revision_refine_crossings( ...
    x,y,q,reliable,evaluatePoint,maxWidth,maxPasses)
% Refine every sampled crossing while sharing new samples across readouts.
x=x(:).'; reliable=logical(reliable(:).');
if size(y,2)~=numel(x) || numel(reliable)~=numel(x)
    error('MVFP:Revision:InvalidRefinementInput','Samples and reliability must match x.');
end
for pass=1:maxPasses
    byReadout=find_all();
    mids=[];
    for j=1:numel(byReadout)
        b=byReadout{j};
        active=strcmp({b.status},'bracketed') & ([b.upper]-[b.lower]>maxWidth);
        mids=[mids,.5*([b(active).lower]+[b(active).upper])]; %#ok<AGROW>
    end
    mids=setdiff(unique(mids),x);
    if isempty(mids), break; end
    for mid=mids
        [values,isReliable]=evaluatePoint(mid);
        values=values(:);
        if numel(values)~=size(y,1) || ~isscalar(isReliable)
            error('MVFP:Revision:InvalidRefinementEvaluation', ...
                'Evaluator must return one value per readout and scalar reliability.');
        end
        x(end+1)=mid; y(:,end+1)=values; reliable(end+1)=logical(isReliable); %#ok<AGROW>
    end
    [x,order]=sort(x); y=y(:,order); reliable=reliable(order);
end
byReadout=find_all();
for j=1:numel(byReadout)
    for k=1:numel(byReadout{j})
        if strcmp(byReadout{j}(k).status,'bracketed') && ...
                byReadout{j}(k).upper-byReadout{j}(k).lower>maxWidth
            byReadout{j}(k).status='refinement-limit';
        end
    end
end

    function allBrackets=find_all()
        allBrackets=cell(1,size(y,1));
        for jj=1:size(y,1)
            allBrackets{jj}=mvfp_revision_crossing_intervals(x,y(jj,:),q,reliable);
        end
    end
end
