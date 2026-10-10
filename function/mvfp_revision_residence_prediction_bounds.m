function bounds=mvfp_revision_residence_prediction_bounds(pcal,phold,tcal,eventBounds,pb,rate)
% Shared systematic inputs occur once in the entire calibrated prediction.
% Bound g(Pb) using its derivative and a midpoint Lipschitz enclosure.
pcal=pcal(:); tcal=tcal(:); eventBounds=eventBounds(:);
if numel(pcal)~=numel(tcal) || numel(pcal)~=numel(eventBounds) || isempty(pcal) || ...
        numel(pb)~=2 || numel(rate)~=2 || any(~isfinite([pcal;phold;tcal;eventBounds;pb(:);rate(:)])) || ...
        pb(2)<pb(1) || rate(1)<=0 || rate(2)<rate(1) || any(eventBounds<0)
    error('MVFP:RevisionResidence:InvalidPredictionInputs','Invalid prediction intervals.');
end
pc=mean(pb); half=diff(pb)/2;
dc=abs(pcal-pc)-half; dh=abs(phold-pc)-half;
if any(dc<=0) || dh<=0, bounds=[-Inf Inf]; return; end
g=mean(log(abs(pcal-pc)))-log(abs(phold-pc));
% Conservative derivative absolute bound; no sqrt(n) reduction of Pb.
variation=half*(mean(1./dc)+1/dh);
gi=g+[-variation variation]; ratios=[gi(1)/rate(1) gi(1)/rate(2) gi(2)/rate(1) gi(2)/rate(2)];
ev=mean(eventBounds); bounds=mean(tcal)+[min(ratios)-ev max(ratios)+ev];
end
