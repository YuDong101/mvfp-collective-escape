function intervals = mvfp_revision_crossing_intervals(x,y,q,reliable)
% Keep every sampled crossing without joining across missing/unreliable data.
x=x(:); y=y(:);
if nargin<4, reliable=true(size(x)); end
reliable=logical(reliable(:));
if numel(x)~=numel(y) || numel(x)~=numel(reliable) || any(~isfinite(x)) || ...
        any(diff(x)<=0) || ~isscalar(q) || ~isfinite(q)
    error('MVFP:Revision:InvalidCrossingInput','Use increasing finite x and matched samples.');
end
valid=reliable & isfinite(y);
intervals=struct('lower',{},'upper',{},'yLower',{},'yUpper',{},'status',{});
for k=1:numel(x)
    if ~valid(k), continue; end
    if y(k)==q
        if k>1 && valid(k-1) && y(k-1)==q, continue; end
        j=k;
        while j<numel(x) && valid(j+1) && y(j+1)==q, j=j+1; end
        status='sampled-exact';
        if j>k, status='sampled-plateau'; end
        intervals(end+1)=struct('lower',x(k),'upper',x(j),'yLower',q,'yUpper',q,'status',status); %#ok<AGROW>
    elseif k<numel(x) && valid(k+1) && (y(k)-q)*(y(k+1)-q)<0
        intervals(end+1)=struct('lower',x(k),'upper',x(k+1), ...
            'yLower',y(k),'yUpper',y(k+1),'status','bracketed'); %#ok<AGROW>
    end
end
end
