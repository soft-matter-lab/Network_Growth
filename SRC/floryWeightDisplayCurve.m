function [x,y] = floryWeightDisplayCurve(fit,edges,p)
% Smooth analytic weight curve. Optional p draws a prediction from number fit.
if nargin<3, p=fit.p; end
x=[]; y=[];
if ~isfinite(p), return; end
width=edges(2)-edges(1); c=fit.supportMin;
validateattributes(width,{'numeric'},{'scalar','integer','positive'});
if any(abs(diff(edges)-width)>1e-9) || any(abs(edges+.5-round(edges+.5))>1e-9)
    error('NetworkGrowth:FloryCurveBins','Weight display requires equal integer-width bins with half-integer edges.');
end
centers=(edges(1:end-1)+edges(2:end))/2;
if p==0
    [~,~,bin]=histcounts(c,edges);
    if bin>0, x=centers(bin); y=1; end
    return;
end
valid=edges(2:end)>c;
if ~any(valid), return; end
xmin=centers(find(valid,1));
x=unique([linspace(xmin,edges(end),1000),centers(valid)]);
if width==1
    y=floryWeightModel(x,p,c)';
else
    % Extend the finite geometric first-moment sum continuously between bins.
    lo=max(c,x-(width-1)/2); hi=x+(width-1)/2; len=hi-lo+1;
    q=1-p; logp=log(p); rho=exp(len*logp);
    y=exp((lo-c)*logp)./(c*q+p).* ...
        ((lo*q+p).*(-expm1(len*logp))-len*q.*rho);
    y=max(0,y);
    % Exact integer-bin sums also avoid cancellation in the p->1 limit.
    exact=floryWeightBinProbabilities(p,c,edges);
    for k=find(valid), y(x==centers(k))=exact(k); end
end
end
