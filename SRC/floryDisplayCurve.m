function [x,y] = floryDisplayCurve(fit,edges)
% Smooth curve of the SELECTED normalized model; no backward extrapolation.
% Unit bins: (1-p)*p^(x-supportMin) on its support only.
% Grouped bins: a continuous envelope through exact integer-bin masses.
% A boundary bin may straddle supportMin; its center represents the whole bin.
x=[]; y=[];
if ~isfinite(fit.p), return; end
width=edges(2)-edges(1);
validateattributes(width,{'numeric'},{'scalar','integer','positive'});
if any(abs(diff(edges)-width)>1e-9) || any(abs(edges+.5-round(edges+.5))>1e-9)
    error('NetworkGrowth:FloryCurveBins','Flory display requires equally spaced integer-width bins with half-integer edges.');
end
% Legacy saved fits used minCount as the conditional tail's support.
if isfield(fit,'supportMin'), support=fit.supportMin; else, support=fit.minCount; end
centers=(edges(1:end-1)+edges(2:end))/2;
if fit.p==0
    [~,~,bin]=histcounts(support,edges);
    if bin>0, x=centers(bin); y=1; end
    return;
end
valid=edges(2:end)>support;
if ~any(valid), return; end
xmin=centers(find(valid,1)); xmax=edges(end);
x=unique([linspace(xmin,xmax,1000),centers(valid)]);
lo=max(support,x-(width-1)/2); hi=x+(width-1)/2;
logp=log(fit.p);
y=exp((lo-support)*logp).*(-expm1((hi-lo+1)*logp));
end
