function [probability,survival] = floryWeightModel(n,p,minCount)
% Weight fractions obtained by weighting a shifted geometric NUMBER model by
% actual n, not by excess n-minCount: W(n)=n*P(n)/E[n].
% W(n)=n*(1-p)^2*p^(n-c)/(c*(1-p)+p), n>=c.
% c=1 is original Flory. c=0 has the same weight law as c=1 for p>0;
% at c=0,p=0 we return its p->0 limit (a unit mass at n=1).
% Rows index n, columns index p. Real n may be used to draw a smooth curve.
% survival(n)=Pr(N>=ceil(n)) is for the WEIGHT distribution.
validateattributes(n,{'numeric'},{'vector','real','finite'});
validateattributes(p,{'numeric'},{'vector','real','finite','>=',0,'<',1});
validateattributes(minCount,{'numeric'},{'scalar','integer','nonnegative','finite'});
n=n(:); p=p(:)'; c=max(1,minCount); q=1-p; denominator=c*q+p;
probability=zeros(numel(n),numel(p)); valid=n>=c;
probability(valid,:)=n(valid).*q.^2./denominator.*p.^(n(valid)-c);
if nargout>1
    k=max(c,ceil(n));
    survival=p.^(k-c).*(k.*q+p)./denominator;
end
end
