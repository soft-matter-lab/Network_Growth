function fit = fitOriginalFloryNumberDistribution(probabilities, counts, options)
% Fit the ORIGINAL normalized Flory PMF (1-p)*p^(n-1), n=1,2,... .
% Rows are full, unit-count PMFs, normalized separately for each replicate.
% Fit their equal-replicate mean on n>=minCount without tail renormalization.
% Objective: sum_{n>=c}(data(n)-model(n))^2 +
%            (anchorWeight-1)*(data(a)-model(a))^2.
% Includes the infinite unobserved tail analytically, so display range cannot
% change p. anchorWeight=1 is ordinary least squares; larger is a SOFT anchor.
if nargin<3, options=struct(); end
options=mergeConfig(struct('minCount',4,'anchorCount',4,'anchorWeight',10, ...
    'bootstrapSamples',2000,'bootstrapSeed',1729),options);
validateattributes(options.minCount,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(options.anchorCount,{'numeric'},{'scalar','integer','>=',options.minCount,'finite'});
validateattributes(options.anchorWeight,{'numeric'},{'scalar','real','>=',1,'finite'});
validateattributes(options.bootstrapSamples,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(options.bootstrapSeed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
validateattributes(counts,{'numeric'},{'vector','nonempty','integer','nonnegative','finite','increasing'});
validateattributes(probabilities,{'numeric'},{'2d','nonempty','real','finite','nonnegative'});
counts=counts(:)';
if size(probabilities,2)~=numel(counts) || counts(1)>1 || any(diff(counts)~=1)
    error('NetworkGrowth:OriginalFloryBins','Supply consecutive unit-count bins starting at 0 or 1.');
end
if any(abs(sum(probabilities,2)-1)>1e-8)
    error('NetworkGrowth:OriginalFloryNormalization','Each row must sum to one over the FULL observed population.');
end
R=size(probabilities,1); empirical=mean(probabilities,1);
c=options.minCount; a=options.anchorCount; w=options.anchorWeight;
tail=counts>=c; n=counts(tail); data=empirical(tail);
if ~any(data>0)
    error('NetworkGrowth:EmptyOriginalFloryTail','There are no observations at or above minCount.');
end
anchor=zeros(R,1);
if any(counts==a), anchor=probabilities(:,counts==a); end
observed=mean(anchor);
[p,objective]=solveFit(data,n,observed,c,a,w);
[plainP,plainSSE]=solveFit(data,n,observed,c,a,1);
if a==1, peakP=0; peak=1; else, peakP=(a-1)/a; peak=(1-peakP)*peakP^(a-1); end
predicted=floryPMF(p,a);
fit=struct('p',p,'pCI95',[NaN NaN],'minCount',c,'anchorCount',a, ...
    'anchorWeight',w,'anchorObserved',observed,'anchorPredicted',predicted, ...
    'anchorResidual',predicted-observed,'anchorMaximum',peak,'pAtAnchorMaximum',peakP, ...
    'exactAnchorFeasible',observed<=peak+1e-12, ...
    'tailSSE',tailSSE(data,n,p,c),'weightedSSE',objective, ...
    'unanchoredP',plainP,'unanchoredTailSSE',plainSSE, ...
    'observedTailMass',sum(data),'modelTailMass',p^(c-1), ...
    'replicateCount',R,'counts',counts,'meanProbability',empirical, ...
    'modelProbability',floryPMF(p,counts),'bootstrapP',[], ...
    'bootstrapSamples',options.bootstrapSamples,'bootstrapSeed',options.bootstrapSeed, ...
    'model','P(n)=(1-p)*p^(n-1), n=1,2,...; zero mass at n=0', ...
    'method','equal-replicate mean PMF, infinite-tail least squares with a soft anchor', ...
    'ciMethod','95% percentile bootstrap of whole replicates; anchor re-estimated each draw');
if R>=2 && options.bootstrapSamples>0
    stream=RandStream('mt19937ar','Seed',options.bootstrapSeed);
    indices=randi(stream,R,R,options.bootstrapSamples);
    fit.bootstrapP=zeros(1,options.bootstrapSamples);
    for b=1:options.bootstrapSamples
        rows=indices(:,b);
        fit.bootstrapP(b)=solveFit(mean(probabilities(rows,tail),1),n,mean(anchor(rows)),c,a,w);
    end
    sorted=sort(fit.bootstrapP);
    position=1+(numel(sorted)-1)*[.025 .975];
    low=floor(position); high=ceil(position); alpha=position-low;
    fit.pCI95=sorted(low).*(1-alpha)+sorted(high).*alpha;
end
end

function [p,best]=solveFit(data,n,anchor,c,a,w)
% The anchor can produce two branches. Search the whole physical interval,
% refine EVERY sampled local minimum, and compare both boundaries.
grid=unique([linspace(0,1-1e-12,257),1-logspace(-12,-2,80),(a-1)/a]);
fun=@(q)tailSSE(data,n,q,c)+(w-1)*(floryPMF(q,a)-anchor).^2;
values=fun(grid);
[best,k]=min(values); p=grid(k);
local=find(values(2:end-1)<=values(1:end-2) & values(2:end-1)<=values(3:end))+1;
opts=optimset('TolX',1e-11,'Display','off');
for j=local
    [candidate,value]=fminbnd(fun,grid(j-1),grid(j+1),opts);
    if value<best, p=candidate; best=value; end
end
end

function value=tailSSE(data,n,p,c)
% Sum model(n)^2 from c to infinity is (1-p)/(1+p)*p^(2*(c-1)).
model=floryPMF(p,n);
value=sum(data.^2)+(1-p)./(1+p).*p.^(2*(c-1))-2*(data*model);
value=max(0,value); % roundoff at exact agreement
end

function value=floryPMF(p,n)
% Output rows index n; columns index p. n=0 has zero probability.
n=n(:); p=p(:)';
value=zeros(numel(n),numel(p)); valid=n>=1;
value(valid,:)=(1-p).*p.^(n(valid)-1);
end
