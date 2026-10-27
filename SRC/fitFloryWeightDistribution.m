function fit = fitFloryWeightDistribution(replicates,edges,options)
% Independently fit equal-replicate mean Kuhn WEIGHT fractions by least squares.
% The model is W(n)=n*P(n)/E[n], where P is the selected number-fraction model.
% Use full per-replicate Kuhn totals, no tail renormalization or peak anchor.
% Fit unit-count probabilities regardless of plot bins and include the infinite
% unobserved tail analytically. Bootstrap whole replicates, not Kuhn segments.
if nargin<3, options=struct(); end
options=mergeConfig(struct('model','shifted','minCount',0,'bootstrapSamples',2000, ...
    'bootstrapSeed',1729,'numberFit',[]),options);
model=validatestring(options.model,{'shifted','original'},mfilename,'model');
validateattributes(options.minCount,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(options.bootstrapSamples,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(options.bootstrapSeed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
validateattributes(edges,{'numeric'},{'vector','real','finite','increasing'});
if numel(edges)<2, error('NetworkGrowth:FloryWeightBins','At least two edges are required.'); end
if ~iscell(replicates), error('NetworkGrowth:FlorySamples','Supply one raw-count vector per replicate.'); end
edges=edges(:)'; numberMin=options.minCount;
if strcmp(model,'original'), numberMin=1; end
c=max(1,numberMin); cutoff=max(c,options.minCount);
replicates=replicates(:); replicates=replicates(~cellfun(@isempty,replicates));
dataR=numel(replicates); totals=zeros(dataR,1); tailCounts=zeros(dataR,1); maxN=0;
for r=1:dataR
    v=replicates{r}(:);
    validateattributes(v,{'numeric'},{'vector','real','integer','finite','nonnegative'});
    if strcmp(model,'shifted') && any(v<numberMin)
        error('NetworkGrowth:FloryBelowMinimum','Shifted minimum %g exceeds an observed count in replicate %d.',numberMin,r);
    end
    replicates{r}=v; totals(r)=sum(v); tailCounts(r)=sum(v>=cutoff); maxN=max(maxN,max(v));
end
positive=totals>0; included=positive & tailCounts>0; R=sum(included);
fit=struct('p',NaN,'pCI95',[NaN NaN],'modelStyle',model,'supportMin',c, ...
    'numberSupportMin',numberMin,'minCount',cutoff,'replicateCount',R, ...
    'dataReplicateCount',dataR,'dataWeightReplicateCount',sum(positive), ...
    'includedReplicates',find(included),'chainCount',sum(tailCounts(included)), ...
    'totalKuhnCount',sum(totals(included)),'binProbability',nan(1,numel(edges)-1), ...
    'outsideBinProbability',NaN,'probabilityBelowCutoff',NaN,'observedTailFraction',NaN, ...
    'cdfDistance',NaN,'binRMSE',NaN,'tailSSE',NaN, ...
    'bootstrapP',[],'bootstrapSamples',options.bootstrapSamples,'bootstrapSeed',options.bootstrapSeed, ...
    'pFromNumber',NaN,'pDifference',NaN,'pDifferenceCI95',[NaN NaN], ...
    'numberPredictionBinProbability',nan(1,numel(edges)-1),'numberPredictionTailSSE',NaN, ...
    'method','least squares of equal-replicate mean full weight fractions, with infinite tail', ...
    'model','W(n)=n*(1-p)^2*p^(n-c)/(c*(1-p)+p), n>=c', ...
    'ciMethod','whole-replicate percentile bootstrap; paired draws for difference from number fit', ...
    'status','no positive Kuhn-weight population at or above cutoff');
if R==0, return; end
unitCounts=0:maxN; unitEdges=-.5:1:(maxN+.5);
allRows=zeros(dataR,numel(unitCounts));
for r=find(positive)'
    allRows(r,:)=chainWeightFractions(replicates{r},unitEdges);
end
rows=allRows(included,:); fittingMean=mean(rows,1); tail=unitCounts>=cutoff;
n=unitCounts(tail); target=fittingMean(tail);
[fit.p,fit.tailSSE]=solveWeightFit(target,n,c,cutoff);
fit.status='fitted';
if numberMin==0 && fit.p==0, fit.status='p->0 limit; number model at p=0 has zero Kuhn mass'; end
fit.binProbability=floryWeightBinProbabilities(fit.p,c,edges);
fit.outsideBinProbability=max(0,1-sum(fit.binProbability));
[~,tailMass]=floryWeightModel(cutoff,fit.p,c);
fit.probabilityBelowCutoff=1-tailMass; fit.observedTailFraction=sum(target);
[fit.cdfDistance,fit.binRMSE]=discrepancy(replicates(positive),edges,fit.p,c,fit.binProbability);
if R>=2 && options.bootstrapSamples>0
    stream=RandStream('mt19937ar','Seed',options.bootstrapSeed);
    indices=randi(stream,R,R,options.bootstrapSamples);
    % Five replicates have only 126 distinct bootstrap compositions. Reusing
    % identical draws preserves the exact bootstrap while avoiding repeated fits.
    [uniqueDraws,~,drawMap]=unique(sort(indices,1)','rows');
    uniqueP=zeros(size(uniqueDraws,1),1);
    for b=1:size(uniqueDraws,1)
        uniqueP(b)=solveWeightFit(mean(rows(uniqueDraws(b,:),tail),1),n,c,cutoff);
    end
    fit.bootstrapP=reshape(uniqueP(drawMap),1,[]);
    fit.pCI95=percentileCI(fit.bootstrapP);
end
number=options.numberFit;
if ~isempty(number) && isfinite(number.p) && strcmp(number.modelStyle,model) && number.supportMin==numberMin
    fit.pFromNumber=number.p; fit.pDifference=fit.p-number.p;
    fit.numberPredictionBinProbability=floryWeightBinProbabilities(number.p,c,edges);
    fit.numberPredictionTailSSE=weightSSE(target,n,number.p,c,cutoff);
    % A difference interval requires the SAME cohort and identical bootstrap draws.
    if isequal(fit.includedReplicates,number.includedReplicates) && ...
            fit.dataReplicateCount==number.dataReplicateCount && ...
            fit.bootstrapSeed==number.bootstrapSeed && ~isempty(fit.bootstrapP) && ...
            numel(fit.bootstrapP)==numel(number.bootstrapP)
        fit.pDifferenceCI95=percentileCI(fit.bootstrapP-number.bootstrapP);
    end
end
end

function [p,best]=solveWeightFit(data,n,c,cutoff)
grid=unique([linspace(0,1-1e-12,257),1-logspace(-12,-2,80)]);
fun=@(q)weightSSE(data,n,q,c,cutoff);
values=fun(grid); [best,k]=min(values); p=grid(k);
middle=values(2:end-1); left=values(1:end-2); right=values(3:end);
local=find(middle<=left & middle<=right & (middle<left | middle<right))+1;
opts=optimset('TolX',1e-11,'Display','off');
for j=local
    [candidate,value]=fminbnd(fun,grid(j-1),grid(j+1),opts);
    if value<best, p=candidate; best=value; end
end
end

function value=weightSSE(data,n,p,c,k)
q=1-p; denominator=c*q+p;
% Exact sum_{n=k}^infinity W(n)^2; factored to avoid 1-p^2 cancellation.
squaredTail=p.^(2*(k-c)).*q./denominator.^2.* ...
    (k^2*q.^2./(1+p)+2*k*p.^2.*q./(1+p).^2+p.^2.*(1+p.^2)./(1+p).^3);
value=max(0,sum(data.^2)+squaredTail-2*(data*floryWeightModel(n,p,c)));
end

function ci=percentileCI(values)
values=sort(values); position=1+(numel(values)-1)*[.025 .975];
lo=floor(position); hi=ceil(position); alpha=position-lo;
ci=values(lo).*(1-alpha)+values(hi).*alpha;
end

function [distance,rmse]=discrepancy(replicates,edges,p,c,prob)
R=numel(replicates); empirical=zeros(size(prob)); values=cell(R,1); weights=cell(R,1);
for r=1:R
    n=replicates{r}; empirical=empirical+chainWeightFractions(n,edges)/R;
    values{r}=n; weights{r}=n/(R*sum(n));
end
rmse=sqrt(mean((empirical-prob).^2));
[k,~,group]=unique(vertcat(values{:})); mass=accumarray(group,vertcat(weights{:}));
cumulative=cumsum(mass);
[~,survivalRight]=floryWeightModel(k+1,p,c);
[~,survivalLeft]=floryWeightModel(k,p,c);
distance=max([abs(cumulative-(1-survivalRight)); abs(cumulative-mass-(1-survivalLeft))]);
end
