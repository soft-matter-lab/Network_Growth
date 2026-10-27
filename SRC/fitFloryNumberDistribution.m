function fit = fitFloryNumberDistribution(replicates, edges, options)
% Select a normalized model (zero probability below its support):
% shifted:  (1-p)*p^(n-c), n>=c; equal-replicate raw-count likelihood.
% original: (1-p)*p^(n-1), n>=1; absolute-PMF least squares on n>=max(1,c).
% c=minCount. For shifted it is a physical minimum, so shorter samples are
% rejected. For original it is only a fitting cutoff. No backward extension.
% Original fitting uses unit-count PMFs internally, independent of plot bins.
if nargin<3, options=struct(); end
options=mergeConfig(struct('model','shifted','minCount',0,'originalAnchorWeight',1, ...
    'bootstrapSamples',2000,'bootstrapSeed',1729),options);
model=validatestring(options.model,{'shifted','original'},mfilename,'model');
validateattributes(options.minCount,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(options.originalAnchorWeight,{'numeric'},{'scalar','real','>=',1,'finite'});
validateattributes(options.bootstrapSamples,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(options.bootstrapSeed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
validateattributes(edges,{'numeric'},{'vector','real','finite','increasing'});
if numel(edges)<2, error('NetworkGrowth:FloryBins','At least two bin edges are required.'); end
if ~iscell(replicates), error('NetworkGrowth:FlorySamples','Supply one sample vector per replicate in a cell array.'); end
edges=edges(:)'; c=options.minCount; support=c;
if strcmp(model,'original'), c=max(1,c); support=1; end
replicates=replicates(:); replicates=replicates(~cellfun(@isempty,replicates));
dataR=numel(replicates); retained=cell(dataR,1); counts=zeros(dataR,1);
for r=1:dataR
    v=replicates{r}(:);
    validateattributes(v,{'numeric'},{'vector','real','finite','integer','nonnegative'});
    if strcmp(model,'shifted') && any(v<c)
        error('NetworkGrowth:FloryBelowMinimum', ...
            'Shifted model minimum is %g, but replicate %d contains n=%g. Lower floryMinCount; no short chains are silently discarded.',c,r,min(v));
    end
    replicates{r}=v; retained{r}=v(v>=c); counts(r)=numel(retained{r});
end
included=counts>0; retained=retained(included);
R=numel(retained); means=cellfun(@mean,retained);
fit=struct('p',NaN,'pCI95',[NaN NaN],'meanCount',NaN,'minCount',c, ...
    'modelStyle',model,'supportMin',support,'replicateCount',R,'chainCount',sum(counts), ...
    'replicateMeans',means,'dataReplicateCount',dataR,'totalChainCount',sum(cellfun(@numel,replicates)), ...
    'retainedCounts',counts,'includedReplicates',find(included),'replicateP',nan(R,1), ...
    'binProbability',nan(1,numel(edges)-1),'conditionalBinProbability',nan(1,numel(edges)-1), ...
    'probabilityBelowCutoff',NaN,'displayBinWeight',nan(1,numel(edges)-1), ...
    'observedTailFraction',NaN,'cdfDistance',NaN,'binRMSE',NaN,'outsideBinProbability',NaN, ...
    'fitCdfDistance',NaN,'fitBinRMSE',NaN,'bootstrapP',[], ...
    'bootstrapSamples',options.bootstrapSamples,'bootstrapSeed',options.bootstrapSeed, ...
    'originalAnchorWeight',NaN,'originalDiagnostics',struct(), ...
    'method','equal-replicate raw-count log likelihood', ...
    'model',sprintf('P(n)=(1-p)*p^(n-%g), n>=%g; zero below support',support,support), ...
    'ciMethod','95% percentile bootstrap of whole replicates; no within-run resampling', ...
    'status','no chains at or above fitting cutoff');
if R==0, return; end
fit.meanCount=mean(means); fit.status='fitted';
if strcmp(model,'shifted')
    excess=fit.meanCount-c;
    fit.p=excess/(1+excess);
    fit.replicateP=(means-c)./(1+means-c);
    if R>=2 && options.bootstrapSamples>0
        stream=RandStream('mt19937ar','Seed',options.bootstrapSeed);
        indices=randi(stream,R,R,options.bootstrapSamples);
        bootMean=mean(reshape(means(indices),size(indices)),1)-c;
        fit.bootstrapP=bootMean./(1+bootMean);
        sorted=sort(fit.bootstrapP);
        position=1+(numel(sorted)-1)*[.025 .975];
        low=floor(position); high=ceil(position); alpha=position-low;
        fit.pCI95=sorted(low).*(1-alpha)+sorted(high).*alpha;
    end
else
    allCounts=vertcat(replicates{:}); unitCounts=0:max(allCounts);
    unitEdges=-.5:1:(unitCounts(end)+.5);
    rows=zeros(R,numel(unitCounts)); indices=find(included);
    for r=1:R
        v=replicates{indices(r)};
        rows(r,:)=histcounts(v,unitEdges)/numel(v); % normalize FULL replicate
    end
    original=fitOriginalFloryNumberDistribution(rows,unitCounts, ...
        struct('minCount',c,'anchorCount',c,'anchorWeight',options.originalAnchorWeight, ...
        'bootstrapSamples',options.bootstrapSamples,'bootstrapSeed',options.bootstrapSeed));
    fit.p=original.p; fit.pCI95=original.pCI95; fit.bootstrapP=original.bootstrapP;
    fit.originalAnchorWeight=options.originalAnchorWeight;
    fit.originalDiagnostics=original; fit.method=original.method; fit.ciMethod=original.ciMethod;
end
fit.binProbability=geometricBins(edges,fit.p,support);
fit.conditionalBinProbability=geometricBins(edges,fit.p,c);
fit.displayBinWeight=fit.binProbability;
fit.observedTailFraction=mean(counts./cellfun(@numel,replicates));
fit.probabilityBelowCutoff=geometricCDF(c-1,fit.p,support);
fit.outsideBinProbability=max(0,1-sum(fit.binProbability));
% Full diagnostics describe the SELECTED normalized model and all data.
% Conditional diagnostics describe tail SHAPE, not the original LS objective.
[fit.cdfDistance,fit.binRMSE]=discrepancy(replicates,edges,fit.p,support,fit.binProbability);
[fit.fitCdfDistance,fit.fitBinRMSE]=discrepancy(retained,edges,fit.p,c,fit.conditionalBinProbability);
end

function prob=geometricBins(edges,p,n0)
% Match histcounts: [left,right), except the final right edge is included.
lo=max(n0,ceil(edges(1:end-1))); hi=ceil(edges(2:end))-1;
hi(end)=floor(edges(end));
prob=zeros(size(lo)); included=hi>=lo;
if p==0
    prob(included & lo<=n0 & hi>=n0)=1;
else
    logp=log(p);
    prob(included)=exp((lo(included)-n0)*logp).* ...
        (-expm1((hi(included)-lo(included)+1)*logp));
end
end

function [distance,rmse]=discrepancy(replicates,edges,p,n0,prob)
R=numel(replicates);
empirical=zeros(size(prob)); values=cell(R,1); weights=cell(R,1);
for r=1:R
    count=numel(replicates{r});
    empirical=empirical+histcounts(replicates{r},edges)/(R*count);
    values{r}=replicates{r}; weights{r}=ones(count,1)/(R*count);
end
rmse=sqrt(mean((empirical-prob).^2));
[k,~,group]=unique(vertcat(values{:}));
mass=accumarray(group,vertcat(weights{:})); cumulative=cumsum(mass);
distance=max([abs(cumulative-geometricCDF(k,p,n0)); ...
    abs(cumulative-mass-geometricCDF(k-1,p,n0))]);
end

function value=geometricCDF(k,p,n0)
value=zeros(size(k)); valid=k>=n0;
if p==0
    value(valid)=1;
else
    value(valid)=-expm1((floor(k(valid))-n0+1)*log(p));
end
end
