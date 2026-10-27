function probability = floryWeightBinProbabilities(p,minCount,edges)
% Sum the analytic weight PMF at actual integer lengths in each display bin.
% Includes the final right edge, like histcounts. Never renormalize cropping.
validateattributes(edges,{'numeric'},{'vector','real','finite','increasing'});
validateattributes(p,{'numeric'},{'scalar','real','finite','>=',0,'<',1});
validateattributes(minCount,{'numeric'},{'scalar','integer','nonnegative','finite'});
if numel(edges)<2, error('NetworkGrowth:FloryWeightBins','At least two edges are required.'); end
edges=edges(:)';
lo=max(max(1,minCount),ceil(edges(1:end-1))); hi=ceil(edges(2:end))-1;
hi(end)=floor(edges(end)); probability=zeros(size(lo));
for k=find(hi>=lo)
    probability(k)=sum(floryWeightModel(lo(k):hi(k),p,minCount));
end
end
