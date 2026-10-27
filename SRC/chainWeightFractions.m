function [fraction, totalSegments] = chainWeightFractions(kuhnCount, edges)
% Fraction of all completed-chain Kuhn segments carried by each length bin.
% Sum actual n values within a bin (not chain count times the bin center).
% Out-of-range chains remain in the denominator; an empty/zero total is undefined.
kuhnCount = kuhnCount(:);
validateattributes(kuhnCount,{'numeric'},{'real','finite','nonnegative'});
validateattributes(edges,{'numeric'},{'vector','real','finite','increasing'});
if numel(edges)<2, error('NetworkGrowth:WeightBins','At least two edges are required.'); end
totalSegments = sum(kuhnCount);
fraction = zeros(1,numel(edges)-1);
if totalSegments==0, return; end
[~,~,bin] = histcounts(kuhnCount,edges);
inside = bin>0;
fraction = accumarray(bin(inside),kuhnCount(inside),[numel(edges)-1 1],@sum,0)'/totalSegments;
end
