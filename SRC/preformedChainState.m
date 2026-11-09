function state = preformedChainState(chainCounts,nodes,bonds,activeEnds,cfg)
% Reservoir and incorporated-chain accounting; loops are completed chains.
feed = cfg.preformedFeed;
lengths = feed.chainLengths;
completed = arrayfun(@(n)sum(bonds(:,3)==n),lengths);
dangling = arrayfun(@(n)sum(activeEnds(:,4)==n),lengths);
if ~isequal(chainCounts+completed+dangling,feed.initialChainCounts) || ...
        any(~ismember(bonds(:,3),lengths)) || ...
        any(~ismember(activeEnds(activeEnds(:,4)>0,4),lengths))
    error('NetworkGrowth:ChainConservation','Preformed chain identities/counts were not conserved.');
end
state = struct('chainLengths',lengths,'initialChainCounts',feed.initialChainCounts, ...
    'freeChainCounts',chainCounts,'danglingChainCounts',dangling, ...
    'completedChainCounts',completed,'freeKuhnSegments',sum(chainCounts.*lengths), ...
    'incorporatedKuhnSegments',sum((completed+dangling).*lengths), ...
    'formedCrosslinkers',size(nodes,1),'unusedKuhnBudget',feed.unusedKuhnBudget);
denom = feed.initialChainCounts; denom(denom==0) = NaN;
state.capturedPercentBySpecies = 100*(completed+dangling)./denom;
state.completedPercentBySpecies = 100*completed./denom;
state.capturedPercent = 100*sum(completed+dangling)/sum(feed.initialChainCounts);
state.completedPercent = 100*sum(completed)/sum(feed.initialChainCounts);
end
