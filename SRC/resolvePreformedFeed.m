function feed = resolvePreformedFeed(cfg)
% xA retains A/(A + Kuhn segments); chains are never treated as single segments.
cfg = mergeConfig(defaultNetworkBConfig(),cfg);
validateattributes(cfg.xA,{'numeric'},{'scalar','real','finite','>',0,'<',1});
validateattributes(cfg.fA,{'numeric'},{'scalar','integer','finite','>=',2});
lengths = cfg.chainLengths(:)';
validateattributes(lengths,{'numeric'},{'vector','nonempty','real','finite','integer','positive'});
if numel(lengths)>3 || numel(unique(lengths))~=numel(lengths)
    error('NetworkGrowth:ChainSpecies','Supply one to three distinct chain lengths.');
end
total = resolveNtotal(cfg);
validateattributes(total,{'numeric'},{'scalar','real','finite','integer','positive'});
initialA = round(cfg.xA*total);
if initialA<1, error('NetworkGrowth:NoSeed','The feed must contain a seed crosslinker.'); end
budget = total-initialA;
fractions = cfg.chainNumberFractions(:)';
if isempty(cfg.chainCounts)
    validateattributes(fractions,{'numeric'},{'vector','real','finite','nonnegative','numel',numel(lengths)});
    if abs(sum(fractions)-1)>1e-10
        error('NetworkGrowth:ChainFractions','Chain-number fractions must sum to one.');
    end
    ideal = budget/sum(fractions.*lengths)*fractions;
    counts = floor(ideal);
    remainder = budget-sum(counts.*lengths);
    % Largest fractional remainders; add at most one per species, never
    % exceed the budget or silently replace long chains with short ones.
    [~,order] = sort(ideal-counts,'descend');
    for s = order
        if ideal(s)>counts(s) && lengths(s)<=remainder
            counts(s) = counts(s)+1; remainder = remainder-lengths(s);
        end
    end
    allocation = 'Kuhn budget with largest-remainder rounding';
else
    counts = cfg.chainCounts(:)';
    validateattributes(counts,{'numeric'},{'vector','real','finite','integer','nonnegative','numel',numel(lengths)});
    if sum(counts.*lengths)>budget
        error('NetworkGrowth:ChainBudget','Explicit chainCounts exceed the available Kuhn budget (%d).',budget);
    end
    fractions = counts/max(sum(counts),1);
    allocation = 'explicit chain counts';
end
if sum(counts)==0
    error('NetworkGrowth:NoPreformedChains','The feed contains no whole preformed chains; increase the budget or reduce chain lengths.');
end
feed = struct('chainLengths',lengths,'requestedNumberFractions',fractions, ...
    'initialChainCounts',counts,'actualNumberFractions',counts/sum(counts), ...
    'initialCrosslinkers',initialA,'kuhnBudget',budget, ...
    'initialKuhnCount',sum(counts.*lengths), ...
    'unusedKuhnBudget',budget-sum(counts.*lengths), ...
    'chainEndArmRatio',2*sum(counts)/(cfg.fA*initialA),'allocation',allocation);
end
