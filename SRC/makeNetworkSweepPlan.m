function [plan, sweep] = makeNetworkSweepPlan(baseCfg, sweep)
% Build a deterministic run plan without launching simulations.
% Replicate blocks reuse their seed across compositions; fA stays fixed.
if nargin < 2, sweep = struct(); end
baseCfg = resolveNetworkConfig(baseCfg);
defaults = struct('mode','replicated_composition','xA',[.01 .025 .05 .075 .10], ...
    'seeds',[1 2 3],'NtotalScaleByXA',[],'continueOnError',false,'saveResults',true, ...
    'parameterSets',{{}});
sweep = mergeConfig(defaults,sweep);
% Explicit parameter sets avoid an accidental Cartesian product of mixtures.
% Each set is crossed with the selected compositions and replicate seeds.
if ~isempty(sweep.parameterSets)
    if ~iscell(sweep.parameterSets) || ~all(cellfun(@(s)isstruct(s)&&isscalar(s),sweep.parameterSets))
        error('NetworkGrowth:SweepParameters','parameterSets must be a cell array of scalar configuration structs.');
    end
    plan = table(); inner = sweep; inner.parameterSets = {};
    for j = 1:numel(sweep.parameterSets)
        overrides = sweep.parameterSets{j};
        if any(~ismember(fieldnames(overrides),{'chainLengths','chainNumberFractions','chainCounts','fA','kAddition','kChem','kEndEnd','kSelfLoop'}))
            error('NetworkGrowth:SweepParameters','Parameter sets support chain mixtures, functionality, and reaction rates only.');
        end
        if ~strcmp(baseCfg.growthMode,'preformed') && any(startsWith(fieldnames(overrides),'chain'))
            error('NetworkGrowth:SweepParameters','Preformed-chain parameter sets require growthMode=preformed.');
        end
        part = makeNetworkSweepPlan(mergeConfig(baseCfg,overrides),inner);
        part.parameterSetIndex = repmat(j,height(part),1);
        plan = [plan;part]; %#ok<AGROW>
    end
    plan.runIndex = (1:height(plan))';
    return;
end
sweep.mode = validatestring(sweep.mode,{'replicates','composition','replicated_composition'});
validateattributes(sweep.seeds,{'numeric'},{'vector','nonempty','real','finite', ...
    'integer','>=',0,'<=',2^32-1});
if numel(unique(sweep.seeds)) ~= numel(sweep.seeds)
    error('NetworkGrowth:DuplicateSeeds','Replicate seeds must be distinct.');
end
validateattributes(baseCfg.fA,{'numeric'},{'scalar','integer','>=',2});
if strcmp(sweep.mode,'replicates')
    fractions = baseCfg.xA;
    if ~isempty(sweep.NtotalScaleByXA)
        error('NetworkGrowth:SweepScale','Replicates use baseCfg.NtotalScale; clear NtotalScaleByXA.');
    end
else
    fractions = sweep.xA(:)';
end
validateattributes(fractions,{'numeric'},{'vector','nonempty','real','finite','>',0,'<',1});
if numel(unique(fractions)) ~= numel(fractions)
    error('NetworkGrowth:DuplicateCompositions','Composition values must be distinct.');
end
seeds = sweep.seeds(:)';
if strcmp(sweep.mode,'composition'), seeds = seeds(1); end
scales = repmat(baseCfg.NtotalScale,size(fractions));
if ~isempty(sweep.NtotalScaleByXA)
    if ~isempty(baseCfg.Ntotal)
        error('NetworkGrowth:SweepScale','NtotalScaleByXA requires baseCfg.Ntotal = [] (computed feed).');
    end
    scales = sweep.NtotalScaleByXA(:)';
    if numel(scales) ~= numel(fractions)
        error('NetworkGrowth:SweepScale','Supply one NtotalScaleByXA value per sweep.xA.');
    end
end
validateattributes(scales,{'numeric'},{'vector','real','finite','positive'});
nRuns = numel(fractions)*numel(seeds);
plan = table('Size',[nRuns 10],'VariableTypes',repmat({'double'},1,10), ...
    'VariableNames',{'runIndex','replicateIndex','compositionIndex','xA','randSeed', ...
    'NtotalScale','Ntotal','initialCrosslinkers','fA','fB'});
k = 0;
for rep = 1:numel(seeds)
    for comp = 1:numel(fractions)
        k = k+1;
        cfg = baseCfg; cfg.xA = fractions(comp); cfg.NtotalScale = scales(comp);
        total = resolveNtotal(cfg);
        validateattributes(total,{'numeric'},{'scalar','finite','integer','positive'});
        initialA = round(cfg.xA*total);
        if initialA < 1
            error('NetworkGrowth:NoSeed','xA = %g and Ntotal = %d provide no seed crosslinker.',cfg.xA,total);
        end
        plan{k,:} = [k rep comp cfg.xA seeds(rep) cfg.NtotalScale total initialA cfg.fA cfg.fB];
    end
end
if strcmp(baseCfg.growthMode,'preformed')
    plan.initialChains = zeros(nRuns,1);
    plan.initialKuhnSegments = zeros(nRuns,1);
    plan.unusedKuhnBudget = zeros(nRuns,1);
    plan.chainEndArmRatio = zeros(nRuns,1);
    for k = 1:nRuns
        cfg = baseCfg; cfg.xA = plan.xA(k); cfg.NtotalScale = plan.NtotalScale(k);
        feed = resolvePreformedFeed(cfg);
        plan.initialChains(k) = sum(feed.initialChainCounts);
        plan.initialKuhnSegments(k) = feed.initialKuhnCount;
        plan.unusedKuhnBudget(k) = feed.unusedKuhnBudget;
        plan.chainEndArmRatio(k) = feed.chainEndArmRatio;
        for s = 1:3
            n = 0; count = 0;
            if s<=numel(feed.chainLengths), n = feed.chainLengths(s); count = feed.initialChainCounts(s); end
            plan.(sprintf('chainLength%d',s))(k) = n;
            plan.(sprintf('initialChainCount%d',s))(k) = count;
        end
    end
end
end
