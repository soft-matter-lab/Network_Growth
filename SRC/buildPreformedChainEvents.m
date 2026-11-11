function events = buildPreformedChainEvents(nodes,bonds,activeEnds,nA,chainCounts,cfg)
% [type,end,species-or-destination,propensity]. No growth or dangling joining.
events = zeros(max(1,numel(chainCounts))*size(activeEnds,1),4);
count = 0;
additionRates = cfg.kAddition*chainCounts/cfg.L^2;
freshRate = cfg.fA*cfg.kChem*nA/cfg.L^2;
for e = 1:size(activeEnds,1)
    if activeEnds(e,4)==0
        for s = 1:numel(chainCounts)
            rate = additionRates(s);
            if rate>0
                count = count+1; events(count,:) = [1,e,s,rate];
            end
        end
    elseif nA>0 && cfg.kChem>0
        count = count+1; events(count,:) = [2,e,0,freshRate];
    end
end
events = events(1:count,:);
adj = makeAdjacency(bonds,size(nodes,1));
bare = find(activeEnds(:,4)==0);
dangling = find(activeEnds(:,4)>0);
chunks = cell(numel(dangling),1);
for index = 1:numel(dangling)
    e = dangling(index);
    p = activeEnds(e,3); n = activeEnds(e,4);
    q = activeEnds(bare,3);
    isLoop = p==q;
    keep = ~adj(p,q)'; % one strand per pair, one loop per anchor
    if ~cfg.allowSelfLoops, keep = keep & ~isLoop; end
    o = bare(keep); q = q(keep); isLoop = isLoop(keep);
    if isempty(o), chunks{index} = zeros(0,4); continue; end
    dv = nodes(q,1:2)-nodes(p,1:2);
    dv = dv-cfg.L*round(dv/cfg.L);
    d2 = sum(dv.^2,2);
    multiplier = repmat(cfg.kEndEnd,size(o));
    if any(isLoop)
        d2(isLoop) = cfg.selfLoopDistanceFactor^2*n*cfg.b^2;
        multiplier(isLoop) = cfg.kSelfLoop;
    end
    % Case B deliberately retains all positive rates (no Case-A cutoff).
    rate = multiplier*cfg.kChem/(pi*n*cfg.b^2).*exp(-d2/(n*cfg.b^2));
    keep = rate>0; o = o(keep); rate = rate(keep);
    chunks{index} = [repmat(3,size(o)),repmat(e,size(o)),o,rate];
end
% Source order, destination order and channel order remain unchanged.
events = [events;vertcat(chunks{:})];
end
