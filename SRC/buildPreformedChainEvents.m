function events = buildPreformedChainEvents(nodes,bonds,activeEnds,nA,chainCounts,cfg)
% [type,end,species-or-destination,propensity]. No growth or dangling joining.
events = zeros(0,4);
for e = 1:size(activeEnds,1)
    if activeEnds(e,4)==0
        for s = 1:numel(chainCounts)
            rate = cfg.kAddition*chainCounts(s)/cfg.L^2;
            if rate>0, events(end+1,:) = [1,e,s,rate]; end %#ok<AGROW>
        end
    elseif nA>0 && cfg.kChem>0
        events(end+1,:) = [2,e,0,cfg.fA*cfg.kChem*nA/cfg.L^2]; %#ok<AGROW>
    end
end
adj = makeAdjacency(bonds,size(nodes,1));
bare = find(activeEnds(:,4)==0);
dangling = find(activeEnds(:,4)>0);
for e = dangling'
    p = activeEnds(e,3); n = activeEnds(e,4);
    for o = bare'
        q = activeEnds(o,3);
        if adj(p,q), continue; end % one strand per pair, one loop per anchor
        if p==q
            if ~cfg.allowSelfLoops, continue; end
            d2 = cfg.selfLoopDistanceFactor^2*n*cfg.b^2;
            multiplier = cfg.kSelfLoop;
        else
            dv = nodes(q,1:2)-nodes(p,1:2);
            dv = dv-cfg.L*round(dv/cfg.L);
            d2 = sum(dv.^2); multiplier = cfg.kEndEnd;
        end
        rate = multiplier*cfg.kChem/(pi*n*cfg.b^2)*exp(-d2/(n*cfg.b^2));
        if rate>0, events(end+1,:) = [3,e,o,rate]; end %#ok<AGROW>
    end
end
end
