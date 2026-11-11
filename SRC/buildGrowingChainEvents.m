function events = buildGrowingChainEvents(nodes,bonds,activeEnds,nA,nB,cfg)
% Case-A event ordering is preserved for identical seeded trajectories.
% Reserve the maximum reservoir-event count, then trim. Keep addition and
% fresh-crosslink events interleaved in the original active-end order.
events = zeros(2*size(activeEnds,1),4);
count = 0;
additionRate = cfg.kChem*(nB/cfg.L^2);
freshRate = cfg.fA*cfg.kChem*(nA/cfg.L^2);
for e = 1:size(activeEnds,1)
    n = activeEnds(e,4);
    if nB>0
        count = count+1; events(count,:) = [1,e,0,additionRate];
    end
    if nA>0 && n>=cfg.minKuhnCrosslink
        count = count+1; events(count,:) = [2,e,0,freshRate];
    end
end
events = events(1:count,:);
events = [events;buildEndJoiningEvents(nodes,bonds,activeEnds,cfg)];
if cfg.nucleation && nA>0 && nB>0 && cfg.kNuc>0
    events(end+1,:) = [4,0,0,cfg.kNuc*nA*nB/cfg.L^2];
end
end
