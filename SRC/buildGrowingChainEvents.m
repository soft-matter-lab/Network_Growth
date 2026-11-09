function events = buildGrowingChainEvents(nodes,bonds,activeEnds,nA,nB,cfg)
% Case-A event ordering is preserved for identical seeded trajectories.
events = zeros(0,4);
for e = 1:size(activeEnds,1)
    n = activeEnds(e,4);
    if nB>0, events(end+1,:) = [1,e,0,cfg.kChem*(nB/cfg.L^2)]; end %#ok<AGROW>
    if nA>0 && n>=cfg.minKuhnCrosslink
        events(end+1,:) = [2,e,0,cfg.fA*cfg.kChem*(nA/cfg.L^2)]; %#ok<AGROW>
    end
end
events = [events;buildEndJoiningEvents(nodes,bonds,activeEnds,cfg)];
if cfg.nucleation && nA>0 && nB>0 && cfg.kNuc>0
    events(end+1,:) = [4,0,0,cfg.kNuc*nA*nB/cfg.L^2];
end
end
