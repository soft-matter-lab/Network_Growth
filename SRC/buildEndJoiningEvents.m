function events = buildEndJoiningEvents(nodes, bonds, activeEnds, cfg)
% Type-3 SSA propensities, enumerated once per unordered active-end pair.
% cfg is a resolved configuration (defaultNetworkAConfig + overrides).
% Same-anchor modes: growing chain + bare arm, or two growing chains.
% A bare arm has n=0; no loop can use an already occupied anchor site.
events = zeros(0,4); % [type, end1, end2, weight]
M = size(activeEnds,1);
if M < 2, return; end
adjNow = makeAdjacency(bonds,size(nodes,1));
for e = 1:M-1
    p1 = activeEnds(e,3); n1 = activeEnds(e,4);
    for o = e+1:M
        p2 = activeEnds(o,3); n2 = activeEnds(o,4);
        if adjNow(p1,p2), continue; end
        isLoop = p1 == p2;
        if isLoop && (~cfg.allowSelfLoops || cfg.kSelfLoop == 0), continue; end
        nTot = n1+n2;
        if nTot <= 0 || nTot < 2*cfg.minKuhnEndEnd, continue; end
        if isLoop
            % Deliberate kinetic penalty, not the actual anchor displacement.
            % For a bare destination n2=0 this is factor*b*sqrt(n1); for
            % two dangling chains the squared coil sizes add: b^2*(n1+n2).
            d2 = cfg.selfLoopDistanceFactor^2*nTot*cfg.b^2;
            kJoin = cfg.kSelfLoop;
        else
            dv = nodes(p2,1:2)-nodes(p1,1:2);
            dv = dv-cfg.L*round(dv/cfg.L);
            d2 = sum(dv.^2);
            kJoin = cfg.kEndEnd;
        end
        weight = kJoin*cfg.kChem/(pi*nTot*cfg.b^2)*exp(-d2/(nTot*cfg.b^2));
        if weight > 1e-14
            events(end+1,:) = [3,e,o,weight]; %#ok<AGROW>
        end
    end
end
end
