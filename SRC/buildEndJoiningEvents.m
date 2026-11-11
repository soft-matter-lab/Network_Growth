function events = buildEndJoiningEvents(nodes, bonds, activeEnds, cfg)
% Type-3 SSA propensities, enumerated once per unordered active-end pair.
% cfg is a resolved configuration (defaultNetworkAConfig + overrides).
% Same-anchor modes: growing chain + bare arm, or two growing chains.
% A bare arm has n=0; no loop can use an already occupied anchor site.
events = zeros(0,4); % [type, end1, end2, weight]
M = size(activeEnds,1);
if M < 2, return; end
adjNow = makeAdjacency(bonds,size(nodes,1));
% Keep e-major/o-major ordering: the SSA cumulative sum depends on it.
% Vectorize only destinations for one source, avoiding an M-by-M workspace.
chunks = cell(M-1,1);
for e = 1:M-1
    p1 = activeEnds(e,3); n1 = activeEnds(e,4);
    o = (e+1:M)';
    p2 = activeEnds(o,3);
    nTot = n1+activeEnds(o,4);
    isLoop = p1 == p2;
    keep = ~adjNow(p1,p2)' & nTot>0 & nTot>=2*cfg.minKuhnEndEnd;
    if ~cfg.allowSelfLoops || cfg.kSelfLoop == 0
        keep = keep & ~isLoop;
    end
    o = o(keep); p2 = p2(keep); nTot = nTot(keep); isLoop = isLoop(keep);
    if isempty(o), chunks{e} = zeros(0,4); continue; end
    dv = nodes(p2,1:2)-nodes(p1,1:2);
    dv = dv-cfg.L*round(dv/cfg.L);
    d2 = sum(dv.^2,2);
    kJoin = repmat(cfg.kEndEnd,size(o));
    if any(isLoop)
        % Rate-only distance penalty for chain+bare and chain+chain loops.
        d2(isLoop) = cfg.selfLoopDistanceFactor^2*nTot(isLoop)*cfg.b^2;
        kJoin(isLoop) = cfg.kSelfLoop;
    end
    % Preserve the scalar formula's operation order and existing cutoff.
    weight = kJoin*cfg.kChem./(pi*nTot*cfg.b^2).*exp(-d2./(nTot*cfg.b^2));
    keep = weight>1e-14;
    o = o(keep); weight = weight(keep);
    chunks{e} = [repmat(3,size(o)),repmat(e,size(o)),o,weight];
end
events = vertcat(chunks{:});
if isempty(events), events = zeros(0,4); end
end
