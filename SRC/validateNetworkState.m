function validateNetworkState(nodes, bonds, activeEnds, fA)
% Fail immediately if topology violates model rules.

N = size(nodes,1);
deg = computeDegrees(bonds, N);

if any(deg > fA)
    error('Topology error: a crosslinker exceeded functionality %d.', fA);
end

if ~isempty(bonds)
    pairs = sort(bonds(:,1:2), 2);
    if size(unique(pairs, 'rows'),1) ~= size(pairs,1)
        error('Topology error: duplicate/parallel bond detected.');
    end
end

expected = sum(fA - deg);
if size(activeEnds,1) ~= expected
    error('Active-end bookkeeping error: have %d, expected %d.', size(activeEnds,1), expected);
end

if ~isempty(activeEnds)
    parent = activeEnds(:,3);
    if any(parent < 1 | parent > N)
        error('Topology error: active end has invalid parent node.');
    end
    if any(deg(parent) >= fA)
        error('Topology error: saturated node still owns an active end.');
    end
end

if N > 1
    adj = makeAdjacency(bonds, N);
    seen = false(N,1);
    seen(1) = true;
    queue = 1;
    while ~isempty(queue)
        i = queue(1);
        queue(1) = [];
        nbr = find(adj(i,:) & ~seen');
        if ~isempty(nbr)
            seen(nbr) = true;
            queue = [queue, nbr];
        end
    end
    if ~all(seen)
        bad = find(~seen);
        error('Topology error: disconnected node(s) detected. First bad node = %d.', bad(1));
    end
end
end

