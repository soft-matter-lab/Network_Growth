function [ii, jj, rv, d] = repulsionPairs(nodes, bonds, L, rEq, repelBondedPairs, useBinning, fullPairThreshold)
% Return candidate repulsive pairs i<j within cutoff rEq under periodic BC.
% Optionally excludes directly bonded pairs.

N = size(nodes,1);
ii = zeros(0,1);
jj = zeros(0,1);
rv = zeros(0,2);
d = zeros(0,1);

if N < 2 || rEq <= 0
    return;
end

if nargin < 5, repelBondedPairs = true; end
if nargin < 6, useBinning = true; end
if nargin < 7, fullPairThreshold = 120; end

useFullPairs = (~useBinning) || (N <= fullPairThreshold) || (rEq >= 0.5 * L);
if useFullPairs
    [ii, jj] = find(triu(true(N),1));
else
    [ii, jj] = binnedPairsPeriodic(nodes, L, rEq);
end

if isempty(ii)
    return;
end

rv = nodes(jj,1:2) - nodes(ii,1:2);
rv = rv - L * round(rv / L);
d = sqrt(sum(rv.^2,2));

use = d < rEq;
ii = ii(use);
jj = jj(use);
rv = rv(use,:);
d = d(use);

if isempty(ii)
    return;
end

if ~repelBondedPairs
    adj = makeAdjacency(bonds, N);
    keep = ~adj(sub2ind([N,N], ii, jj));
    ii = ii(keep);
    jj = jj(keep);
    rv = rv(keep,:);
    d = d(keep);
end
end

function [ii, jj] = binnedPairsPeriodic(nodes, L, rCut)
% Cell-list pair generation for periodic 2D domain.

N = size(nodes,1);
ii = zeros(0,1);
jj = zeros(0,1);
if N < 2
    return;
end

cellSize = rCut;
nCell = max(1, floor(L / max(cellSize, 1e-12)));

xy = mod(nodes(:,1:2), L);
cx = floor(xy(:,1) / L * nCell) + 1;
cy = floor(xy(:,2) / L * nCell) + 1;
cx = min(max(cx,1), nCell);
cy = min(max(cy,1), nCell);
cellLin = cx + (cy - 1) * nCell;

bucket = accumarray(cellLin, (1:N)', [nCell * nCell, 1], @(v){v}, {[]});

for c = 1:nCell*nCell
    idxA = bucket{c};
    if isempty(idxA)
        continue;
    end
    [ix, iy] = ind2sub([nCell, nCell], c);

    neigh = zeros(9,1);
    kk = 0;
    for dx = -1:1
        jx = mod(ix - 1 + dx, nCell) + 1;
        for dy = -1:1
            jy = mod(iy - 1 + dy, nCell) + 1;
            kk = kk + 1;
            neigh(kk) = sub2ind([nCell, nCell], jx, jy);
        end
    end
    neigh = unique(neigh);

    for n = 1:numel(neigh)
        c2 = neigh(n);
        if c2 < c
            continue;
        end
        idxB = bucket{c2};
        if isempty(idxB)
            continue;
        end

        if c2 == c
            if numel(idxA) < 2
                continue;
            end
            [iLocal, jLocal] = find(triu(true(numel(idxA), numel(idxA)),1));
            p = idxA(iLocal);
            q = idxA(jLocal);
        else
            [P, Q] = ndgrid(idxA, idxB);
            p = P(:);
            q = Q(:);
        end

        ii = [ii; p];
        jj = [jj; q];
    end
end
end
