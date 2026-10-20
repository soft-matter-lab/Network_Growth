function [Etotal, Ebond, Erep] = networkEnergy(nodes, bonds, b, kT, kRep, rEq, L, repelBondedPairs, useBinning, fullPairThreshold, kRepAnh, chainModel)
% Energy for the force model used during equilibration.

Ebond = 0;
Erep = 0;
N = size(nodes,1);
if nargin < 8, repelBondedPairs = true; end
if nargin < 9, useBinning = true; end
if nargin < 10, fullPairThreshold = 120; end
if nargin < 11, kRepAnh = 0; end
if nargin < 12, chainModel = 'cohen'; end
if N < 1
    Etotal = 0;
    return;
end

if ~isempty(bonds)
    sel = bonds(:,1) ~= bonds(:,2);
    if any(sel)
        p = bonds(sel,1);
        q = bonds(sel,2);
        n = max(bonds(sel,3),1);
        img = bonds(sel,4:5);
        rv = nodes(q,1:2) + img * L - nodes(p,1:2);
        r = sqrt(sum(rv.^2,2));
        [~, Eb] = chainSpring(r, n, b, kT, chainModel);
        Ebond = sum(Eb);
    end
end

if kRep > 0 && rEq > 0 && N > 1
    [ii, ~, ~, d] = repulsionPairs(nodes, bonds, L, rEq, repelBondedPairs, useBinning, fullPairThreshold);
    if ~isempty(ii)
        s = rEq - d;
        Erep = sum(0.5 * kRep .* s.^2 + kRepAnh .* s.^3 / 3);
    end
end

Etotal = Ebond + Erep;
end
