function [Etotal, Ebond, Erep] = networkEnergy(nodes, bonds, b, kT, kRep, rEq, L, repelBondedPairs, useBinning, fullPairThreshold)
% Energy for the force model used during equilibration.

Ebond = 0;
Erep = 0;
N = size(nodes,1);
if nargin < 8, repelBondedPairs = true; end
if nargin < 9, useBinning = true; end
if nargin < 10, fullPairThreshold = 120; end
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
        Ebond = sum(kT .* (3 * r.^2 ./ (2 * n * b^2)));
    end
end

if kRep > 0 && rEq > 0 && N > 1
    [ii, ~, ~, d] = repulsionPairs(nodes, bonds, L, rEq, repelBondedPairs, useBinning, fullPairThreshold);
    if ~isempty(ii)
        Erep = sum(0.5 * kRep .* (rEq - d).^2);
    end
end

Etotal = Ebond + Erep;
end
