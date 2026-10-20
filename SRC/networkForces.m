function F = networkForces(nodes, bonds, b, kT, kRep, rEq, L, repelBondedPairs, useBinning, fullPairThreshold, kRepAnh)
% Net force on each junction:
% (1) tensile bonded Gaussian spring
% (2) repulsion within rEq: force magnitude kRep*s + kRepAnh*s^2, s = rEq - d
%     (kRepAnh > 0 adds an anharmonic stiffening; 0 recovers harmonic)
N = size(nodes,1);
F = zeros(N,2);

if nargin < 8, repelBondedPairs = true; end
if nargin < 9, useBinning = true; end
if nargin < 10, fullPairThreshold = 120; end
if nargin < 11, kRepAnh = 0; end

if ~isempty(bonds)
    sel = bonds(:,1) ~= bonds(:,2);
    if any(sel)
        p = bonds(sel,1);
        q = bonds(sel,2);
        n = max(bonds(sel,3),1);
        img = bonds(sel,4:5);
        rv = nodes(q,1:2) + img * L - nodes(p,1:2);
        c = 3 * kT ./ (n * b^2);
        fv = rv .* c;
        F(:,1) = F(:,1) + accumarray([p; q], [ fv(:,1); -fv(:,1)], [N,1]);
        F(:,2) = F(:,2) + accumarray([p; q], [ fv(:,2); -fv(:,2)], [N,1]);
    end
end

if kRep > 0 && rEq > 0 && N > 1
    [ii, jj, rv, d] = repulsionPairs(nodes, bonds, L, rEq, repelBondedPairs, useBinning, fullPairThreshold);
    if ~isempty(ii)
        coincident = d <= 1e-9 * b;
        if any(coincident)
            theta0 = 2 * pi * mod(ii(coincident) + sqrt(2) * jj(coincident), 1);
            epsSep = 1e-9 * b;
            rv(coincident,:) = epsSep * [cos(theta0), sin(theta0)];
            d(coincident) = epsSep;
        end

        s = rEq - d;
        fe = rv .* ((kRep * s + kRepAnh * s.^2) ./ d);
        F(:,1) = F(:,1) + accumarray([ii; jj], [-fe(:,1); fe(:,1)], [N,1]);
        F(:,2) = F(:,2) + accumarray([ii; jj], [-fe(:,2); fe(:,2)], [N,1]);
    end
end
end
