function [lambda, available] = networkDumpPrestretch(s)
% Fraction of contour length: lambda_i = r_i/(n_i*b), for completed bonds.
% Existing magnitude-only dumps are sufficient; network-only dumps use the
% stored bond images, never a new minimum image. Self-loops remain at zero.
% Dangling chains have no measured tip position and supply no samples.
lambda = zeros(0,1); available = false;
hasB = isfield(s,'b') && isnumeric(s.b) && isscalar(s.b) && ...
    isfinite(s.b) && isreal(s.b) && s.b > 0;
if isfield(s,'completed') && all(isfield(s.completed,{'kuhnCount','endToEnd'}))
    n = s.completed.kuhnCount(:);
    r = s.completed.endToEnd(:);
    if isempty(n) && isempty(r), available = true; return; end
    if hasB
        contour = n*s.b;
    elseif isfield(s.completed,'contourLength')
        contour = s.completed.contourLength(:);
    else
        return;
    end
elseif hasB && all(isfield(s,{'nodes','bonds','L'}))
    n = s.bonds(:,3);
    p = s.bonds(:,1); q = s.bonds(:,2);
    dv = s.nodes(q,1:2) + s.bonds(:,4:5)*s.L - s.nodes(p,1:2);
    r = sqrt(sum(dv.^2,2));
    contour = n*s.b;
else
    return;
end
% Do not pair arrays from incompatible populations or invent missing lengths.
if numel(n) ~= numel(r) || numel(n) ~= numel(contour), return; end
available = true;
valid = isfinite(n) & n > 0 & isfinite(contour) & contour > 0 & isfinite(r) & r >= 0;
lambda = r(valid)./contour(valid);
lambda = lambda(isfinite(lambda)); % no clipping at lambda=1
end
