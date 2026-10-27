function [dv, available] = networkDumpDisplacements(s)
% First stored bond endpoint -> second endpoint, with the recorded periodic image.
% Translate only: do not wrap again, rotate, or duplicate vectors with -dv.
available = true;
if isfield(s,'completed') && isfield(s.completed,'endToEndVector')
    dv = s.completed.endToEndVector;
elseif isfield(s,'nodes') && isfield(s,'bonds')
    p = s.bonds(:,1); q = s.bonds(:,2);
    dv = s.nodes(q,1:2) + s.bonds(:,4:5)*s.L - s.nodes(p,1:2);
else
    dv = zeros(0,2); available = false;
end
end
