function Ntotal = resolveNtotal(cfg)
% Total feed count from box size and Kuhn length:
%   areaFraction * L^2 = Ntotal * pi * b^2 / 4
% (crosslinkers and linkers are assumed to occupy the same area).
% A non-empty cfg.Ntotal is used as given.

if isfield(cfg, 'Ntotal') && ~isempty(cfg.Ntotal)
    Ntotal = cfg.Ntotal;
else
    Ntotal = max(1, round(cfg.areaFraction * 4 * cfg.L^2 / (pi * cfg.b^2)));
end
end
