function Ntotal = resolveNtotal(cfg)
% Total feed count from box size and Kuhn length:
%   areaFraction * L^2 = Ntotal * pi * b^2 / 4
% (crosslinkers and linkers are assumed to occupy the same area).
% The computed value is multiplied by cfg.NtotalScale (default 1).
% A non-empty cfg.Ntotal is used as given (no scaling, so repeated calls are safe).

if isfield(cfg, 'Ntotal') && ~isempty(cfg.Ntotal)
    Ntotal = cfg.Ntotal;
else
    scale = 1;
    if isfield(cfg, 'NtotalScale') && ~isempty(cfg.NtotalScale)
        scale = cfg.NtotalScale;
    end
    Ntotal = max(1, round(scale * cfg.areaFraction * 4 * cfg.L^2 / (pi * cfg.b^2)));
end
end
