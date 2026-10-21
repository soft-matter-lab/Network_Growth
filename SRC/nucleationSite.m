function pos = nucleationSite(nodes, cfg)
% Position for a new nucleus in the periodic box.
% With cfg.nucAvoidOccupied, draws cfg.nucCandidates uniform random sites and
% keeps the one farthest (periodic distance) from any existing junction.

K = max(1, round(cfg.nucCandidates));
if ~cfg.nucAvoidOccupied || isempty(nodes)
    K = 1;
end
cand = cfg.L * rand(K, 2);

if K == 1
    pos = cand;
    return;
end

dmin = inf(K, 1);
for k = 1:K
    dv = nodes(:,1:2) - cand(k,:);
    dv = dv - cfg.L * round(dv / cfg.L);
    dmin(k) = min(sum(dv.^2, 2));
end
[~, best] = max(dmin);
pos = cand(best,:);
end
