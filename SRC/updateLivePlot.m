function updateLivePlot(viz, nodes, bonds, activeEnds, cfg, histT, histA, histB)
% Update both left and right panels for SSA-level visualization.

if isfield(viz,'showDanglingClouds')
    % Offline replay uses truecolor so its blue clouds do not recolor junctions.
    if viz.showDanglingClouds
        rgb = danglingCloudColors(activeEnds,nodes,cfg.L,cfg.b,viz.Xg,viz.Yg);
        set(viz.hHeat,'CData',rgb,'Visible','on');
    else
        set(viz.hHeat,'Visible','off');
    end
else
    H = buildHeatmap(activeEnds, bonds, nodes, cfg.L, cfg.b, viz.Xg, viz.Yg);
    set(viz.hHeat,'CData',H);
end
[skX, skY] = buildSkeleton(nodes, bonds, cfg.L);
set(viz.hChain, 'XData', skX, 'YData', skY);

deg = computeDegrees(bonds, size(nodes,1));
ip = find(deg <= 1);
ir = find(deg >= 2 & deg < cfg.fA);
ib = find(deg >= cfg.fA);
nodesWrap = wrapPoint(nodes(:,1:2), cfg.L);

set(viz.hPink, 'XData', nodesWrap(ip,1), 'YData', nodesWrap(ip,2));
set(viz.hRed, 'XData', nodesWrap(ir,1), 'YData', nodesWrap(ir,2));
set(viz.hBlack, 'XData', nodesWrap(ib,1), 'YData', nodesWrap(ib,2));

if isempty(bonds)
    il = [];
else
    il = unique(bonds(bonds(:,1) == bonds(:,2),1));
end
if isempty(il)
    set(viz.hLoop, 'XData', NaN, 'YData', NaN);
else
    set(viz.hLoop, 'XData', nodesWrap(il,1), 'YData', nodesWrap(il,2));
end

set(viz.hA, 'XData', histT, 'YData', histA);
set(viz.hB, 'XData', histT, 'YData', histB);
end

