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
updateNetworkGeometry(viz,nodes,bonds,cfg.L);

set(viz.hA, 'XData', histT, 'YData', histA);
if isfield(cfg,'growthMode') && strcmp(cfg.growthMode,'preformed')
    for s = 1:numel(viz.hB), set(viz.hB(s),'XData',histT,'YData',histB(:,s)); end
else
    set(viz.hB, 'XData', histT, 'YData', histB);
end
end

