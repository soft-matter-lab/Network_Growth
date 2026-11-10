function drawMotion(nodes, bonds, L, activeEnds, viz, fA, b)
% Update left-panel geometry during mechanical substeps.

H = buildHeatmap(activeEnds, bonds, nodes, L, b, viz.Xg, viz.Yg);
set(viz.hHeat, 'CData', H);

updateNetworkGeometry(viz,nodes,bonds,L);
end

