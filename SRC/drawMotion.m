function drawMotion(nodes, bonds, L, activeEnds, viz, fA, b)
% Update left-panel geometry during mechanical substeps.

H = buildHeatmap(activeEnds, bonds, nodes, L, b, viz.Xg, viz.Yg);
set(viz.hHeat, 'CData', H);

[skX, skY] = buildSkeleton(nodes, bonds, L);
set(viz.hChain, 'XData', skX, 'YData', skY);

deg = computeDegrees(bonds, size(nodes,1));
nodesWrap = wrapPoint(nodes(:,1:2), L);
set(viz.hPink, 'XData', nodesWrap(deg <= 1,1), 'YData', nodesWrap(deg <= 1,2));
set(viz.hRed, 'XData', nodesWrap(deg >= 2 & deg < fA,1), 'YData', nodesWrap(deg >= 2 & deg < fA,2));
set(viz.hBlack, 'XData', nodesWrap(deg >= fA,1), 'YData', nodesWrap(deg >= fA,2));

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
end

