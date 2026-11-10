function updateNetworkGeometry(graphics, nodes, bonds, L)
% Preserve image-aware periodic geometry while applying stable visual colors.
style = graphics.colorStyle;
for s = 1:numel(graphics.hChain)
    selected = bonds;
    if ~isempty(style.chainLengths), selected = bonds(bonds(:,3)==style.chainLengths(s),:); end
    [x,y] = buildSkeleton(nodes,selected,L);
    set(graphics.hChain(s),'XData',x,'YData',y);
end
degree = computeDegrees(bonds,size(nodes,1));
xy = wrapPoint(nodes(:,1:2),L);
rgb = networkNodeColors(degree,style.fA,style.showNodeConnectivity);
set(graphics.hNodes,'XData',xy(:,1),'YData',xy(:,2),'CData',rgb);
if style.showSphericalNodes
    updateSphereNodes(graphics.hNodeSpheres,xy,rgb,graphics.nodeSize,L);
end
if isfield(graphics,'hLoop')
    loops = unique(bonds(bonds(:,1)==bonds(:,2),1));
    if isempty(loops), set(graphics.hLoop,'XData',NaN,'YData',NaN);
    else, set(graphics.hLoop,'XData',xy(loops,1),'YData',xy(loops,2)); end
end
end
