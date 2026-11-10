function graphics = initNetworkGraphics(ax, style, nodeSize, bondWidth, bondColor)
% Shared node/bond artists for live, replay, and checkpoint network panels.
graphics.colorStyle = style;
if style.isPreformed, bondColor = [0 0 0]; end
count = max(1,numel(style.chainLengths));
graphics.hChain = gobjects(1,count);
for s = 1:count
    color = bondColor; label = 'Strands';
    if ~isempty(style.chainLengths)
        color = style.bondColors(s,:); label = sprintf('N = %d',style.chainLengths(s));
    end
    graphics.hChain(s) = plot(ax,NaN,NaN,'-','Color',color,'LineWidth',bondWidth, ...
        'DisplayName',label,'Tag','networkBonds');
end
graphics.hNodes = scatter(ax,NaN,NaN,nodeSize,[0 0 0],'filled', ...
    'MarkerEdgeColor','none','HandleVisibility','off','Tag','networkNodes');
graphics.nodeSize = nodeSize;
graphics.hNodeSpheres = patch(ax,'Faces',zeros(0,3),'Vertices',zeros(0,2), ...
    'FaceVertexCData',zeros(0,3),'FaceColor','interp','EdgeColor','none', ...
    'FaceLighting','none','HandleVisibility','off','Tag','networkNodeSpheres');
if style.showSphericalNodes
    set(graphics.hNodes,'Visible','off');
else
    set(graphics.hNodeSpheres,'Visible','off');
end
if numel(style.chainLengths)>1
    legend(ax,graphics.hChain,'Location','northeast','AutoUpdate','off');
end
end
