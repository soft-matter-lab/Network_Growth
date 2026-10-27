function [fig, rdf] = plotNetworkDump(snapshotOrFile, options)
% Plot one saved state: periodic network and per-snapshot probability histograms.
% Pass histogramStyle from networkDumpHistogramStyle for shared bins/axes.
% A filename inside a dumps folder automatically uses that run's shared axes.
% The optional second output returns numerical RDF data only; no RDF is plotted.
if nargin < 2, options = struct(); end
options = mergeConfig(struct('visible','on','outputFile','','rdfEdges',[], ...
    'histogramStyle',[]),options);
if isstruct(snapshotOrFile)
    s = snapshotOrFile;
else
    loaded = load(snapshotOrFile,'snapshot'); s = loaded.snapshot;
end
if isempty(options.histogramStyle)
    source = s;
    if ~isstruct(snapshotOrFile)
        folder = fileparts(char(snapshotOrFile));
        if isempty(folder), folder = pwd; end
        if isfile(fullfile(folder,'metadata.mat')), source = folder; end
    end
    options.histogramStyle = networkDumpHistogramStyle(source,options);
end
style = options.histogramStyle;
fig = figure('Color','w','Visible',options.visible,'Position',[100 100 1100 750]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact');
ax = nexttile(layout,1); rdf = [];
if isfield(s,'nodes')
    [x,y] = buildSkeleton(s.nodes,s.bonds,s.L);
    plot(ax,x,y,'-','Color',[.35 .35 .35]); hold(ax,'on');
    xy = wrapPoint(s.nodes(:,1:2),s.L);
    scatter(ax,xy(:,1),xy(:,2),35,s.nodes(:,3),'filled');
    axis(ax,'equal'); xlim(ax,[0 s.L]); ylim(ax,[0 s.L]);
    cb = colorbar(ax); cb.Label.String = 'junction degree';
    xlabel(ax,'x'); ylabel(ax,'y');
    title(ax,sprintf('Network: %d junctions, %d chains',s.nodeCount,s.bondCount));
    if nargout > 1, rdf = networkRdf(s.nodes,s.L,options.rdfEdges); end
else
    text(ax,.5,.5,'Network payload was disabled','HorizontalAlignment','center'); axis(ax,'off');
end
ax.Tag = 'network';
ax = nexttile(layout,2);
if isfield(s,'completed')
    histogram(ax,s.completed.kuhnCount,style.nEdges,'Normalization','probability');
    xlabel(ax,'n (Kuhn segments)'); ylabel(ax,'probability per bin');
    title(ax,sprintf('Chain number fraction (N = %d)',numel(s.completed.kuhnCount))); grid(ax,'on');
    if isempty(s.completed.kuhnCount)
        text(ax,.5,.5,'No samples','Units','normalized','HorizontalAlignment','center');
    end
    xlim(ax,style.kuhnXLim); ylim(ax,style.kuhnYLim);
else
    text(ax,.5,.5,'Distribution payload was disabled','HorizontalAlignment','center'); axis(ax,'off');
end
ax.Tag = 'kuhnHistogram';
ax = nexttile(layout,3);
if isfield(s,'completed')
    [fraction,totalSegments] = chainWeightFractions(s.completed.kuhnCount,style.weightEdges);
    histogram(ax,'BinEdges',style.weightEdges,'BinCounts',fraction,'Normalization','count');
    xlabel(ax,'n (Kuhn segments)'); ylabel(ax,'fraction of completed-chain Kuhn segments');
    title(ax,sprintf('Chain weight fraction (total n = %g)',totalSegments)); grid(ax,'on');
    if totalSegments==0
        text(ax,.5,.5,'No completed-chain segments','Units','normalized','HorizontalAlignment','center');
    end
    xlim(ax,style.weightXLim); ylim(ax,style.weightYLim);
else
    text(ax,.5,.5,'Distribution payload was disabled','HorizontalAlignment','center'); axis(ax,'off');
end
ax.Tag = 'weightFractionHistogram';
ax = nexttile(layout,4);
[dv,available] = networkDumpDisplacements(s);
if available
    edges = style.displacementEdges;
    centers = (edges(1:end-1)+edges(2:end))/2;
    probability = zeros(numel(centers));
    if ~isempty(dv)
        probability = histcounts2(dv(:,1),dv(:,2),edges,edges)/size(dv,1);
    end
    imagesc(ax,centers,centers,probability');
    set(ax,'YDir','normal'); axis(ax,'equal');
    xlim(ax,style.displacementLim); ylim(ax,style.displacementLim);
    clim(ax,style.displacementCLim); colormap(ax,parula);
    cb = colorbar(ax); cb.Label.String = 'probability per 2D bin';
    xline(ax,0,':','Color',[.7 .7 .7]); yline(ax,0,':','Color',[.7 .7 .7]);
    xlabel(ax,'\Delta x'); ylabel(ax,'\Delta y');
    title(ax,sprintf('End-to-end vectors (N = %d)',size(dv,1)));
    if isempty(dv)
        text(ax,.5,.5,'No bonds','Units','normalized','HorizontalAlignment','center','Color','w');
    end
else
    text(ax,.5,.5,{'Vector data unavailable','This older dump needs the network payload'}, ...
        'Units','normalized','HorizontalAlignment','center'); axis(ax,'off');
end
ax.Tag = 'displacementHistogram';
title(layout,sprintf('A usage %.3g%% | event %d | t = %.5g | %s', ...
    s.usagePercent,s.event,s.ssaTime,s.phase),'Interpreter','none');
if ~isempty(options.outputFile), exportgraphics(fig,options.outputFile,'Resolution',180); end
end
