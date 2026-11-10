function viz = initVisualization(cfg, histT, histA, histB, visible)
% Create the two-panel live figure and return handles/grid.
if nargin < 5, visible = 'on'; end

viz = struct();
viz.fig = figure('Color', 'w', 'Position', [50 50 1200 600], 'Visible', visible);

viz.ax1 = subplot(1,2,1);
hold(viz.ax1, 'on');
viz.hHeat = imagesc(viz.ax1, [0 cfg.L], [0 cfg.L], zeros(cfg.heatGrid));
set(viz.ax1, 'YDir', 'normal');
colormap(viz.ax1, jet);
clim(viz.ax1, [0 1.2]);
graphics = initNetworkGraphics(viz.ax1,networkColorStyle(cfg),90,2,[1 1 1]);
viz.colorStyle = graphics.colorStyle;
viz.hChain = graphics.hChain;
viz.hNodes = graphics.hNodes;
viz.hNodeSpheres = graphics.hNodeSpheres;
viz.nodeSize = graphics.nodeSize;
viz.hLoop = scatter(viz.ax1, NaN, NaN, 200, 'MarkerEdgeColor', 'c', 'LineWidth', 1.8,'HandleVisibility','off');
axis(viz.ax1, [0 cfg.L 0 cfg.L]);
axis(viz.ax1, 'square');
axis(viz.ax1, 'off');
title(viz.ax1, 'Polymer Network Formation (Part A, 2D)');
isPreformed = isfield(cfg,'growthMode') && strcmp(cfg.growthMode,'preformed');
if isPreformed, title(viz.ax1,'Preformed-chain end-linking (Case B, 2D)'); end

viz.ax2 = subplot(1,2,2);
hold(viz.ax2, 'on');
viz.hA = plot(viz.ax2, histT, histA, 'r-', 'LineWidth', 2, 'DisplayName', 'Crosslinkers (A) remaining');
if isPreformed
    viz.hB = gobjects(1,numel(cfg.chainLengths));
    for s = 1:numel(cfg.chainLengths)
        viz.hB(s) = plot(viz.ax2,histT,histB(:,s),'LineWidth',2, ...
            'DisplayName',sprintf('Free chains: N = %d',cfg.chainLengths(s)));
    end
else
    viz.hB = plot(viz.ax2, histT, histB, 'b-', 'LineWidth', 2, 'DisplayName', 'Linkers (B) remaining');
end
xlabel(viz.ax2, 'time');
ylabel(viz.ax2, 'count');
grid(viz.ax2, 'on');
legend(viz.ax2, 'Location', 'northeast');

[viz.Xg, viz.Yg] = meshgrid(linspace(0, cfg.L, cfg.heatGrid), linspace(0, cfg.L, cfg.heatGrid));
end

