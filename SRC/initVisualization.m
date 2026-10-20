function viz = initVisualization(cfg, histT, histA, histB)
% Create the two-panel live figure and return handles/grid.

viz = struct();
viz.fig = figure('Color', 'w', 'Position', [50 50 1200 600]);

viz.ax1 = subplot(1,2,1);
hold(viz.ax1, 'on');
viz.hHeat = imagesc(viz.ax1, [0 cfg.L], [0 cfg.L], zeros(cfg.heatGrid));
set(viz.ax1, 'YDir', 'normal');
colormap(viz.ax1, jet);
clim(viz.ax1, [0 1.2]);
viz.hChain = plot(viz.ax1, NaN, NaN, 'w-', 'LineWidth', 2);
viz.hPink = scatter(viz.ax1, NaN, NaN, 90, [1 0.7 0.9], 'filled', 'MarkerEdgeColor', 'w');
viz.hRed = scatter(viz.ax1, NaN, NaN, 90, 'r', 'filled', 'MarkerEdgeColor', 'w');
viz.hBlack = scatter(viz.ax1, NaN, NaN, 90, 'k', 'filled', 'MarkerEdgeColor', 'w');
viz.hLoop = scatter(viz.ax1, NaN, NaN, 200, 'MarkerEdgeColor', 'c', 'LineWidth', 1.8);
axis(viz.ax1, [0 cfg.L 0 cfg.L]);
axis(viz.ax1, 'square');
axis(viz.ax1, 'off');
title(viz.ax1, 'Polymer Network Formation (Part A, 2D)');

viz.ax2 = subplot(1,2,2);
hold(viz.ax2, 'on');
viz.hA = plot(viz.ax2, histT, histA, 'r-', 'LineWidth', 2, 'DisplayName', 'Crosslinkers (A) remaining');
viz.hB = plot(viz.ax2, histT, histB, 'b-', 'LineWidth', 2, 'DisplayName', 'Linkers (B) remaining');
xlabel(viz.ax2, 'time');
ylabel(viz.ax2, 'count');
grid(viz.ax2, 'on');
legend(viz.ax2, 'Location', 'northeast');

[viz.Xg, viz.Yg] = meshgrid(linspace(0, cfg.L, cfg.heatGrid), linspace(0, cfg.L, cfg.heatGrid));
end

