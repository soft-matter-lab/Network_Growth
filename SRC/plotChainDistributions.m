function stats = plotChainDistributions(nodes, bonds, cfg, outputDir, showFigure)
% Chain contour length (Kuhn segments) and end-to-end distance distributions.
% Saves chain_distributions.png and chain_lengths.csv in outputDir when given.

stats = struct('n', [], 'r', []);
if isempty(bonds)
    return;
end

p = bonds(:,1);
q = bonds(:,2);
n = bonds(:,3);
rv = nodes(q,1:2) + bonds(:,4:5) * cfg.L - nodes(p,1:2);
r = sqrt(sum(rv.^2, 2));

stats.n = n;
stats.r = r;

vis = 'off';
if showFigure
    vis = 'on';
end
fig = figure('Color', 'w', 'Visible', vis, 'Position', [100 100 1000 400]);

subplot(1,2,1);
histogram(n, 'BinMethod', 'integers', 'FaceColor', [0.3 0.6 0.9]);
xlabel('contour length n (Kuhn segments)');
ylabel('count');
grid on;
title(sprintf('Kuhn count (mean = %.1f)', mean(n)));

subplot(1,2,2);
histogram(r, 'FaceColor', [0.9 0.5 0.3]);
xlabel('end-to-end distance r');
ylabel('count');
grid on;
title(sprintf('End-to-end (mean = %.2f, \\langle r^2\\rangle^{1/2} = %.2f)', ...
    mean(r), sqrt(mean(r.^2))));

if nargin >= 4 && ~isempty(outputDir)
    exportgraphics(fig, fullfile(outputDir, 'chain_distributions.png'), 'Resolution', 200);
    writetable(table(p, q, n, r, 'VariableNames', {'node1','node2','kuhnCount','endToEnd'}), ...
        fullfile(outputDir, 'chain_lengths.csv'));
end

if ~showFigure
    close(fig);
end
end
