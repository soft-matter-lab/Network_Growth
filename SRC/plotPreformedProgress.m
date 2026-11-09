function [fig, progress] = plotPreformedProgress(runOrDumpDir, options)
% Per-species background counts and captured/completed-chain conversion.
if nargin<2, options = struct(); end
options = mergeConfig(struct('visible','on','outputFile',''),options);
[files,metadata,dumpDir] = listNetworkDumpFiles(runOrDumpDir);
if ~isfield(metadata,'preformedFeed')
    error('NetworkGrowth:PreformedProgress','This run has no preformed-chain reservoir data.');
end
lengths = metadata.preformedFeed.chainLengths;
count = numel(files); species = numel(lengths);
time = zeros(count,1); event = time; freeA = time;
free = zeros(count,species); captured = free; completed = free;
for k = 1:count
    data = load(fullfile(dumpDir,files(k).name),'snapshot'); s = data.snapshot;
    time(k) = s.ssaTime; event(k) = s.event; freeA(k) = s.freeCrosslinkers;
    free(k,:) = s.preformed.freeChainCounts;
    captured(k,:) = s.preformed.capturedPercentBySpecies;
    completed(k,:) = s.preformed.completedPercentBySpecies;
end
progress = table(event,time,freeA,'VariableNames',{'event','ssaTime','freeCrosslinkers'});
fig = figure('Color','w','Visible',options.visible,'Position',[100 100 1200 450]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact');
ax1 = nexttile(layout); hold(ax1,'on'); ax2 = nexttile(layout); hold(ax2,'on');
colors = lines(species);
for j = 1:species
    label = sprintf('N = %d',lengths(j));
    plot(ax1,time,free(:,j),'Color',colors(j,:),'LineWidth',1.6,'DisplayName',label);
    plot(ax2,time,captured(:,j),'Color',colors(j,:),'LineWidth',1.6,'DisplayName',[label ' captured']);
    plot(ax2,time,completed(:,j),'--','Color',colors(j,:),'LineWidth',1.6,'DisplayName',[label ' completed']);
    progress.(sprintf('freeChains_N%d',lengths(j))) = free(:,j);
    progress.(sprintf('capturedPercent_N%d',lengths(j))) = captured(:,j);
    progress.(sprintf('completedPercent_N%d',lengths(j))) = completed(:,j);
end
xlabel(ax1,'SSA time'); ylabel(ax1,'Background chains'); title(ax1,'Remaining preformed chains');
xlabel(ax2,'SSA time'); ylabel(ax2,'Percent of initial species count'); ylim(ax2,[0 100]);
title(ax2,'Captured: one or both ends attached; completed: both');
grid(ax1,'on'); grid(ax2,'on'); legend(ax1,'Location','best'); legend(ax2,'Location','best');
if ~isempty(options.outputFile), exportgraphics(fig,options.outputFile,'Resolution',180); end
end
