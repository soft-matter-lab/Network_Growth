function [fig, selection] = plotNetworkDumps(runDirs, options)
% Compare distributions at conversion milestones, within/across runs.
% runDirs: one run/dumps directory, or a cell/string array of directories.
% options: usagePercent, population ('completed','dangling','all'), labels,
% nEdges, contourEdges, rEdges, visible, outputFile. Shared probability bins.
% End-to-end distances are defined only for completed chains in this model.
if nargin < 2, options = struct(); end
defaults = struct('usagePercent',[0 25 50 75 100],'population','completed', ...
    'labels',{{}},'nEdges',[],'contourEdges',[],'rEdges',[], ...
    'visible','on','outputFile','');
options = mergeConfig(defaults,options);
population = validatestring(options.population,{'completed','dangling','all'});
validateattributes(options.usagePercent,{'numeric'},{'vector','finite','>=',0,'<=',100});
if ischar(runDirs), runDirs = {runDirs}; else, runDirs = cellstr(runDirs); end
labels = options.labels;
if isempty(labels), labels = compose('run %d',1:numel(runDirs)); end
labels = cellstr(labels);
if numel(labels) ~= numel(runDirs)
    error('NetworkGrowth:PlotLabels','Supply one label per run directory.');
end
selection = table(); samples = cell(0,3); legendLabels = {};
for j = 1:numel(runDirs)
    [snapshots,~,index] = loadNetworkDumps(runDirs{j},options.usagePercent);
    for target = unique(options.usagePercent(:)')
        k = find(cellfun(@(s)any(abs(s.usageTargets-target)<1e-9),snapshots),1);
        if isempty(k)
            warning('NetworkGrowth:MissingMilestone','%s: no saved %.3g%% checkpoint.',runDirs{j},target);
            continue;
        end
        s = snapshots{k};
        if ~isfield(s,'completed')
            error('NetworkGrowth:NoDistributions','Distribution payload is missing in %s.',runDirs{j});
        end
        switch population
            case 'completed'
                n = s.completed.kuhnCount; contour = s.completed.contourLength;
            case 'dangling'
                n = s.dangling.kuhnCount; contour = s.dangling.contourLength;
            case 'all'
                n = [s.completed.kuhnCount;s.dangling.kuhnCount];
                contour = [s.completed.contourLength;s.dangling.contourLength];
        end
        samples(end+1,:) = {n,contour,s.completed.endToEnd}; %#ok<AGROW>
        legendLabels{end+1} = sprintf('%s: %g%% (actual %.3g%%)',labels{j},target,s.usagePercent); %#ok<AGROW>
        row = table(j,string(runDirs{j}),string(index.file(k)),target, ...
            s.usagePercent,s.event,s.ssaTime,numel(n),numel(s.completed.endToEnd), ...
            'VariableNames',{'runIndex','runDir','file','targetPercent', ...
            'actualPercent','event','ssaTime','contourSampleCount','endToEndSampleCount'});
        selection = [selection;row]; %#ok<AGROW>
    end
end
if isempty(selection), error('NetworkGrowth:NoMilestones','No requested checkpoints were saved.'); end
edges = {options.nEdges,options.contourEdges,options.rEdges};
for col = 1:3
    if isempty(edges{col})
        values = vertcat(samples{:,col});
        if col == 1
            if isempty(values), values = 1; end
            lo = floor(min(values))-.5; hi = ceil(max(values))+.5;
            edges{col} = linspace(lo,hi,min(100,round(hi-lo))+1);
        else
            hi = max([vertcat(values);1]);
            edges{col} = linspace(0,hi,31);
        end
    end
    validateattributes(edges{col},{'numeric'},{'vector','real','finite','increasing'});
    if numel(edges{col}) < 2, error('NetworkGrowth:PlotBins','At least two bin edges are required.'); end
end
fig = figure('Color','w','Visible',options.visible,'Position',[100 100 1400 470]);
layout = tiledlayout(fig,1,3,'TileSpacing','compact');
colors = lines(size(samples,1));
xlabels = {'n (Kuhn segments)','contour length n b','end-to-end distance r'};
titles = {[population ' chains'],[population ' chains'],'completed chains only'};
for col = 1:3
    ax = nexttile(layout); hold(ax,'on');
    for k = 1:size(samples,1)
        v = samples{k,col};
        if isempty(v)
            plot(ax,NaN,NaN,'Color',colors(k,:),'DisplayName',[legendLabels{k} ' (empty)']);
        else
            % Normalize by all samples, even when user-specified bins clip tails.
            counts = histcounts(v,edges{col}) / numel(v);
            stairs(ax,edges{col},[counts counts(end)],'LineWidth',1.4, ...
                'Color',colors(k,:),'DisplayName',legendLabels{k});
        end
    end
    xlabel(ax,xlabels{col}); ylabel(ax,'probability per bin');
    title(ax,titles{col}); grid(ax,'on');
    if col == 3, legend(ax,'show','Location','eastoutside','Interpreter','none'); end
end
title(layout,'Chain distributions at free-crosslinker usage checkpoints');
if ~isempty(options.outputFile), exportgraphics(fig,options.outputFile,'Resolution',180); end
end
