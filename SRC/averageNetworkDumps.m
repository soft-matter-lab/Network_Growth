function [fig, result] = averageNetworkDumps(runDirs, options)
% Equal-replicate mean probability histograms of completed chains at milestones.
% Normalize each replicate separately before averaging; never pool chain counts.
if nargin < 2, options = struct(); end
options = mergeConfig(struct('usagePercent',[25 50 100],'nBinWidth',1, ...
    'distanceBins',30,'showSpread',true,'requireAllReplicates',true, ...
    'visible','on','outputFile',''),options);
validateattributes(options.usagePercent,{'numeric'},{'vector','nonempty','finite','>=',0,'<=',100});
validateattributes(options.nBinWidth,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(options.distanceBins,{'numeric'},{'scalar','integer','positive','finite'});
runDirs = cellstr(runDirs);
R = numel(runDirs); targets = unique(options.usagePercent(:)'); T = numel(targets);
if R<2, error('NetworkGrowth:ReplicateFolders','At least two replicate folders are required.'); end
samples = cell(R,T,2); selection = table(); seeds = zeros(R,1);
maxN = 1; maxD = 0;
for rep = 1:R
    [states,metadata,index] = loadNetworkDumps(runDirs{rep},targets);
    cfg = metadata.requestedConfig;
    validateattributes(cfg.randSeed,{'numeric'},{'scalar','integer','nonnegative','finite'});
    seeds(rep) = cfg.randSeed;
    physics = physicsConfig(cfg);
    if rep == 1
        reference = physics; b = cfg.b;
    elseif ~isequaln(reference,physics)
        error('NetworkGrowth:ReplicateSettings','Run %d has different settings; select one composition and parameter set.',rep);
    end
    for t = 1:T
        k = find(cellfun(@(s)any(abs(s.usageTargets-targets(t))<1e-9),states),1);
        actual = NaN; count = 0; file = ""; status = "missing checkpoint";
        if ~isempty(k)
            s = states{k}; actual = s.usagePercent; file = index.file(k);
            if ~isfield(s,'completed')
                error('NetworkGrowth:NoDistributions','Run %d lacks distribution data.',rep);
            end
            count = numel(s.completed.kuhnCount);
            status = "empty population";
            if count>0
                samples{rep,t,1} = s.completed.kuhnCount(:);
                samples{rep,t,2} = s.completed.endToEnd(:);
                maxN = max([maxN;samples{rep,t,1}]);
                maxD = max([maxD;samples{rep,t,2}]);
                status = "included";
            end
        end
        if status~="included" && options.requireAllReplicates
            error('NetworkGrowth:IncompleteReplicates', ...
                'Run %d at %g%%: %s. Use common reached milestones, or explicitly set requireAllReplicates=false.',rep,targets(t),status);
        end
        row = table(rep,string(runDirs{rep}),seeds(rep),targets(t),actual,count,file,status, ...
            'VariableNames',{'replicate','runDir','seed','targetPercent','actualPercent','chainCount','file','status'});
        selection = [selection;row]; %#ok<AGROW>
    end
end
if numel(unique(seeds))~=R
    error('NetworkGrowth:DuplicateSeeds','Replicates must have distinct seeds; check for duplicate folders or repeated realizations.');
end
if any(selection.status~="included")
    warning('NetworkGrowth:IncompleteReplicates','Excluded missing/empty populations; contributing replicate counts are labeled and saved.');
end
nEdges = -.5:options.nBinWidth:(ceil((maxN+.5)/options.nBinWidth)*options.nBinWidth-.5);
if nEdges(end)<maxN, nEdges(end+1)=nEdges(end)+options.nBinWidth; end
if maxD<=0, maxD=1; end
rEdges = linspace(0,maxD*1.05,options.distanceBins+1);
result = struct('targets',targets,'seeds',seeds,'selection',selection, ...
    'edges',{{nEdges,nEdges*b,rEdges}},'mean',{{}},'sd',{{}},'perReplicate',{{}}, ...
    'replicateCount',zeros(1,T),'weighting','equal replicate; probability per bin', ...
    'spread','one sample standard deviation across replicate histograms');
for col = 1:3
    edges = result.edges{col}; values = nan(R,numel(edges)-1,T);
    for t = 1:T
        for rep = 1:R
            if col<3, v=samples{rep,t,1}; else, v=samples{rep,t,2}; end
            if col==2, v=v*b; end
            if ~isempty(v), values(rep,:,t)=histcounts(v,edges)/numel(v); end
        end
        valid = ~isnan(values(:,1,t)); nr = sum(valid);
        result.replicateCount(t) = nr;
        result.mean{col}(t,:) = mean(values(valid,:,t),1);
        if nr>=2
            result.sd{col}(t,:) = std(values(valid,:,t),0,1);
        else
            result.sd{col}(t,:) = nan(1,numel(edges)-1);
        end
    end
    result.perReplicate{col} = values;
end
fig = figure('Color','w','Visible',options.visible,'Position',[100 100 1400 480]);
layout = tiledlayout(fig,1,3,'TileSpacing','compact'); colors = lines(T);
labels = {'n (Kuhn segments)','contour length n b','end-to-end distance r'};
for col = 1:3
    ax = nexttile(layout); hold(ax,'on'); edges = result.edges{col};
    x = reshape([edges(1:end-1);edges(2:end)],1,[]);
    for t = 1:T
        m = repelem(result.mean{col}(t,:),2); sd = repelem(result.sd{col}(t,:),2);
        if options.showSpread && result.replicateCount(t)>=2
            fill(ax,[x fliplr(x)],[max(0,m-sd) fliplr(min(1,m+sd))],colors(t,:), ...
                'FaceAlpha',.12,'EdgeColor','none','HandleVisibility','off');
        end
        plot(ax,x,m,'Color',colors(t,:),'LineWidth',1.8, ...
            'DisplayName',sprintf('%g%% (%d replicates)',targets(t),result.replicateCount(t)));
    end
    xlabel(ax,labels{col}); ylabel(ax,'mean probability per bin'); grid(ax,'on');
    xlim(ax,[edges(1) edges(end)]);
    if col==3, legend(ax,'show','Location','eastoutside'); end
end
title(layout,'Completed-chain distributions: equal-replicate means');
if ~isempty(options.outputFile), exportgraphics(fig,options.outputFile,'Resolution',180); end
end

function physics = physicsConfig(cfg)
% Compare every non-output setting, including chemistry and relaxation settings.
names = fieldnames(cfg);
output = startsWith(names,{'dump','save','show','video','output','data','plot','eqPlot','heatGrid'});
output = output | ismember(names,{'randSeed','makeVideo'});
physics = orderfields(rmfield(cfg,names(output)));
end
