function [fig, result] = averageNetworkDumps(runDirs, options)
% Equal-replicate mean probability histograms of completed chains at milestones.
% Normalize each replicate separately before averaging; never pool chain counts.
% Returns three figures (number fraction, weight fraction, end-to-end distance),
% each with one column per usage target. outputFile is a base name; exports add suffixes.
if nargin < 2, options = struct(); end
options = mergeConfig(struct('usagePercent',[25 50 100],'nBinWidth',1, ...
    'distanceBins',30,'showSpread',true,'requireAllReplicates',true, ...
    'visible','on','outputFile','','kuhnYLim',[],'weightYLim',[], ...
    'fitFlory',true,'fitFloryWeight',[],'showWeightFromNumberFit',true,'floryModel','shifted', ...
    'floryMinCount',0,'floryOriginalAnchorWeight',1, ...
    'floryBootstrapSamples',2000,'floryBootstrapSeed',1729),options);
if isempty(options.fitFloryWeight), options.fitFloryWeight=options.fitFlory; end
validateattributes(options.usagePercent,{'numeric'},{'vector','nonempty','finite','>=',0,'<=',100});
validateattributes(options.nBinWidth,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(options.distanceBins,{'numeric'},{'scalar','integer','positive','finite'});
if ~isempty(options.kuhnYLim)
    validateattributes(options.kuhnYLim,{'numeric'},{'vector','numel',2,'real','finite','nonnegative','increasing'});
end
if ~isempty(options.weightYLim)
    validateattributes(options.weightYLim,{'numeric'},{'vector','numel',2,'real','finite','nonnegative','increasing'});
end
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
        reference = physics;
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
    'distributionNames',{{'kuhn_number','kuhn_weight','end_to_end'}}, ...
    'edges',{{nEdges,nEdges,rEdges}},'mean',{{}},'sd',{{}},'perReplicate',{{}}, ...
    'replicateCount',zeros(1,T),'weighting','equal replicate; probability per bin', ...
    'distributionReplicateCount',zeros(3,T),'outputFiles',{{}}, ...
    'spread','one sample standard deviation across replicate histograms');
for col = 1:3
    edges = result.edges{col}; values = nan(R,numel(edges)-1,T);
    for t = 1:T
        for rep = 1:R
            if col<3, v=samples{rep,t,1}; else, v=samples{rep,t,2}; end
            if ~isempty(v)
                if col==2
                    [fraction,totalSegments] = chainWeightFractions(v,edges);
                    if totalSegments>0
                        values(rep,:,t)=fraction;
                    elseif options.requireAllReplicates
                        error('NetworkGrowth:EmptyWeightPopulation', ...
                            'Run %d at %g%% has no completed-chain Kuhn segments.',rep,targets(t));
                    else
                        warning('NetworkGrowth:EmptyWeightPopulation', ...
                            'Excluding run %d at %g%% from the weight fraction: zero total Kuhn segments.',rep,targets(t));
                    end
                else
                    values(rep,:,t)=histcounts(v,edges)/numel(v);
                end
            end
        end
        valid = ~isnan(values(:,1,t)); nr = sum(valid);
        if col==1, result.replicateCount(t) = nr; end
        result.distributionReplicateCount(col,t) = nr;
        result.mean{col}(t,:) = mean(values(valid,:,t),1);
        if nr>=2
            result.sd{col}(t,:) = std(values(valid,:,t),0,1);
        else
            result.sd{col}(t,:) = nan(1,numel(edges)-1);
        end
    end
    result.perReplicate{col} = values;
end
result.floryFit = cell(1,T); result.floryFitTable = table();
if options.fitFlory
    fitOptions=struct('model',options.floryModel,'minCount',options.floryMinCount, ...
        'originalAnchorWeight',options.floryOriginalAnchorWeight, ...
        'bootstrapSamples',options.floryBootstrapSamples,'bootstrapSeed',options.floryBootstrapSeed);
    for t=1:T
        f=fitFloryNumberDistribution(samples(:,t,1),nEdges,fitOptions);
        if f.replicateCount<result.replicateCount(t)
            if options.requireAllReplicates
                error('NetworkGrowth:EmptyFloryTail', ...
                    'At %g%%, some replicates have no chains with n >= %g. Lower floryMinCount or set requireAllReplicates=false.',targets(t),options.floryMinCount);
            else
                warning('NetworkGrowth:EmptyFloryTail', ...
                    'At %g%%, fitting uses %d of %d nonempty replicates with n >= %g.', ...
                    targets(t),f.replicateCount,result.replicateCount(t),options.floryMinCount);
            end
        end
        result.floryFit{t}=f;
        row=table(targets(t),f.p,f.pCI95(1),f.pCI95(2),f.meanCount,f.cdfDistance, ...
            f.fitCdfDistance,f.binRMSE,f.outsideBinProbability,f.probabilityBelowCutoff, ...
            f.minCount,f.replicateCount,f.chainCount,f.totalChainCount,f.observedTailFraction, ...
            string(f.modelStyle),f.supportMin,f.originalAnchorWeight, ...
            'VariableNames',{'targetPercent','p','pCI95Low','pCI95High','meanFittedKuhnCount', ...
            'fullCdfDistance','tailCdfDistance','fullBinRMSE','outsideBinProbability', ...
            'probabilityBelowCutoff','fitMinCount','fitReplicateCount','fitChainCount','totalChainCount', ...
            'observedTailFraction','model','modelMinCount','originalAnchorWeight'});
        result.floryFitTable=[result.floryFitTable;row]; %#ok<AGROW>
    end
end
result.floryWeightFit=cell(1,T); result.floryWeightFitTable=table();
if options.fitFloryWeight
    weightOptions=struct('model',options.floryModel,'minCount',options.floryMinCount, ...
        'bootstrapSamples',options.floryBootstrapSamples,'bootstrapSeed',options.floryBootstrapSeed);
    for t=1:T
        weightOptions.numberFit=result.floryFit{t};
        f=fitFloryWeightDistribution(samples(:,t,1),nEdges,weightOptions);
        if f.replicateCount<result.distributionReplicateCount(2,t)
            if options.requireAllReplicates
                error('NetworkGrowth:EmptyFloryWeightTail', ...
                    'At %g%%, some replicates have no Kuhn weight at n >= %g.',targets(t),f.minCount);
            else
                warning('NetworkGrowth:EmptyFloryWeightTail', ...
                    'At %g%%, weight fitting uses %d of %d weight-distribution replicates.', ...
                    targets(t),f.replicateCount,result.distributionReplicateCount(2,t));
            end
        end
        result.floryWeightFit{t}=f;
        row=table(targets(t),string(f.modelStyle),f.numberSupportMin,f.supportMin,f.minCount, ...
            f.p,f.pCI95(1),f.pCI95(2),f.pFromNumber,f.pDifference, ...
            f.pDifferenceCI95(1),f.pDifferenceCI95(2),f.cdfDistance,f.binRMSE, ...
            f.tailSSE,f.numberPredictionTailSSE,f.outsideBinProbability, ...
            f.probabilityBelowCutoff,f.replicateCount,f.chainCount,f.totalKuhnCount, ...
            'VariableNames',{'targetPercent','model','modelMinCount','weightSupportMin','fitMinCount', ...
            'p','pCI95Low','pCI95High','pFromNumber','pDifference','pDifferenceCI95Low', ...
            'pDifferenceCI95High','fullCdfDistance','fullBinRMSE','tailSSE', ...
            'numberPredictionTailSSE','outsideBinProbability','probabilityBelowCutoff', ...
            'fitReplicateCount','fitChainCount','totalKuhnCount'});
        result.floryWeightFitTable=[result.floryWeightFitTable;row]; %#ok<AGROW>
    end
end
fig = gobjects(1,3); colors = [0 .38 .68; .78 .29 .12; .15 .50 .43];
labels = {'n (Kuhn segments)','n (Kuhn segments)','end-to-end distance r'};
names = {'Kuhn-count number fraction','Kuhn-count weight fraction','End-to-end distance'};
ylabels = {'Chain number fraction','Kuhn-segment weight fraction','Chain number fraction'};
curveX=cell(2,T); curveY=cell(2,T);
predictionX=cell(1,T); predictionY=cell(1,T);
if options.fitFlory
    for t=1:T
        [curveX{1,t},curveY{1,t}]=floryDisplayCurve(result.floryFit{t},result.edges{1});
    end
end
if options.fitFloryWeight
    for t=1:T
        f=result.floryWeightFit{t};
        [curveX{2,t},curveY{2,t}]=floryWeightDisplayCurve(f,result.edges{2});
        if options.showWeightFromNumberFit && isfinite(f.pFromNumber)
            [predictionX{t},predictionY{t}]=floryWeightDisplayCurve(f,result.edges{2},f.pFromNumber);
        end
    end
end
fitEnabled=[options.fitFlory options.fitFloryWeight false];
fits={result.floryFit,result.floryWeightFit};
for col = 1:3
    fig(col) = figure('Color','w','Visible',options.visible,'Name',names{col}, ...
        'Position',[60+20*(col-1) 80 min(1800,180+420*T) 540]);
    layout = tiledlayout(fig(col),1,T,'TileSpacing','compact','Padding','compact');
    edges = result.edges{col};
    x = reshape([edges(1:end-1);edges(2:end)],1,[]);
    % Include the displayed SD bands when setting this figure's shared limits.
    upper = result.mean{col};
    if options.showSpread
        spread = result.sd{col}; spread(~isfinite(spread)) = 0;
        upper = min(1,upper+spread);
    end
    upper = upper(:);
    if fitEnabled(col)
        for t=1:T, upper=[upper;curveY{col,t}(:)]; end %#ok<AGROW>
    end
    if col==2
        for t=1:T, upper=[upper;predictionY{t}(:)]; end %#ok<AGROW>
    end
    peak = max([0;upper(isfinite(upper))]);
    sharedYLim = [0 max(.05,1.1*peak)];
    if col==1 && ~isempty(options.kuhnYLim), sharedYLim=options.kuhnYLim(:)'; end
    if col==2 && ~isempty(options.weightYLim), sharedYLim=options.weightYLim(:)'; end
    legendHandles=gobjects(0); legendAxes=[];
    for t = 1:T
        ax = nexttile(layout); hold(ax,'on');
        set(ax,'FontName','Arial','FontSize',11,'LineWidth',.8,'Box','off', ...
            'TickDir','out','XGrid','off','YGrid','on','GridAlpha',.1,'Layer','top');
        m = repelem(result.mean{col}(t,:),2); sd = repelem(result.sd{col}(t,:),2);
        nr = result.distributionReplicateCount(col,t);
        hSpread=gobjects(0); hFit=gobjects(0); hPrediction=gobjects(0);
        if options.showSpread && nr>=2
            hSpread=fill(ax,[x fliplr(x)],[max(0,m-sd) fliplr(min(1,m+sd))],colors(col,:), ...
                'FaceAlpha',.1,'EdgeColor','none','DisplayName','+/- 1 SD');
        end
        hMean=plot(ax,x,m,'Color',colors(col,:),'LineWidth',2,'DisplayName','Replicate mean');
        title(ax,sprintf('%g%% usage',targets(t)),'FontSize',14,'FontWeight','bold');
        if fitEnabled(col) && isfinite(fits{col}{t}.p)
            f=fits{col}{t};
            if strcmp(f.modelStyle,'shifted')
                minimum=f.supportMin;
                if col==2, minimum=f.numberSupportMin; end
                fitLabel=sprintf('Shifted Flory (n_{min} = %g)',minimum);
            else
                fitLabel=sprintf('Original Flory (fit n >= %g)',f.minCount);
            end
            if col==2, fitLabel=[fitLabel ' weight fit']; end
            if col==2 && ~isempty(predictionX{t})
                hPrediction=plot(ax,predictionX{t},predictionY{t},'-.', ...
                    'Color',[.55 .55 .55],'LineWidth',1.6, ...
                    'DisplayName','Prediction using number-fit p','Tag','weightFromNumberCurve');
                if f.pFromNumber==0, set(hPrediction,'Marker','o','MarkerSize',4); end
            end
            hFit=plot(ax,curveX{col,t},curveY{col,t},'k--','LineWidth',2, ...
                'DisplayName',fitLabel,'Tag','floryFitCurve');
            if f.p==0, set(hFit,'Marker','o','MarkerSize',4); end
            if col==2 && ~isempty(predictionX{t})
                cohortLabel=sprintf('%d replicates',f.replicateCount);
                if f.replicateCount~=result.floryFit{t}.replicateCount
                    cohortLabel=sprintf('%d weight / %d number replicates',f.replicateCount,result.floryFit{t}.replicateCount);
                end
                subtitle(ax,{sprintf('p_W = %.4f | p_N = %.4f',f.p,f.pFromNumber),cohortLabel},'FontSize',11);
            elseif col==2
                subtitle(ax,sprintf('p_W = %.4f | %d replicates',f.p,f.replicateCount),'FontSize',11);
            else
                subtitle(ax,sprintf('p = %.4f  |  %d replicates',f.p,f.replicateCount),'FontSize',11);
            end
        elseif fitEnabled(col) && nr>0
            subtitle(ax,'No chains at or above fitting cutoff');
        else
            subtitle(ax,sprintf('%d replicates',nr),'FontSize',11);
        end
        xlim(ax,[edges(1) edges(end)]); ylim(ax,sharedYLim);
        handles=[hMean hSpread hFit hPrediction];
        if numel(handles)>numel(legendHandles), legendHandles=handles; legendAxes=ax; end
        if nr==0
            text(ax,.5,.5,'No contributing replicates','Units','normalized','HorizontalAlignment','center');
        end
    end
    title(layout,names{col},'FontName','Arial','FontSize',17,'FontWeight','bold');
    xlabel(layout,labels{col},'FontName','Arial','FontSize',12);
    ylabel(layout,ylabels{col},'FontName','Arial','FontSize',12);
    lgd=legend(legendAxes,legendHandles,'Orientation','horizontal','Box','off','FontSize',10);
    lgd.Layout.Tile='south';
end
if ~isempty(options.outputFile)
    [folder,base,ext] = fileparts(options.outputFile);
    if isempty(ext), ext='.png'; end
    for col = 1:3
        result.outputFiles{col} = fullfile(folder,[base '_' result.distributionNames{col} ext]);
        % Hover controls can otherwise be captured when figures are visible.
        axesList=findall(fig(col),'Type','axes');
        for ax=axesList', ax.Toolbar.Visible='off'; end
        exportgraphics(fig(col),result.outputFiles{col},'Resolution',180);
    end
end
end

function physics = physicsConfig(cfg)
% Compare every non-output setting, including chemistry and relaxation settings.
names = fieldnames(cfg);
output = startsWith(names,{'dump','save','show','video','output','data','plot','eqPlot','heatGrid'});
output = output | ismember(names,{'randSeed','makeVideo'});
physics = orderfields(rmfield(cfg,names(output)));
end
