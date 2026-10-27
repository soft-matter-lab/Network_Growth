function report = postprocessNetworkDumps(runOrDumpDir, options)
% Export a saved run's video, individual checkpoints and distribution comparison.
% User-facing entry point: visualize_pkg_cleaned.m in the project root.
if nargin < 2, options = struct(); end
defaults = struct('outputRoot','','makeVideo',true,'videoFrameEvery',1, ...
    'videoFPS',6,'videoHoldLastSeconds',2,'videoPreview',false,'heatGrid',150, ...
    'saveCheckpointPlots',true,'checkpointUsagePercent',0:5:100,'saveFinalPlot',true, ...
    'saveComparison',true,'comparisonUsagePercent',[0 25 50 75 100], ...
    'population','completed','showComparisonFigure',true, ...
    'kuhnXLim',[],'endToEndXLim',[],'histogramYLim',[], ...
    'kuhnYLim',[],'endToEndYLim',[],'endToEndBins',30, ...
    'displacementLimit',[],'displacementBins',40,'displacementCLim',[], ...
    'displacementColorPercentile',98, ...
    'weightXLim',[],'weightYLim',[]);
options = mergeConfig(defaults,options);
validateattributes(options.checkpointUsagePercent,{'numeric'},{'vector','nonempty','finite','>=',0,'<=',100});
validateattributes(options.comparisonUsagePercent,{'numeric'},{'vector','nonempty','finite','>=',0,'<=',100});
[files,metadata,dumpDir] = listNetworkDumpFiles(runOrDumpDir);
if isempty(files), error('NetworkGrowth:NoSnapshots','There are no completed snapshot files.'); end
first = load(fullfile(dumpDir,files(1).name),'snapshot');
hasNetwork = isfield(first.snapshot,'nodes');
hasDistributions = isfield(first.snapshot,'completed');
outputRoot = char(options.outputRoot);
if isempty(outputRoot), outputRoot = fullfile(fileparts(dumpDir),'postprocessing'); end
[outputDir,~] = createUniqueRunFolder(outputRoot,'visualization');
report = struct('outputDir',outputDir,'dumpDir',dumpDir,'videoFile','', ...
    'comparisonFile','','checkpointFiles',strings(0,1),'finalFile','', ...
    'histogramStyle',[]);
if options.saveCheckpointPlots || options.saveFinalPlot
    % Include every stored state and final relaxation when fixing bins/axes.
    report.histogramStyle = networkDumpHistogramStyle(dumpDir,options);
end
save(fullfile(outputDir,'visualization_config.mat'),'options','dumpDir','metadata');

if options.makeVideo
    if hasNetwork
        report.videoFile = fullfile(outputDir,'network_evolution.mp4');
        frameIndex = exportNetworkDumpVideo(dumpDir,report.videoFile,options);
        writetable(frameIndex,fullfile(outputDir,'video_frames.csv'));
    else
        warning('NetworkGrowth:NoNetworkVideo','Network payload was disabled; skipping the video.');
    end
end

if options.saveCheckpointPlots
    [snapshots,~,index] = loadNetworkDumps(dumpDir,options.checkpointUsagePercent);
    targets = cellfun(@(s)s.usageTargets,snapshots,'UniformOutput',false);
    missing = setdiff(options.checkpointUsagePercent,[targets{:}]);
    if ~isempty(missing)
        warning('NetworkGrowth:MissingMilestone','No saved checkpoints for usage %s%%; skipped.',mat2str(missing));
    end
    checkpointDir = fullfile(outputDir,'checkpoints'); mkdir(checkpointDir);
    for k = 1:numel(snapshots)
        [~,stem] = fileparts(char(index.file(k)));
        outputFile = fullfile(checkpointDir,[stem '.png']);
        fig = plotNetworkDump(snapshots{k},struct('visible','off','outputFile',outputFile, ...
            'histogramStyle',report.histogramStyle));
        close(fig);
        report.checkpointFiles(end+1,1) = string(outputFile);
    end
    if ~isempty(index)
        index.usageTargets = string(cellfun(@mat2str,index.usageTargets,'UniformOutput',false));
        writetable(index,fullfile(outputDir,'checkpoint_index.csv'));
    end
    fprintf('Saved %d checkpoint figures.\n',numel(snapshots));
    clear snapshots;
end

if options.saveFinalPlot
    final = load(fullfile(dumpDir,files(end).name),'snapshot');
    if strcmp(final.snapshot.phase,'final')
        report.finalFile = fullfile(outputDir,'final_state.png');
        fig = plotNetworkDump(final.snapshot,struct('visible','off','outputFile',report.finalFile, ...
            'histogramStyle',report.histogramStyle));
        close(fig);
    else
        warning('NetworkGrowth:NoFinalSnapshot','No final snapshot was saved; final_state.png is skipped.');
    end
end

if options.saveComparison
    if hasDistributions
        % Skip absent targets without aborting an otherwise successful video export.
        snapshots = loadNetworkDumps(dumpDir,options.comparisonUsagePercent);
        targets = cellfun(@(s)s.usageTargets,snapshots,'UniformOutput',false);
        available = intersect(options.comparisonUsagePercent,[targets{:}]);
        missing = setdiff(options.comparisonUsagePercent,available);
        if ~isempty(missing)
            warning('NetworkGrowth:MissingMilestone','Comparison has no saved checkpoints for usage %s%%; skipped.',mat2str(missing));
        end
        clear snapshots;
        if ~isempty(available)
            visible = 'off'; if options.showComparisonFigure, visible = 'on'; end
            report.comparisonFile = fullfile(outputDir,'distribution_comparison.png');
            [fig,selection] = plotNetworkDumps(dumpDir,struct('usagePercent',available, ...
                'population',options.population,'visible',visible,'outputFile',report.comparisonFile));
            if ~options.showComparisonFigure, close(fig); end
            writetable(selection,fullfile(outputDir,'comparison_index.csv'));
        end
    else
        warning('NetworkGrowth:NoDistributions','Distribution payload was disabled; skipping the comparison.');
    end
end
save(fullfile(outputDir,'visualization_report.mat'),'report');
end
