% Three single-run or replicate-average figures: Kuhn number fraction, Kuhn weight fraction,
% and end-to-end distance. Usage snapshots run left to right in each.
% Each figure shares x/y limits across its own snapshots only.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(projectRoot,'SRC'));
% Supply one run folder, or multiple replicates with identical physics and distinct seeds.
% A single run has no replicate SD shading; prescribed mixtures use no Flory fits.
runFolders = {}; % supply saved run folders, or select a batch below
% Alternatively select a batch and xA; leave runFolders empty to use these:
batchFolder = ''; % contains batch_summary.mat from run_sweep_preformed
selectedXA = .025;
selectedParameterSet = 1; % used when the sweep has multiple parameter sets

averageOptions = struct();
averageOptions.usagePercent = [20 40 60 80 100];
averageOptions.nBinWidth = 1;       % integer Kuhn bin width; larger = coarser histogram
averageOptions.distanceBins = 30;
averageOptions.kuhnYLim = [];   % fixed number-fraction y range; [] = automatic
averageOptions.weightYLim = [];    % weight figure's own shared range; independent of number figure
averageOptions.showSpread = true;  % shaded +/- one replicate SD, not confidence intervals
averageOptions.requireAllReplicates = true; % avoid changing replicate cohorts at later usage
averageOptions.visible = 'on';
averageOptions.fitFlory = false; % prescribed chain lengths: no Flory growth fit
averageOptions.fitFloryWeight = false;
averageOptions.showWeightFromNumberFit = false;

if isempty(runFolders) && ~isempty(batchFolder)
    data = load(fullfile(batchFolder,'batch_summary.mat'),'summary');
    rows = abs(data.summary.xA-selectedXA)<1e-10 & data.summary.status=="completed";
    if ismember('parameterSetIndex',data.summary.Properties.VariableNames)
        rows = rows & data.summary.parameterSetIndex==selectedParameterSet;
    end
    runFolders = cellstr(data.summary.runDir(rows));
end
runFolders = cellstr(runFolders);
if isempty(runFolders)
    error('NetworkGrowth:ReplicateFolders','Set at least one run folder, or batchFolder and selectedXA.');
end
outputRoot = fullfile(projectRoot,'replicate_averages','preformed_chains');
[averageOutputDir,~] = createUniqueRunFolder(outputRoot,'average');
% Base filename: exports append _kuhn_number, _kuhn_weight, and _end_to_end.
averageOptions.outputFile = fullfile(averageOutputDir,'replicate_average.png');
[averageFigure,averageData] = averageNetworkDumps(runFolders,averageOptions);
save(fullfile(averageOutputDir,'replicate_average.mat'),'averageData','averageOptions','runFolders');
writetable(averageData.selection,fullfile(averageOutputDir,'replicate_selection.csv'));
fprintf('Distribution plots saved to: %s\n',averageOutputDir);



