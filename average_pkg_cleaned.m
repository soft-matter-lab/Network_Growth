% Replicate-averaged evolution plots. Edit folders and options, then run.
addpath(fullfile(fileparts(mfilename('fullpath')),'SRC'));
% Supply all replicate run folders with identical physics/composition and distinct seeds.
runFolders = {
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_125518_906\run_001_rep_01_xA_0.05_seed_11231_20261026_125518_950'
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_125518_906\run_002_rep_02_xA_0.05_seed_2015_20261026_125542_837'
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_125518_906\run_003_rep_03_xA_0.05_seed_3215464_20261026_125612_807'
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_125518_906\run_004_rep_04_xA_0.05_seed_5_20261026_125640_519'
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_125518_906\run_005_rep_05_xA_0.05_seed_10_20261026_125711_143'
};
% Alternatively select a batch and xA; leave runFolders empty to use these:
batchFolder = ''; % contains batch_summary.mat from run_pkg_cleaned_sweeps
selectedXA = .05;

averageOptions = struct();
averageOptions.usagePercent = [25 50 100];
averageOptions.nBinWidth = 1;       % integer Kuhn bin width; larger = coarser histogram
averageOptions.distanceBins = 30;
averageOptions.showSpread = true;  % shaded +/- one replicate SD, not confidence intervals
averageOptions.requireAllReplicates = true; % avoid changing replicate cohorts at later usage
averageOptions.visible = 'on';

if isempty(runFolders) && ~isempty(batchFolder)
    data = load(fullfile(batchFolder,'batch_summary.mat'),'summary');
    rows = abs(data.summary.xA-selectedXA)<1e-10 & data.summary.status=="completed";
    runFolders = cellstr(data.summary.runDir(rows));
end
if numel(runFolders)<2
    error('NetworkGrowth:ReplicateFolders','Set at least two runFolders, or batchFolder and selectedXA.');
end
outputRoot = fullfile(fileparts(mfilename('fullpath')),'replicate_averages');
[averageOutputDir,~] = createUniqueRunFolder(outputRoot,'average');
averageOptions.outputFile = fullfile(averageOutputDir,'replicate_average.png');
[averageFigure,averageData] = averageNetworkDumps(runFolders,averageOptions);
save(fullfile(averageOutputDir,'replicate_average.mat'),'averageData','averageOptions','runFolders');
writetable(averageData.selection,fullfile(averageOutputDir,'replicate_selection.csv'));
fprintf('Averaged distributions saved to: %s\n',averageOutputDir);
