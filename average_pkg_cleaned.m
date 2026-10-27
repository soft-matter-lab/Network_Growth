% Three replicate-average figures: Kuhn number fraction, Kuhn weight fraction,
% and end-to-end distance. Usage snapshots run left to right in each.
% Each figure shares x/y limits across its own snapshots only.
addpath(fullfile(fileparts(mfilename('fullpath')),'SRC'));
% Supply all replicate run folders with identical physics/composition and distinct seeds.
runFolders = {
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_164546_823\run_001_rep_01_xA_0.05_seed_11231_20261026_164548_463'
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_164546_823\run_002_rep_02_xA_0.05_seed_2015_20261026_170757_626'
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_164546_823\run_003_rep_03_xA_0.05_seed_3215464_20261026_173121_381'
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_164546_823\run_004_rep_04_xA_0.05_seed_5_20261026_175125_518'
    'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_164546_823\run_005_rep_05_xA_0.05_seed_10_20261026_181604_459'
};
% Alternatively select a batch and xA; leave runFolders empty to use these:
batchFolder = ''; % contains batch_summary.mat from run_pkg_cleaned_sweeps
selectedXA = .05;

averageOptions = struct();
averageOptions.usagePercent = [5 25 50 75 100];
averageOptions.nBinWidth = 1;       % integer Kuhn bin width; larger = coarser histogram
averageOptions.distanceBins = 30;
averageOptions.kuhnYLim = [0 .35];   % fixed number-fraction y range; [] = automatic
averageOptions.weightYLim = [];    % weight figure's own shared range; independent of number figure
averageOptions.showSpread = true;  % shaded +/- one replicate SD, not confidence intervals
averageOptions.requireAllReplicates = true; % avoid changing replicate cohorts at later usage
averageOptions.visible = 'on';
averageOptions.fitFlory = true;    % overlay a fit on each Kuhn number-fraction panel
averageOptions.fitFloryWeight = true; % independent fit to each weight-fraction panel
averageOptions.showWeightFromNumberFit = true; % also predict weight using number-fit p (no refit)
averageOptions.floryModel = 'shifted'; % 'shifted' or 'original'
averageOptions.floryMinCount = 4;  % shifted: n_min in (1-p)*p^(n-n_min), n>=n_min
% original: (1-p)*p^(n-1), n>=1; floryMinCount only sets the fitting cutoff.
% Shifted minimum must not exceed any observed count; neither model is extrapolated.
averageOptions.floryOriginalAnchorWeight = 1; % original only: 1=ordinary LS; 10=soft anchor at cutoff
averageOptions.floryBootstrapSamples = 2000; % resample whole replicates for 95% CI
averageOptions.floryBootstrapSeed = 1729;    % reproducible; does not alter global RNG

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
% Base filename: exports append _kuhn_number, _kuhn_weight, and _end_to_end.
averageOptions.outputFile = fullfile(averageOutputDir,'replicate_average.png');
[averageFigure,averageData] = averageNetworkDumps(runFolders,averageOptions);
save(fullfile(averageOutputDir,'replicate_average.mat'),'averageData','averageOptions','runFolders');
writetable(averageData.selection,fullfile(averageOutputDir,'replicate_selection.csv'));
if averageOptions.fitFlory
    writetable(averageData.floryFitTable,fullfile(averageOutputDir,'flory_fits.csv'));
    fprintf('Flory model: %s; normalized on n >= %g; fit uses n >= %g.\n', ...
        averageData.floryFit{1}.modelStyle,averageData.floryFit{1}.supportMin,averageData.floryFit{1}.minCount);
    disp(averageData.floryFitTable(:,{'targetPercent','model','p','pCI95Low','pCI95High','fullCdfDistance'}));
end
if averageOptions.fitFloryWeight
    writetable(averageData.floryWeightFitTable,fullfile(averageOutputDir,'flory_weight_fits.csv'));
    fprintf('Weight fits: W(n)=n*P(n)/mean(n); pDifference = p_weight - p_number.\n');
    disp(averageData.floryWeightFitTable(:,{'targetPercent','model','p','pCI95Low','pCI95High', ...
        'pFromNumber','pDifference','pDifferenceCI95Low','pDifferenceCI95High'}));
end
fprintf('Averaged distributions saved to: %s\n',averageOutputDir);
