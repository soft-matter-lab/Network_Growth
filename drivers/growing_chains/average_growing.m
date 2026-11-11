% Three single-run or replicate-average figures: Kuhn number fraction, Kuhn weight fraction,
% and end-to-end distance. Usage snapshots run left to right in each.
% Each figure shares x/y limits across its own snapshots only.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(projectRoot,'SRC'));
% Supply one run folder, or multiple replicates with identical physics and distinct seeds.
% A single run has fitted curves but no replicate SD shading or bootstrap CI.
runFolders = {'E:\PhD\My Research\Polydisperse_network_growth\replicated_composition_20261102_201637_118_f4\replicated_composition_20261102_201637_118\run_001_rep_01_xA_0.01_seed_11231_20261102_201638_431',
               'E:\PhD\My Research\Polydisperse_network_growth\replicated_composition_20261102_201637_118_f4\replicated_composition_20261102_201637_118\run_005_rep_02_xA_0.01_seed_9_20261103_001430_152',
               'E:\PhD\My Research\Polydisperse_network_growth\replicated_composition_20261102_201637_118_f4\replicated_composition_20261102_201637_118\run_009_rep_03_xA_0.01_seed_7_20261103_025431_731'}; % supply saved run folders, or select a batch below
% Alternatively select a batch and xA; leave runFolders empty to use these:
batchFolder = ''; % contains batch_summary.mat from run_sweep_growing
selectedXA = .05;

averageOptions = struct();
averageOptions.usagePercent = [20 40 60 80 100];
averageOptions.nBinWidth = 1;       % integer Kuhn bin width; larger = coarser histogram
averageOptions.distanceBins = 30;
averageOptions.kuhnYLim = [0 .45];   % fixed number-fraction y range; [] = automatic
averageOptions.weightYLim = [];    % weight figure's own shared range; independent of number figure
averageOptions.showSpread = true;  % shaded +/- one replicate SD, not confidence intervals
averageOptions.requireAllReplicates = true; % avoid changing replicate cohorts at later usage
averageOptions.visible = 'on';
averageOptions.fitFlory = true;    % fit number fractions; required for number-predicted weight curves
averageOptions.fitFloryWeight = false; % true adds an independent weight fit; false skips that fit
averageOptions.showWeightFromNumberFit = true; % show weight theory using number-fit p, even when fitFloryWeight=false
averageOptions.floryModel = 'shifted'; % 'shifted' or 'original'
averageOptions.floryMinCount = 2;  % matches this run's minimum; shifted: n_min in (1-p)*p^(n-n_min)
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
runFolders = cellstr(runFolders);
if isempty(runFolders)
    error('NetworkGrowth:ReplicateFolders','Set at least one run folder, or batchFolder and selectedXA.');
end
outputRoot = fullfile(projectRoot,'replicate_averages','growing_chains');
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
fprintf('Distribution plots saved to: %s\n',averageOutputDir);


