% Offline visualization driver. Edit inputFolder, then run this file.
% Accepts ONE simulation run folder or its dumps subfolder. No simulation runs.
addpath(fullfile(fileparts(mfilename('fullpath')),'SRC'));

inputFolder = 'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\batch_runs_pkg_cleaned_sweeps\replicates_20261026_164546_823\run_001_rep_01_xA_0.05_seed_11231_20261026_164548_463'; % e.g. 'E:\...\runs_pkg_cleaned\single_xA_...'; '' opens a folder picker

vizOptions = struct();
vizOptions.outputRoot = ''; % '' = a new visualization folder under the selected run

% Network movie: clouds + bonds + junctions, alongside reservoir histories.
vizOptions.makeVideo = true;
vizOptions.videoFrameEvery = 1; % use every Nth saved snapshot; always include the last
vizOptions.videoFPS = 6;
vizOptions.videoHoldLastSeconds = 2;
vizOptions.videoPreview = false; % show the figure while rendering
vizOptions.heatGrid = 150;
% One saved state per video frame, with equal frame duration. No interpolation.

% Checkpoint panels: network, chain number fraction, weight fraction, 2D displacement.
vizOptions.saveCheckpointPlots = true;
vizOptions.checkpointUsagePercent = 0:5:100;
vizOptions.saveFinalPlot = true; % final state after relaxation, separate from 100%
% Histogram bins and axes are fixed across ALL saved states, including the final.
% [] chooses a shared range automatically; set limits to compare different runs.
vizOptions.kuhnXLim = [];        % contour length n, e.g. [0 100]
vizOptions.weightXLim = [];     % [] = same n range as kuhnXLim, or e.g. [0 100]
vizOptions.kuhnYLim = [];        % fixed probability-per-bin y range, e.g. [0 .5]
vizOptions.weightYLim = [];     % independent fixed weight-fraction y range, e.g. [0 .2]
vizOptions.displacementBins = 40; % bins per axis in the 2D endpoint histogram
vizOptions.displacementLimit = []; % [] = run-wide auto; R means x,y span [-R,R]
vizOptions.displacementCLim = [];  % [] = run-wide probability color scale; e.g. [0 .05]

% Overlaid distributions: Kuhn counts, contour lengths n*b, end-to-end lengths.
vizOptions.saveComparison = true;
vizOptions.comparisonUsagePercent = [0 25 50 75 100]; % or 0:5:100 for every milestone
vizOptions.population = 'completed'; % 'completed', 'dangling', or 'all' for contour panels
vizOptions.showComparisonFigure = true;

if isempty(inputFolder)
    inputFolder = uigetdir(pwd,'Select one run folder or its dumps folder');
    if isequal(inputFolder,0)
        disp('Visualization cancelled.');
        return;
    end
end

visualizationReport = postprocessNetworkDumps(inputFolder,vizOptions);
fprintf('\nVisualization output: %s\n',visualizationReport.outputDir);
disp(visualizationReport);
