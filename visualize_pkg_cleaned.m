% Offline visualization driver. Edit inputFolder, then run this file.
% Accepts ONE simulation run folder or its dumps subfolder. No simulation runs.
addpath(fullfile(fileparts(mfilename('fullpath')),'SRC'));

inputFolder = 'E:\PhD\My Research\Polydisperse_network_growth\Code\Network_Growth\runs_pkg_cleaned\single_xA_0p0500_N_6714_20261026_122029_480'; % e.g. 'E:\...\runs_pkg_cleaned\single_xA_...'; '' opens a folder picker

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

% Individual checkpoint panels: network, Kuhn counts, end-to-end lengths.
vizOptions.saveCheckpointPlots = true;
vizOptions.checkpointUsagePercent = 0:5:100;
vizOptions.saveFinalPlot = true; % final state after relaxation, separate from 100%
% Histogram bins and axes are fixed across ALL saved states, including the final.
% [] chooses a shared range automatically; set limits to compare different runs.
vizOptions.kuhnXLim = [];        % contour length n, e.g. [0 100]
vizOptions.endToEndXLim = [];    % end-to-end distance, e.g. [0 30]
vizOptions.kuhnYLim = [];        % fixed Kuhn-count y range, e.g. [0 50]
vizOptions.endToEndYLim = [];    % separate fixed end-to-end y range, e.g. [0 20]
vizOptions.endToEndBins = 30;    % shared end-to-end bins

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
