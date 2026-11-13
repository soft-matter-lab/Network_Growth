% Final-state paths + rigid clusters. Accepts one run/dumps folder or a batch.
% No growth or relaxation. Shared analysis also supports preformed chains.
projectRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(projectRoot,'SRC'));
inputFolder=''; % a run, sweep/batch, or unpacked simulation handoff folder

connectivityOptions=struct();
connectivityOptions.expectedGrowthMode='growing';
connectivityOptions.outputRoot=''; % '' = new connectivity_final_* folder inside inputFolder
connectivityOptions.saveSnapshots=true; % paths_x.png, paths_y.png, rigid_clusters.png
connectivityOptions.snapshotResolution=600;
connectivityOptions.showNodeConnectivity=true; % gray -> black in path snapshots
connectivityOptions.continueOnError=true; % batch failures reported explicitly, never zeroed
% Same numerical settings as the paper analysis; normally leave unchanged.
connectivityOptions.rankTolerance=1e-9;
connectivityOptions.pairTolerance=1e-8;
connectivityOptions.imageRange=2;
connectivityOptions.checkRankStability=true;

if isempty(inputFolder)
    inputFolder=uigetdir(projectRoot,'Select a growing-chain run or batch folder');
    if isequal(inputFolder,0),disp('Connectivity analysis cancelled.');return;end
end
connectivityReport=postprocessNetworkConnectivity(inputFolder,connectivityOptions);
disp(connectivityReport.summary);
