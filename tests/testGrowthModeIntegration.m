function tests = testGrowthModeIntegration
tests = functiontests(localfunctions);
end
function setupOnce(tc)
tc.TestData.oldPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'SRC'));
tc.TestData.root = tempname; mkdir(tc.TestData.root);
fprintf('Case B integration artifacts: %s\n',tc.TestData.root);
cfg = defaultNetworkBConfig(); cfg.Ntotal = 200; cfg.xA = .05; cfg.fA = 4;
cfg.L = 12; cfg.b = 1.5; cfg.randSeed = 17; cfg.chainModel = 'gaussian';
cfg.chainLengths = [2 5 10]; cfg.chainNumberFractions = [.5 .3 .2];
cfg.allowSelfLoops = true;
cfg.showLivePlot = false; cfg.makeVideo = false; cfg.showFinalHistogram = false;
cfg.saveFinalDistributions = false; cfg.saveData = true; cfg.saveDumps = true;
cfg.dumpEveryEvents = 1; cfg.dumpUsagePercent = 0:25:100;
cfg.saveLammps = true; cfg.lammpsNRho = 101;
cfg.outputDir = fullfile(tc.TestData.root,'single');
evalc('r = runNetworkB_cleaned(cfg);');
tc.TestData.cfg = cfg; tc.TestData.run = r;
end
function teardownOnce(tc)
close all; path(tc.TestData.oldPath);
end
function testLammpsAndSavedStateBuilder(tc)
r = tc.TestData.run; cfg = tc.TestData.cfg;
verifyTrue(tc,startsWith(r.lammps.status,'written'));
verifyEqual(tc,r.lammps.totalKuhn,r.preformed.incorporatedKuhnSegments);
verifyEqual(tc,r.lammps.pairRepulsionCutoff,cfg.b*sqrt(r.preformed.incorporatedKuhnSegments/size(r.nodes,1)),'RelTol',1e-14);
verifyEqual(tc,r.lammps.bondCount,sum(r.bonds(:,1)~=r.bonds(:,2)));
destination = fullfile(tc.TestData.root,'rebuilt');
evalc('info = buildLammpsIngredients(cfg.outputDir,struct(''outputDir'',destination));');
verifyTrue(tc,isfile(fullfile(destination,'network.data')));
% The raw MAT source is sufficient; rebuilding requires no growth execution.
saved = load(r.dataFile);
verifyEqual(tc,saved.preformed,r.preformed); verifyEqual(tc,saved.chainCountHistory,r.chainCountHistory);
end
function testSweepParameterSetsAndSavedMetadata(tc)
cfg = tc.TestData.cfg; cfg.saveData = false; cfg.saveLammps = false;
sweep = struct('mode','replicates','seeds',[17 28],'saveResults',true);
sweep.parameterSets = {struct('chainLengths',2,'chainNumberFractions',1,'fA',3), ...
    struct('chainLengths',[2 5],'chainNumberFractions',[.5 .5],'fA',4)};
plan = makeNetworkSweepPlan(cfg,sweep);
verifyEqual(tc,height(plan),4); verifyEqual(tc,plan.parameterSetIndex,[1;1;2;2]);
verifyEqual(tc,plan.fA,[3;3;4;4]);
evalc('[summary,batchDir] = runNetworkSweep(cfg,sweep,fullfile(tc.TestData.root,''sweep''));');
verifyTrue(tc,all(summary.status=="completed")); verifyTrue(tc,isfile(fullfile(batchDir,'run_plan.csv')));
for k = 1:height(summary)
    data = load(fullfile(summary.runDir(k),'results_struct.mat'),'results'); r = data.results;
    verifyEqual(tc,r.growthMode,'preformed'); verifyEqual(tc,r.config.fA,plan.fA(k));
    verifyEqual(tc,summary.freeMonomers(k),0);
    verifyEqual(tc,summary.freeChains(k),sum(r.preformed.freeChainCounts));
    verifyEqual(tc,sum(r.config.preformedFeed.initialChainCounts),plan.initialChains(k));
end
end
function testPostprocessingAndVideo(tc)
options = struct('makeVideo',true,'videoFrameEvery',1000,'videoHoldLastSeconds',0, ...
    'videoFPS',2,'heatGrid',25,'saveCheckpointPlots',false, ...
    'saveFinalPlot',true,'saveComparison',false,'showComparisonFigure',false);
evalc('report = postprocessNetworkDumps(tc.TestData.cfg.outputDir,options);');
verifyTrue(tc,isfile(report.finalFile)); verifyTrue(tc,isfile(report.videoFile));
verifyTrue(tc,isfile(report.preformedProgressFile));
progress = readtable(fullfile(report.outputDir,'preformed_progress.csv'));
verifyEqual(tc,progress.freeChains_N2(end),tc.TestData.run.preformed.freeChainCounts(1));
states = loadNetworkDumps(tc.TestData.cfg.outputDir);
lambda = networkDumpPrestretch(states{end});
verifyEqual(tc,numel(lambda),size(tc.TestData.run.bonds,1));
% Exercise the route-aware plot handles independently of the movie encoder.
cfg = tc.TestData.run.config; cfg.heatGrid = 25;
viz = initVisualization(cfg,0,9,cfg.preformedFeed.initialChainCounts,'off');
cleanup = onCleanup(@()close(viz.fig));
verifyEqual(tc,numel(viz.hB),3);
verifyTrue(tc,contains(viz.hB(2).DisplayName,'N = 5'));
end
function testAverageDefaultsAndMixtureGuard(tc)
options = struct('usagePercent',25,'visible','off');
evalc('[fig,result] = averageNetworkDumps({tc.TestData.cfg.outputDir},options);');
cleanup = onCleanup(@()close(fig));
verifyEmpty(tc,result.floryFitTable); verifyEmpty(tc,result.floryWeightFitTable);
cfg = tc.TestData.cfg; cfg.randSeed = 28; cfg.chainNumberFractions = [.3 .5 .2];
cfg.outputDir = fullfile(tc.TestData.root,'other_mixture'); cfg.saveLammps = false;
evalc('runNetworkB_cleaned(cfg);');
verifyError(tc,@()averageNetworkDumps({tc.TestData.cfg.outputDir,cfg.outputDir},options), ...
    'NetworkGrowth:ReplicateSettings');
end
