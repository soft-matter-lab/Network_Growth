function tests = testSelfLoopKinetics
tests = functiontests(localfunctions);
end

function setupOnce(tc)
tc.TestData.oldPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'SRC'));
tc.TestData.root = tempname; mkdir(tc.TestData.root);
fprintf('Self-loop test artifacts: %s\n',tc.TestData.root);
end

function teardownOnce(tc)
path(tc.TestData.oldPath);
end

function cfg = config()
cfg = defaultNetworkAConfig(); cfg.L = 10; cfg.b = 2;
cfg.kChem = 3; cfg.kEndEnd = 2; cfg.allowSelfLoops = true;
cfg.kSelfLoop = .01; cfg.selfLoopDistanceFactor = 1;
cfg.minKuhnEndEnd = 3;
end

function testBothMechanismsHaveCombinedLengthPenalty(tc)
cfg = config(); nodes = [5 5 0]; bonds = zeros(0,5);
oneChain = [5 5 1 6;5 5 1 0];
twoChains = [5 5 1 2;5 5 1 4];
a = buildEndJoiningEvents(nodes,bonds,oneChain,cfg);
b = buildEndJoiningEvents(nodes,bonds,twoChains,cfg);
verifyEqual(tc,a(:,1:3),[3 1 2]); verifyEqual(tc,b(:,1:3),[3 1 2]);
expected = .01*3/(pi*6*4)*exp(-1);
verifyEqual(tc,a(4),expected,'RelTol',1e-14);
verifyEqual(tc,b(4),expected,'RelTol',1e-14);
% The requested RMS-distance prescription introduces a constant e^-1 factor.
cfg.selfLoopDistanceFactor = 0;
zeroDistance = buildEndJoiningEvents(nodes,bonds,oneChain,cfg);
verifyEqual(tc,a(4)/zeroDistance(4),exp(-1),'RelTol',1e-14);
cfg.selfLoopDistanceFactor = 2;
larger = buildEndJoiningEvents(nodes,bonds,twoChains,cfg);
verifyEqual(tc,larger(4)/zeroDistance(4),exp(-4),'RelTol',1e-14);
end

function testIndependentRateAndEnableSwitch(tc)
cfg = config(); nodes = [1 1 0;9 1 0]; bonds = zeros(0,5);
ends = [1 1 1 6;1 1 1 0;9 1 2 6];
allEvents = buildEndJoiningEvents(nodes,bonds,ends,cfg);
isLoop = ends(allEvents(:,2),3)==ends(allEvents(:,3),3);
cfg.kSelfLoop = .02;
doubled = buildEndJoiningEvents(nodes,bonds,ends,cfg);
verifyEqual(tc,doubled(isLoop,4),2*allEvents(isLoop,4),'RelTol',1e-14);
verifyEqual(tc,doubled(~isLoop,:),allEvents(~isLoop,:));
% End joining on different anchors still uses the periodic distance of 2.
expected = cfg.kEndEnd*cfg.kChem/(pi*12*cfg.b^2)*exp(-4/(12*cfg.b^2));
row = allEvents(:,2)==1 & allEvents(:,3)==3;
verifyEqual(tc,allEvents(row,4),expected,'RelTol',1e-14);
cfg.kEndEnd = 0;
loopsOnly = buildEndJoiningEvents(nodes,bonds,ends,cfg);
verifyEqual(tc,loopsOnly,doubled(isLoop,:));
cfg.kEndEnd = 2; cfg.kSelfLoop = 0;
verifyEqual(tc,buildEndJoiningEvents(nodes,bonds,ends,cfg),allEvents(~isLoop,:));
cfg.kSelfLoop = .01; cfg.allowSelfLoops = false;
verifyEqual(tc,buildEndJoiningEvents(nodes,bonds,ends,cfg),allEvents(~isLoop,:));
end

function testAvailableArmsAndEligibility(tc)
cfg = config(); nodes = [5 5 0]; bonds = zeros(0,5);
% Each spare bare arm is one destination; unordered pairs are counted once.
ends = [5 5 1 6;5 5 1 0;5 5 1 0];
events = buildEndJoiningEvents(nodes,bonds,ends,cfg);
verifyEqual(tc,events(:,1:3),[3 1 2;3 1 3]);
verifyEqual(tc,events(1,4),events(2,4));
verifyEmpty(tc,buildEndJoiningEvents(nodes,bonds,ends(1,:),cfg));
ends(1,4) = 5;
verifyEmpty(tc,buildEndJoiningEvents(nodes,bonds,ends,cfg));
cfg.minKuhnEndEnd = 0;
verifyEmpty(tc,buildEndJoiningEvents(nodes,bonds,ends(2:3,:),cfg));
% The existing no-duplicate-bond rule still permits only one loop per anchor.
ends(1,4) = 6;
verifyEmpty(tc,buildEndJoiningEvents(nodes,[1 1 6 0 0],ends,cfg));
end

function testBothLoopModesExecuteAndConserveMaterial(tc)
cfg = smallConfig(); cfg.kSelfLoop = 1e6; cfg.saveDumps = true;
cfg.dumpUsagePercent = []; cfg.dumpEveryEvents = 1;
observed = false(1,2);
for seed = 1:20
    cfg.randSeed = seed; cfg.outputDir = fullfile(tc.TestData.root,sprintf('mode_%d',seed));
    evalc('r = runNetworkA_cleaned(cfg);');
    verifyEqual(tc,size(r.bonds,1),1);
    verifyEqual(tc,r.bonds(1,[1 2 4 5]),[1 1 0 0]);
    verifyEqual(tc,r.nodes(:,3),2); verifyEmpty(tc,r.activeEnds);
    verifyEqual(tc,sum(r.bonds(:,3))+r.nB,12); verifyEqual(tc,r.nA,0);
    states = loadNetworkDumps(cfg.outputDir);
    beforeLoop = states{end-2}; afterLoop = states{end-1};
    n = beforeLoop.activeEnds(:,4);
    verifyEqual(tc,numel(n),2); verifyEqual(tc,sum(n),r.bonds(1,3));
    verifyEqual(tc,beforeLoop.freeMonomers,afterLoop.freeMonomers);
    verifyEqual(tc,beforeLoop.freeCrosslinkers,afterLoop.freeCrosslinkers);
    verifyEqual(tc,afterLoop.nodeCount,beforeLoop.nodeCount);
    verifyEqual(tc,r.history(end,4),3);
    if any(n==0), observed(1) = true; else, observed(2) = true; end
    if all(observed), break; end
end
verifyTrue(tc,all(observed),'The seeded runs must exercise both loop mechanisms.');
end

function testZeroRateMatchesLoopsDisabled(tc)
cfg = smallConfig(); cfg.randSeed = 17; cfg.kSelfLoop = 0;
evalc('zeroRate = runNetworkA_cleaned(cfg);'); after = rng;
cfg.allowSelfLoops = false; cfg.kSelfLoop = .01;
evalc('disabled = runNetworkA_cleaned(cfg);');
verifyEqual(tc,rng,after);
verifyEmpty(tc,zeroRate.bonds);
fields = {'nodes','bonds','activeEnds','history','t','nA','nB'};
for k = 1:numel(fields), verifyEqual(tc,zeroRate.(fields{k}),disabled.(fields{k})); end
end

function testSweepUsesLoopSettings(tc)
cfg = smallConfig(); cfg.kSelfLoop = .01; cfg.selfLoopDistanceFactor = 1.5;
sweep = struct('mode','replicates','seeds',[17 28],'saveResults',true);
evalc('[summary,batchDir] = runNetworkSweep(cfg,sweep,fullfile(tc.TestData.root,''sweep''));');
verifyEqual(tc,summary.status,repmat("completed",2,1));
saved = load(fullfile(batchDir,'batch_config.mat'),'baseCfg');
verifyEqual(tc,saved.baseCfg.kSelfLoop,.01);
for k = 1:height(summary)
    saved = load(fullfile(summary.runDir(k),'results_struct.mat'),'results');
    r = saved.results;
    verifyTrue(tc,r.config.allowSelfLoops);
    verifyEqual(tc,r.config.kSelfLoop,.01);
    verifyEqual(tc,r.config.selfLoopDistanceFactor,1.5);
    verifyEqual(tc,size(r.bonds,1),1);
    verifyEqual(tc,r.bonds(:,1:2),[1 1]);
end
end

function cfg = smallConfig()
cfg = defaultNetworkAConfig();
cfg.Ntotal = 13; cfg.xA = 1/13; cfg.fA = 2; cfg.L = 10; cfg.b = 2;
cfg.allowSelfLoops = true; cfg.kEndEnd = 0; cfg.minKuhnEndEnd = 1;
cfg.showLivePlot = false; cfg.makeVideo = false; cfg.showFinalHistogram = false;
cfg.saveFinalDistributions = false; cfg.saveData = false; cfg.saveDumps = false;
cfg.saveLammps = false;
end
