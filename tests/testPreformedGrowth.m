function tests = testPreformedGrowth
tests = functiontests(localfunctions);
end
function setupOnce(tc)
tc.TestData.oldPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'SRC'));
tc.TestData.root = tempname; mkdir(tc.TestData.root);
end
function teardownOnce(tc)
path(tc.TestData.oldPath);
end
function cfg = smallConfig()
cfg = defaultNetworkBConfig();
cfg.Ntotal = 200; cfg.xA = .05; cfg.fA = 4; cfg.L = 12; cfg.b = 1.5;
cfg.randSeed = 17; cfg.chainModel = 'gaussian';
cfg.showLivePlot = false; cfg.makeVideo = false; cfg.showFinalHistogram = false;
cfg.saveFinalDistributions = false; cfg.saveData = false; cfg.saveDumps = false;
cfg.saveLammps = false;
end
function testFeedAllocationAndExplicitCounts(tc)
cfg = smallConfig(); cfg.chainLengths = [2 5 10]; cfg.chainNumberFractions = [.5 .3 .2];
feed = resolvePreformedFeed(cfg);
verifyEqual(tc,feed.initialCrosslinkers,10);
verifyEqual(tc,feed.initialChainCounts,[22 13 8]);
verifyEqual(tc,feed.initialKuhnCount,189); verifyEqual(tc,feed.unusedKuhnBudget,1);
verifyEqual(tc,feed.chainEndArmRatio,86/40);
cfg.chainCounts = [3 2 1]; feed = resolvePreformedFeed(cfg);
verifyEqual(tc,feed.initialChainCounts,[3 2 1]); verifyEqual(tc,feed.initialKuhnCount,26);
cfg.chainCounts = []; cfg.chainNumberFractions = [1 0 0]; feed = resolvePreformedFeed(cfg);
verifyEqual(tc,feed.initialChainCounts,[95 0 0]);
cfg.chainLengths = 20; cfg.chainNumberFractions = 1; feed = resolvePreformedFeed(cfg);
verifyEqual(tc,feed.initialChainCounts,9); verifyEqual(tc,feed.unusedKuhnBudget,10);
end
function testFeedValidation(tc)
cfg = smallConfig(); cfg.chainLengths = [1 2 3 4]; cfg.chainNumberFractions = [.25 .25 .25 .25];
verifyError(tc,@()resolvePreformedFeed(cfg),'NetworkGrowth:ChainSpecies');
cfg.chainLengths = [2 2]; cfg.chainNumberFractions = [.5 .5];
verifyError(tc,@()resolvePreformedFeed(cfg),'NetworkGrowth:ChainSpecies');
cfg.chainLengths = [2 5]; cfg.chainNumberFractions = [.2 .3];
verifyError(tc,@()resolvePreformedFeed(cfg),'NetworkGrowth:ChainFractions');
cfg.chainCounts = [100 100];
verifyError(tc,@()resolvePreformedFeed(cfg),'NetworkGrowth:ChainBudget');
cfg.chainCounts = [0 0];
verifyError(tc,@()resolvePreformedFeed(cfg),'NetworkGrowth:NoPreformedChains');
cfg = smallConfig(); cfg.nucleation = true;
verifyError(tc,@()runNetworkB_cleaned(cfg),'NetworkGrowth:PreformedNucleation');
end
function testEventPropensitiesAndReactants(tc)
cfg = smallConfig(); cfg.L = 10; cfg.b = 2; cfg.kAddition = 3;
cfg.kChem = 5; cfg.kEndEnd = 7; cfg.allowSelfLoops = true;
nodes = [1 1 0;9 1 0]; bonds = zeros(0,5);
ends = [1 1 1 0;1 1 1 0;1 1 1 6;9 1 2 0;9 1 2 9];
r = buildPreformedChainEvents(nodes,bonds,ends,4,[8 2 0],cfg);
addition = r(r(:,1)==1,:); fresh = r(r(:,1)==2,:); existing = r(r(:,1)==3,:);
verifyEqual(tc,size(addition,1),6); verifyEqual(tc,sum(addition(:,4)),3*3*10/100,'AbsTol',1e-14);
verifyTrue(tc,all(ends(addition(:,2),4)==0)); verifyFalse(tc,any(addition(:,3)==3));
verifyEqual(tc,fresh(:,2),[3;5]); verifyEqual(tc,fresh(:,4),repmat(4*5*4/100,2,1));
verifyTrue(tc,all(ends(existing(:,2),4)>0 & ends(existing(:,3),4)==0));
inter = existing(existing(:,2)==3 & existing(:,3)==4,:);
verifyEqual(tc,inter(4),7*5/(pi*6*4)*exp(-4/(6*4)),'RelTol',1e-14);
loops = existing(existing(:,2)==3 & existing(:,3)<=2,:);
verifyEqual(tc,loops(:,4),repmat(.01*5/(pi*6*4)*exp(-1),2,1),'RelTol',1e-14);
cfg.kSelfLoop = 0;
r = buildPreformedChainEvents(nodes,bonds,ends,0,[0 0 0],cfg);
verifyTrue(tc,all(r(:,1)==3)); verifyTrue(tc,all(ends(r(:,2),3)~=ends(r(:,3),3)));
end
function testNoParallelStrandsNoSecondLoopNoDanglingJoining(tc)
cfg = smallConfig(); cfg.allowSelfLoops = true;
nodes = [1 1 0;9 1 0]; ends = [1 1 1 6;1 1 1 0;9 1 2 5;9 1 2 0];
bonds = [1 2 20 0 0;1 1 20 0 0;2 2 20 0 0];
verifyEmpty(tc,buildPreformedChainEvents(nodes,bonds,ends,0,0,cfg));
verifyEmpty(tc,buildPreformedChainEvents(nodes,zeros(0,5),ends([1 3],:),0,0,cfg));
% N=1 is valid in Case B; Case A's joined-length threshold must not suppress it.
ends(1,4) = 1; cfg.allowSelfLoops = false;
r = buildPreformedChainEvents(nodes,zeros(0,5),ends([1 4],:),0,0,cfg);
verifyEqual(tc,r(:,1:3),[3 1 2]);
end
function testMonoBiTriTrajectoriesAndConservation(tc)
for species = 1:3
    cfg = smallConfig(); lengths = [2 5 10]; cfg.chainLengths = lengths(1:species);
    cfg.chainNumberFractions = ones(1,species)/species;
    cfg.allowSelfLoops = true; cfg.saveDumps = true; cfg.dumpEveryEvents = 1;
    cfg.outputDir = fullfile(tc.TestData.root,sprintf('species_%d',species));
    evalc('r = runNetworkB_cleaned(cfg);');
    [states,metadata] = loadNetworkDumps(cfg.outputDir);
    verifyEqual(tc,metadata.initialMonomers,0); verifyEqual(tc,metadata.schemaVersion,2);
    verifyEqual(tc,r.chainCountHistory(1,:),r.config.preformedFeed.initialChainCounts);
    verifyEqual(tc,size(r.chainCountHistory,1),size(r.history,1)+1);
    verifyTrue(tc,all(ismember(r.history(:,4),[1 2 3])));
    for k = 1:numel(states)
        s = states{k}; p = s.preformed;
        verifyEqual(tc,p.freeChainCounts+p.danglingChainCounts+p.completedChainCounts,p.initialChainCounts);
        verifyEqual(tc,p.freeKuhnSegments+p.incorporatedKuhnSegments,metadata.preformedFeed.initialKuhnCount);
        verifyEqual(tc,s.nodeCount+s.freeCrosslinkers,metadata.initialCrosslinkers);
        verifyEqual(tc,s.freeMonomers,0);
        verifyTrue(tc,all(ismember(s.bonds(:,3),cfg.chainLengths)));
        d = s.activeEnds(s.activeEnds(:,4)>0,:);
        verifyTrue(tc,all(ismember(d(:,4),cfg.chainLengths)));
        for node = 1:s.nodeCount
            verifyEqual(tc,s.nodes(node,3)+sum(s.activeEnds(:,3)==node),cfg.fA);
        end
        validateNetworkState(s.nodes,s.bonds,s.activeEnds,cfg.fA,false);
        if k>1 && strcmp(s.phase,'event')
            prev = states{k-1}; type = r.history(s.event,4);
            switch type
                case 1
                    verifyEqual(tc,sum(prev.preformed.freeChainCounts-p.freeChainCounts),1);
                    verifyEqual(tc,s.bondCount,prev.bondCount);
                    verifyEqual(tc,sum(p.danglingChainCounts-prev.preformed.danglingChainCounts),1);
                case {2,3}
                    verifyEqual(tc,p.freeChainCounts,prev.preformed.freeChainCounts);
                    verifyEqual(tc,s.bondCount,prev.bondCount+1);
                    verifyEqual(tc,sum(p.danglingChainCounts-prev.preformed.danglingChainCounts),-1);
            end
        end
    end
end
end
function testLoopConsumesOneChainAndTwoArms(tc)
cfg = smallConfig(); cfg.Ntotal = 21; cfg.xA = 1/21; cfg.fA = 2;
cfg.chainLengths = 10; cfg.kSelfLoop = 1e10; cfg.allowSelfLoops = true;
evalc('r = runNetworkB_cleaned(cfg);');
verifyEqual(tc,r.bonds,[1 1 10 0 0]); verifyEqual(tc,r.nodes(:,3),2);
verifyEqual(tc,r.preformed.completedChainCounts,1); verifyEqual(tc,r.preformed.freeChainCounts,1);
verifyEmpty(tc,r.activeEnds);
cfg.allowSelfLoops = false;
evalc('disabled = runNetworkB_cleaned(cfg);');
verifyEmpty(tc,disabled.bonds); verifyEqual(tc,disabled.preformed.danglingChainCounts,2);
verifyEqual(tc,disabled.stopReason,'no_allowed_events');
end
function testFillRangeDoesNotShrinkWithScale(tc)
cfg = smallConfig(); cfg.Ntotal = []; cfg.L = 12; cfg.b = 1.5;
cfg.chainLengths = 2; cfg.rEqMode = 'fill'; cfg.NtotalScale = 1;
evalc('one = runNetworkB_cleaned(cfg);');
cfg.NtotalScale = 2;
evalc('two = runNetworkB_cleaned(cfg);');
verifyEqual(tc,one.rEq,two.rEq);
verifyGreaterThan(tc,two.config.preformedFeed.initialKuhnCount,one.config.preformedFeed.initialKuhnCount);
end
function testFreshPlacementUsesFixedGaussianCoil(tc)
cfg = smallConfig(); cfg.Ntotal = 18; cfg.xA = 2/18; cfg.fA = 2;
cfg.chainLengths = 16; cfg.chainCounts = 1;
cfg.allowSelfLoops = false; cfg.kT = 0; cfg.kRep = 0; cfg.kRepAnh = 0;
% Two SSA draws per event, then the radial draw and angle at fresh capture.
rng(cfg.randSeed,'twister'); draws = rand(1,6);
expected = cfg.b*sqrt(16)*sqrt(-log(draws(5)))*[cos(2*pi*draws(6)),sin(2*pi*draws(6))];
evalc('r = runNetworkB_cleaned(cfg);');
verifyEqual(tc,r.history(:,4),[1;2]); verifyEqual(tc,r.bonds(:,3),16);
verifyEqual(tc,r.nodes(2,1:2)-r.nodes(1,1:2),expected,'AbsTol',1e-13);
end
