function tests = testPerformanceEquivalence
% Exact comparison against the frozen pre-optimization implementation.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
tc.TestData.oldPath = path;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
tc.TestData.root = tempname; mkdir(tc.TestData.root);
tc.TestData.baseline = fullfile(tc.TestData.root,'baseline'); mkdir(tc.TestData.baseline);
source = fullfile(root,'tests','reference','performance_baseline');
files = dir(fullfile(source,'*.m.txt'));
for k = 1:numel(files)
    copyfile(fullfile(source,files(k).name),fullfile(tc.TestData.baseline,files(k).name(1:end-4)));
end
end

function teardownOnce(tc)
path(tc.TestData.oldPath); clearImplementations();
end

function clearImplementations()
clear buildEndJoiningEvents buildGrowingChainEvents buildPreformedChainEvents
clear chainSpring repulsionPairs equilibrateNetwork
end

function useBaseline(tc)
addpath(tc.TestData.baseline,'-begin'); clearImplementations();
assert(strcmp(which('buildEndJoiningEvents'),fullfile(tc.TestData.baseline,'buildEndJoiningEvents.m')));
end

function useCurrent(tc)
if contains([path pathsep],[tc.TestData.baseline pathsep]), rmpath(tc.TestData.baseline); end
clearImplementations();
end

function testEventTablesAndRngExactlyMatch(tc)
% Invalid/unused reservoir entries are deliberately avoided; zero, tiny and
% underflowed physical rates, empty reactants and occupied anchors are covered.
restoreRng = rng; cleanupRng = onCleanup(@()rng(restoreRng)); %#ok<NASGU>
rng(104729,'twister'); inputs = cell(1,240);
for k = 1:numel(inputs)
    c = defaultNetworkBConfig();
    c.L = 7+mod(k,9); c.b = .3+mod(k,7)/3;
    c.allowSelfLoops = mod(k,2)==0;
    rates = [0 1e-16 .01 1 7];
    c.kSelfLoop = rates(1+mod(k,5)); c.kEndEnd = rates(1+mod(floor(k/5),5));
    c.kChem = rates(1+mod(floor(k/25),5)); c.kAddition = rates(1+mod(k+2,5));
    c.selfLoopDistanceFactor = mod(k,4); c.minKuhnEndEnd = mod(k,5);
    c.minKuhnCrosslink = mod(k,3); c.nucleation = mod(k,3)==0;
    c.kNuc = rates(1+mod(k,5)); c.fA = 2+mod(k,5);
    N = 1+mod(k,19); M = mod(k,47);
    nodes = [c.L*(4*rand(N,2)-1.5),zeros(N,1)];
    parent = randi(N,M,1); lengths = randi([0 30],M,1);
    lengths(1:3:end) = 0;
    ends = [nodes(parent,1:2),parent,lengths];
    nb = mod(k,12);
    bonds = [randi(N,nb,2),ones(nb,1),zeros(nb,2)];
    counts = [0 3 11]; counts = counts(1:1+mod(k,3));
    if mod(k,2), counts = counts'; end
    inputs{k} = {nodes,bonds,ends,mod(k,7),mod(k,11),counts,c};
end
% Periodic half-box ties, a zero-distance inter-anchor pair, both self-loop
% mechanisms, N=1 Case-B strands, and rates around Case A's 1e-14 cutoff.
for factor = [.5 1 2]
    c = defaultNetworkBConfig(); c.L = 10; c.b = 1;
    c.allowSelfLoops = true; c.kSelfLoop = 1; c.selfLoopDistanceFactor = 0;
    c.kChem = factor*1e-14*(2*pi); c.minKuhnEndEnd = 0;
    nodes = [0 0 0;5 -5 0;10 0 0];
    ends = [0 0 1 1;0 0 1 0;5 -5 2 1;10 0 3 1;10 0 3 0];
    inputs{end+1} = {nodes,zeros(0,5),ends,0,0,[0 0 0],c}; %#ok<AGROW>
end
useBaseline(tc); cleanup = onCleanup(@()useCurrent(tc)); %#ok<NASGU>
expected = evaluateTables(inputs); useCurrent(tc);
before = rng; actual = evaluateTables(inputs);
for k = 1:numel(actual)
    verifyEqual(tc,actual{k},expected{k},sprintf('Event-table comparison %d',k));
end
verifyEqual(tc,rng,before);
end

function outputs = evaluateTables(inputs)
outputs = cell(numel(inputs),3);
for k = 1:numel(inputs)
    a = inputs{k}; nodes=a{1}; bonds=a{2}; ends=a{3}; c=a{7};
    outputs{k,1} = buildEndJoiningEvents(nodes,bonds,ends,c);
    outputs{k,2} = buildGrowingChainEvents(nodes,bonds,ends,a{4},a{5},c);
    outputs{k,3} = buildPreformedChainEvents(nodes,bonds,ends,a{4},a{6},c);
end
end

function testForcesEnergiesPairsAndRelaxationExactlyMatch(tc)
rngBefore = rng; cleanupRng = onCleanup(@()rng(rngBefore)); %#ok<NASGU>
rng(577,'twister');
c = quietConfig(); c.L=12; c.b=1.2; c.kRepAnh=.03;
nodes = [12*(2*rand(35,2)-.5),zeros(35,1)]; nodes(2,1:2)=nodes(1,1:2);
bonds = [(1:34)',(2:35)',randi([2 20],34,1),zeros(34,2);1 1 5 0 0];
ends = zeros(0,4);
useBaseline(tc); cleanup = onCleanup(@()useCurrent(tc)); %#ok<NASGU>
expected = mechanicsSamples(nodes,bonds,ends,c);
useCurrent(tc); before=rng;
actual = mechanicsSamples(nodes,bonds,ends,c);
verifyEqual(tc,actual,expected); verifyEqual(tc,rng,before);
end

function out = mechanicsSamples(nodes,bonds,ends,c)
out = {};
for model = {'gaussian','cohen'}
    c.chainModel = model{1};
    r = [0;1e-14;.95-eps(.95);.95;.95+eps(.95);1;10]*c.b;
    n = ones(size(r)); [coefficient,energy]=chainSpring(r,n,c.b,c.kT,c.chainModel);
    forceOnly=chainSpring(r,n,c.b,c.kT,c.chainModel);
    out{end+1}={coefficient,energy,forceOnly}; %#ok<AGROW>
    for binning = [false true]
        for repelBonded = [false true]
            for radius = [0 2 7]
                args={c.b,c.kT,c.kRep,radius,c.L,repelBonded,binning,0,c.kRepAnh,c.chainModel};
                force=networkForces(nodes,bonds,args{:});
                [total,eb,er]=networkEnergy(nodes,bonds,args{:});
                [ii,jj,rv,d]=repulsionPairs(nodes,bonds,c.L,radius,repelBonded,binning,0);
                out{end+1}={force,total,eb,er,ii,jj,rv,d}; %#ok<AGROW>
            end
        end
    end
    for integrator = {'damped2nd','overdamped'}
        for damping = [0 .5]
            c.eqIntegrator=integrator{1}; c.gammaEq=damping;
            c.maxEq=25; c.tolV=1e-5;
            [relaxed,speed,it]=equilibrateNetwork(nodes,bonds,ends,c,2,[]);
            out{end+1}={relaxed,speed,it}; %#ok<AGROW>
        end
    end
end
% Exercise cached pair indices across changing N/coordinates/configuration.
for N = [1 2 35 7 35 35]
    xyz=nodes(1:N,:);
    xyz(:,1)=xyz(:,1)+(.05*N)*(1:N)';
    [i,j,v,d]=repulsionPairs(xyz,zeros(0,5),12,3,true,false,120);
    out{end+1}={i,j,v,d}; %#ok<AGROW>
end
% Above the memory cap, then return to a previously cached size.
for N = [1420 35]
    xyz=[mod((1:N)'*.731,12),mod((1:N)'*.419,12),zeros(N,1)];
    [i,j,v,d]=repulsionPairs(xyz,[1 2 3 0 0],12,.6,false,false,120);
    out{end+1}={i,j,v,d}; %#ok<AGROW>
end
end

function testSeededTrajectoriesAndEverySavedStateExactlyMatch(tc)
cleanup = onCleanup(@()useCurrent(tc)); %#ok<NASGU>
savedRng=rng; cleanupRng=onCleanup(@()rng(savedRng)); %#ok<NASGU>
for variant = 1:8
  for replicate = 1:2
    fprintf('Comparing full trajectory variant %d/8, seed %d/2...\n',variant,replicate);
    c=quietConfig(); c.Ntotal=160; c.xA=.075; c.fA=4;
    c.L=12; c.b=1.4; c.randSeed=100+variant+1000*(replicate-1); c.rEqMode='fill';
    c.allowSelfLoops=mod(variant,2)==0;
    if mod(variant,2)==0, c.chainModel='cohen'; else, c.chainModel='gaussian'; end
    if variant==3
        c.nucleation=true; c.kNuc=.1; c.eqIntegrator='overdamped';
        c.repelBondedPairs=false; c.useRepulsionBinning=true; c.repulsionFullPairThreshold=0;
        % Also exercise the existing iteration-cap behavior in this stiff
        % stress case, without making the regression suite take minutes.
        c.maxEq=80; c.maxEqFinal=200;
    elseif variant==4
        c.kSelfLoop=0; c.fA=3; c.rEqMode='feed';
    end
    if variant>=5
        c.growthMode='preformed'; ns=min(variant-4,3);
        c.chainLengths=[2 5 9]; c.chainLengths=c.chainLengths(1:ns);
        c.chainNumberFractions=ones(1,ns)/ns; c.chainCounts=[]; c.kAddition=1;
    end
    if variant==8
        c.Ntotal=21; c.xA=1/21; c.fA=2;
        c.chainLengths=10; c.chainNumberFractions=1; c.kSelfLoop=1e10;
    end
    c.saveDumps=ismember(variant,[2 6]); c.dumpEveryEvents=1;
    c.outputDir=fullfile(tc.TestData.root,sprintf('baseline_%d_%d',variant,replicate));
    useBaseline(tc); evalc('expected=runNetwork_cleaned(c);'); expectedRng=rng;
    if c.saveDumps, expectedDumps=loadNetworkDumps(c.outputDir); end
    c.outputDir=fullfile(tc.TestData.root,sprintf('current_%d_%d',variant,replicate));
    useCurrent(tc); evalc('actual=runNetwork_cleaned(c);');
    verifyEqual(tc,rng,expectedRng,sprintf('RNG differs for variant %d',variant));
    verifyEqual(tc,physicalResult(actual),physicalResult(expected),sprintf('Variant %d',variant));
    if c.saveDumps
        verifyEqual(tc,loadNetworkDumps(c.outputDir),expectedDumps);
    end
  end
end
end

function r = physicalResult(r)
% Only wall-clock timing and the two deliberately distinct output locations
% are excluded. Histories, energies, relaxation iterations and config remain.
r=rmfield(r,{'runtimeSeconds','outputDir','dataFile','dumpDir'});
r.config.outputDir=''; r.requestedConfig.outputDir='';
end

function c = quietConfig()
c=defaultNetworkAConfig();
c.showLivePlot=false; c.makeVideo=false; c.showFinalHistogram=false;
c.saveFinalDistributions=false; c.saveData=false; c.saveDumps=false; c.saveLammps=false;
end
