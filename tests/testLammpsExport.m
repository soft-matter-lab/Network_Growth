function tests = testLammpsExport
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath = path;
addpath(fullfile(root,'SRC'));
testCase.TestData.outputRoot = tempname;
mkdir(testCase.TestData.outputRoot);
fprintf('LAMMPS test artifacts: %s\n',testCase.TestData.outputRoot);
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function [cfg,nodes,bonds,ends] = fixture(testCase,label)
cfg = defaultNetworkAConfig();
cfg.outputDir = fullfile(testCase.TestData.outputRoot,label);
cfg.L = 10; cfg.b = 2;
cfg.lammpsNRho = 21; cfg.lammpsRhoRange = [0 20];
% First bond crosses a periodic boundary; first atom is unwrapped twice.
nodes = [11 -9 1;9 1 2;5 4 1];
bonds = [1 2 4 0 -1;2 3 8 0 0];
ends = [1 1 1 6;5 4 3 0];
end

function testFormatsAndRadius(testCase)
[cfg,nodes,bonds,ends] = fixture(testCase,'formats');
before = rng;
info = writeLammpsNetwork(nodes,bonds,ends,cfg);
verifyEqual(testCase,rng,before);
verifyEqual(testCase,info.status,'written');
verifyEqual(testCase,info.totalKuhn,18);
verifyEqual(testCase,info.atomCount,3);
verifyEqual(testCase,info.areaPerCrosslinker,6*pi,'AbsTol',1e-12);
verifyEqual(testCase,info.rRepulsive,sqrt(6),'AbsTol',1e-12);
verifyEqual(testCase,[info.RLower info.RUpper info.pairRepulsionCutoff],[.8 4 2]*sqrt(6),'AbsTol',1e-12);
verifyEmpty(testCase,info.nonMinimumImageBondRows);
files = dir(cfg.outputDir); files = files(~[files.isdir]);
verifyEqual(testCase,sort({files.name}),sort({'network.data','bond.table','network.localdensity.table'}));

data = fileread(info.dataFile);
verifyTrue(testCase,contains(data,sprintf('3 atoms\n2 bonds')));
verifyTrue(testCase,contains(data,'Atoms # bpm/sphere'));
atoms = sectionRows(data,'Atoms # bpm/sphere',3,8);
verifyEqual(testCase,atoms(:,1:5),[(1:3)',ones(3,4)]);
verifyEqual(testCase,atoms(:,6:8),[1 1 0;9 1 0;5 4 0]);
links = sectionRows(data,'Bonds',2,4);
verifyEqual(testCase,links,[1 1 1 2;2 1 2 3]);
% Reconstruct the periodic bond vectors from the serialized coordinates.
delta = atoms(links(:,4),6:7)-atoms(links(:,3),6:7);
delta = delta-cfg.L*round(delta/cfg.L);
verifyEqual(testCase,delta,[-2 0;-4 3]);

table = fileread(info.bondTableFile);
verifyTrue(testCase,contains(table,sprintf('KEY\nN 2\n\n')));
rows = sectionRows(table,'N 2',2,5);
verifyEqual(testCase,rows,[1 1 2 4 2;2 2 3 8 2]);
% Independently parse exactly the header-declared count, as LAMMPS does.
lines = splitlines(string(fileread(info.localDensityFile)));
verifyEqual(testCase,sscanf(lines(3),'%f')',[1 21]);
verifyEqual(testCase,sscanf(lines(5),'%f')',[.8 4]*sqrt(6),'AbsTol',1e-12);
verifyEqual(testCase,lines(6:7),["1";"1"]);
verifyEqual(testCase,sscanf(lines(8),'%f')',[0 20 1]);
values = str2double(lines(9:29));
verifyEqual(testCase,values,.414*((0:20)'-12.8).^2,'AbsTol',1e-12);
verifyTrue(testCase,all(strlength(strtrim(lines(30:end)))==0));
end

function testDanglingAndSelfLoopAccounting(testCase)
[cfg,nodes,bonds,ends] = fixture(testCase,'loops');
bonds = [bonds;2 2 3 0 0];
old = warning('off','NetworkGrowth:LammpsSelfLoops');
cleanup = onCleanup(@()warning(old)); %#ok<NASGU>
info = writeLammpsNetwork(nodes,bonds,ends,cfg);
verifyEqual(testCase,info.bondCount,2);
verifyEqual(testCase,info.omittedSelfLoopRows,3);
verifyEqual(testCase,info.exportedBondRows,[1;2]);
verifyEqual(testCase,info.completedKuhn,15);
verifyEqual(testCase,info.danglingKuhn,6);
verifyEqual(testCase,info.totalKuhn,21);
verifyEqual(testCase,info.rRepulsive,sqrt(7),'AbsTol',1e-12);
verifyEqual(testCase,sectionRows(fileread(info.dataFile),'Bonds',2,4),[1 1 1 2;2 1 2 3]);
cfg.outputDir = fullfile(testCase.TestData.outputRoot,'completed_only');
cfg.lammpsIncludeDanglingKuhn = false;
info = writeLammpsNetwork(nodes,bonds,ends,cfg);
verifyEqual(testCase,info.totalKuhn,15);
verifyEqual(testCase,info.rRepulsive,sqrt(5),'AbsTol',1e-12);
end

function testUnusedBackgroundExcluded(testCase)
[cfg,nodes,bonds,ends] = fixture(testCase,'background');
% 10 feed crosslinkers and 90 feed Kuhn segments, but only 3 formed
% crosslinkers and 18 chain-bound segments (12 completed + 6 dangling).
cfg.Ntotal = 100; cfg.xA = .1;
info = writeLammpsNetwork(nodes,bonds,ends,cfg);
verifyEqual(testCase,info.atomCount,3);
verifyEqual(testCase,info.totalKuhn,18);
verifyEqual(testCase,info.rRepulsive,sqrt(6),'AbsTol',1e-12);
% Changing only the unused feed reservoirs must leave the radius unchanged.
cfg.Ntotal = 200; cfg.xA = .2;
cfg.outputDir = fullfile(testCase.TestData.outputRoot,'more_background');
more = writeLammpsNetwork(nodes,bonds,ends,cfg);
verifyEqual(testCase,more.atomCount,info.atomCount);
verifyEqual(testCase,more.totalKuhn,info.totalKuhn);
verifyEqual(testCase,more.rRepulsive,info.rRepulsive);
end

function testPeriodicMismatchAndValidation(testCase)
[cfg,nodes,bonds,ends] = fixture(testCase,'mismatch');
bonds(1,4) = 1;
cfg.lammpsPeriodicPolicy = 'error';
verifyError(testCase,@()writeLammpsNetwork(nodes,bonds,ends,cfg),'NetworkGrowth:LammpsPeriodicImages');
verifyFalse(testCase,isfolder(cfg.outputDir));
cfg.lammpsPeriodicPolicy = 'warn';
old = warning('off','NetworkGrowth:LammpsPeriodicImages');
cleanup = onCleanup(@()warning(old)); %#ok<NASGU>
info = writeLammpsNetwork(nodes,bonds,ends,cfg);
verifyEqual(testCase,info.status,'written_with_image_changes');
verifyEqual(testCase,info.nonMinimumImageBondRows,1);
before = fileread(info.dataFile);
verifyError(testCase,@()writeLammpsNetwork(nodes,bonds,ends,cfg),'NetworkGrowth:ExistingLammpsFiles');
verifyEqual(testCase,fileread(info.dataFile),before);

[cfg,nodes,bonds,ends] = fixture(testCase,'invalid');
cfg.lammpsRhoRange = [10 0];
verifyError(testCase,@()writeLammpsNetwork(nodes,bonds,ends,cfg),'NetworkGrowth:LammpsDensityRange');
cfg.lammpsRhoRange = [0 20]; cfg.outputDir = '';
verifyError(testCase,@()writeLammpsNetwork(nodes,bonds,ends,cfg),'NetworkGrowth:LammpsOutputDir');
end

function testEmptyNetworksSkipExplicitly(testCase)
[cfg,nodes,~,ends] = fixture(testCase,'empty');
old = warning('off','NetworkGrowth:LammpsEmptyNetwork');
cleanup = onCleanup(@()warning(old)); %#ok<NASGU>
info = writeLammpsNetwork(nodes,zeros(0,5),ends,cfg);
verifyEqual(testCase,info.status,'skipped_no_interjunction_bonds');
verifyEqual(testCase,info.totalKuhn,6);
verifyEqual(testCase,info.rRepulsive,sqrt(2),'AbsTol',1e-12);
verifyEmpty(testCase,info.dataFile);
verifyFalse(testCase,isfolder(cfg.outputDir));
info = writeLammpsNetwork(nodes,[1 1 4 0 0],ends,cfg);
verifyEqual(testCase,info.status,'skipped_no_interjunction_bonds');
verifyEqual(testCase,info.omittedSelfLoopRows,1);
end

function testSingleRunPreservesTrajectory(testCase)
cfg = smallConfig(testCase,'single_on');
cfg.saveData = true; cfg.saveLammps = true;
evalc('withExport = runNetworkA_cleaned(cfg);');
after = rng;
cfg.outputDir = fullfile(testCase.TestData.outputRoot,'single_off');
cfg.saveLammps = false;
evalc('withoutExport = runNetworkA_cleaned(cfg);');
verifyEqual(testCase,rng,after);
fields = {'nodes','bonds','activeEnds','history','t','nA','nB','rEq','energy'};
for k = 1:numel(fields)
    verifyEqual(testCase,withExport.(fields{k}),withoutExport.(fields{k}));
end
verifyTrue(testCase,isfile(withExport.lammps.dataFile));
verifyFalse(testCase,isfile(fullfile(cfg.outputDir,'network.data')));
verifyEqual(testCase,withoutExport.lammps.status,'disabled');
verifyEqual(testCase,withExport.lammps.totalKuhn, ...
    sum(withExport.bonds(:,3))+sum(withExport.activeEnds(:,4)));
initialA = round(withExport.config.xA*withExport.config.Ntotal);
initialB = withExport.config.Ntotal-initialA;
verifyEqual(testCase,withExport.lammps.atomCount,initialA-withExport.nA);
verifyEqual(testCase,withExport.lammps.totalKuhn,initialB-withExport.nB);
saved = load(withExport.dataFile,'lammps');
verifyEqual(testCase,saved.lammps,withExport.lammps);
end

function testSweepWritesEachSample(testCase)
cfg = smallConfig(testCase,'unused'); cfg.saveLammps = true;
sweep = struct('mode','replicated_composition','xA',[.2 .3],'seeds',[17 28]);
out = fullfile(testCase.TestData.outputRoot,'batch');
evalc('[summary,batchDir] = runNetworkSweep(cfg,sweep,out);');
verifyEqual(testCase,summary.status,repmat("completed",4,1));
verifyEqual(testCase,numel(unique(summary.runDir)),4);
verifyFalse(testCase,isfile(fullfile(batchDir,'network.data')));
for k = 1:height(summary)
    folder = summary.runDir(k);
    verifyTrue(testCase,isfile(fullfile(folder,'network.data')));
    verifyTrue(testCase,isfile(fullfile(folder,'bond.table')));
    verifyTrue(testCase,isfile(fullfile(folder,'network.localdensity.table')));
    saved = load(fullfile(folder,'results_struct.mat'),'results');
    info = saved.results.lammps;
    verifyEqual(testCase,summary.lammpsStatus(k),string(info.status));
    verifyEqual(testCase,summary.lammpsRepulsiveRadius(k),info.rRepulsive);
    verifyEqual(testCase,info.atomCount,size(saved.results.nodes,1));
    verifyEqual(testCase,info.totalKuhn, ...
        sum(saved.results.bonds(:,3))+sum(saved.results.activeEnds(:,4)));
end
end

function cfg = smallConfig(testCase,label)
[cfg,~,~,~] = fixture(testCase,label);
cfg.Ntotal = 30; cfg.xA = .2; cfg.fA = 4;
cfg.randSeed = 17; cfg.chainModel = 'gaussian';
cfg.showLivePlot = false; cfg.makeVideo = false;
cfg.showFinalHistogram = false; cfg.saveFinalDistributions = false;
cfg.saveData = false; cfg.saveDumps = false;
end

function rows = sectionRows(text,heading,n,width)
lines = splitlines(string(text));
start = find(strtrim(lines)==heading,1)+1;
while strlength(strtrim(lines(start)))==0, start = start+1; end
rows = zeros(n,width);
for k = 1:n
    rows(k,:) = sscanf(lines(start+k-1),'%f')';
end
end
