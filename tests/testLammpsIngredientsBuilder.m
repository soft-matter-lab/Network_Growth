function tests = testLammpsIngredientsBuilder
tests = functiontests(localfunctions);
end

function setupOnce(tc)
tc.TestData.oldPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'SRC'));
tc.TestData.root = tempname; mkdir(tc.TestData.root);
tc.TestData.warning = warning('off','NetworkGrowth:LammpsSelfLoops');
fprintf('Builder test artifacts: %s\n',tc.TestData.root);
end

function teardownOnce(tc)
path(tc.TestData.oldPath); warning(tc.TestData.warning);
end

function [folder,r,opts] = fixture(tc,label)
folder = fullfile(tc.TestData.root,label); mkdir(folder);
cfg = struct('L',10,'b',2,'Ntotal',100,'xA',.5,'outputDir','old/moved/run', ...
    'lammpsIncludeDanglingKuhn',false); % stale settings must not exclude dangling segments
r = struct('nodes',[1 1 3;3 1 2;3 4 1], ...
    'bonds',[1 2 2 0 0;2 3 3 0 0;1 1 4 0 0], ...
    'activeEnds',[3 4 3 5;1 1 1 0],'config',cfg);
opts = struct('exportConfig',struct('lammpsNRho',21));
end

function testLegacyResultsAndInputPreservation(tc)
[folder,results,opts] = fixture(tc,'results');
file = fullfile(folder,'results_struct.mat'); save(file,'results');
% A competing raw final file must not override the preferred results struct.
nodes = results.nodes; bonds = results.bonds; bonds(:,3) = 50;
activeEnds = results.activeEnds; cfg = results.config;
save(fullfile(folder,'Network_PartA_v2_Cleaned_Data.mat'),'nodes','bonds','activeEnds','cfg');
before = dir(file); beforeData = filereadBytes(file);
info = buildLammpsIngredients(folder,opts);
verifyEqual(tc,info.sourceKind,'results'); verifyEqual(tc,info.sourceFile,file);
verifyEqual(tc,info.sourcePhase,'final'); verifyEqual(tc,info.atomCount,3);
verifyEqual(tc,info.bondCount,2); verifyEqual(tc,info.omittedSelfLoopRows,3);
verifyEqual(tc,info.totalKuhn,14); verifyEqual(tc,info.danglingKuhn,5);
verifyEqual(tc,info.rRepulsive,sqrt(14/3),'AbsTol',1e-12);
verifyEqual(tc,info.outputDir,folder); verifyTrue(tc,isfile(info.reportFile));
saved = load(info.reportFile,'info'); verifyEqual(tc,saved.info,info);
verifyEqual(tc,filereadBytes(file),beforeData);
after = dir(file); verifyEqual(tc,after.datenum,before.datenum);
end

function testRawMatWithCustomName(tc)
[folder,r,opts] = fixture(tc,'raw');
nodes = r.nodes; bonds = r.bonds; activeEnds = r.activeEnds; cfg = r.config;
file = fullfile(folder,'old_custom_final.mat');
save(file,'nodes','bonds','activeEnds','cfg');
info = buildLammpsIngredients(folder,opts);
verifyEqual(tc,info.sourceKind,'raw'); verifyEqual(tc,info.sourceFile,file);
verifyEqual(tc,info.totalKuhn,14); verifyEqual(tc,info.config.b,2);
end

function testFinalDumpAndMetadataFallback(tc)
[folder,r,opts] = fixture(tc,'dumps'); dumpDir = fullfile(folder,'dumps'); mkdir(dumpDir);
metadata = struct('config',r.config); save(fullfile(dumpDir,'metadata.mat'),'metadata');
snapshot = rmfield(r,'config'); snapshot.phase = 'final';
save(fullfile(dumpDir,'snapshot_000002.mat'),'snapshot');
snapshot.phase = 'event'; snapshot.bonds(:,3) = 20;
save(fullfile(dumpDir,'snapshot_000003.mat'),'snapshot');
opts.allowLatestSnapshot = true;
info = buildLammpsIngredients(dumpDir,opts);
verifyEqual(tc,info.sourceKind,'snapshot'); verifyEqual(tc,info.sourcePhase,'final');
verifyEqual(tc,info.sourceFile,fullfile(dumpDir,'snapshot_000002.mat'));
verifyEqual(tc,info.totalKuhn,14); verifyEqual(tc,info.outputDir,folder);
verifyFalse(tc,isfile(fullfile(dumpDir,'network.data')));
end

function testUnfinishedRunRequiresOptIn(tc)
[folder,r,opts] = fixture(tc,'checkpoint'); dumpDir = fullfile(folder,'dumps'); mkdir(dumpDir);
snapshot = rmfield(r,'config'); snapshot.L = 10; snapshot.b = 2; snapshot.phase = 'event';
save(fullfile(dumpDir,'snapshot_000001.mat'),'snapshot');
snapshot.activeEnds(:,4) = [7;0];
save(fullfile(dumpDir,'snapshot_000002.mat'),'snapshot');
verifyError(tc,@()buildLammpsIngredients(folder,opts),'NetworkGrowth:LammpsBuilderNoFinal');
verifyFalse(tc,isfile(fullfile(folder,'network.data')));
opts.allowLatestSnapshot = true;
old = warning('off','NetworkGrowth:LammpsBuilderCheckpoint'); cleanup = onCleanup(@()warning(old));
info = buildLammpsIngredients(folder,opts);
verifyEqual(tc,info.sourcePhase,'event'); verifyEqual(tc,info.totalKuhn,16);
verifyTrue(tc,endsWith(info.sourceFile,'snapshot_000002.mat'));
verifyTrue(tc,startsWith(fileread(info.dataFile),'Network Growth saved state'));
end

function testMissingGeometryAndDistributionOnly(tc)
[folder,results,opts] = fixture(tc,'missing_geometry'); results = rmfield(results,'config');
save(fullfile(folder,'results_struct.mat'),'results');
verifyError(tc,@()buildLammpsIngredients(folder,opts),'NetworkGrowth:LammpsBuilderMissingData');
verifyFalse(tc,isfile(fullfile(folder,'network.data')));
[folder,~,opts] = fixture(tc,'distribution_only'); mkdir(fullfile(folder,'dumps'));
snapshot = struct('phase','final','L',10,'b',2, ...
    'completed',struct('kuhnCount',[2;3;4],'endToEnd',[2;3;0]));
save(fullfile(folder,'dumps','snapshot_000002.mat'),'snapshot');
verifyError(tc,@()buildLammpsIngredients(folder,opts),'NetworkGrowth:LammpsBuilderMissingData');
verifyFalse(tc,isfile(fullfile(folder,'network.data')));
end

function testExplicitSourceDestinationAndOverwriteProtection(tc)
[folder,results,opts] = fixture(tc,'explicit');
save(fullfile(folder,'alternative.mat'),'results');
opts.sourceFile = 'alternative.mat'; opts.outputDir = fullfile(folder,'lammps_rebuild');
info = buildLammpsIngredients(folder,opts);
verifyTrue(tc,isfile(fullfile(opts.outputDir,'network.data')));
verifyFalse(tc,isfile(fullfile(folder,'network.data')));
before = fileread(info.dataFile);
verifyError(tc,@()buildLammpsIngredients(folder,opts),'NetworkGrowth:ExistingLammpsFiles');
verifyEqual(tc,fileread(info.dataFile),before);
opts.exportConfig.b = 99;
verifyError(tc,@()buildLammpsIngredients(folder,opts),'NetworkGrowth:LammpsBuilderOptions');
end

function testAmbiguousRawFilesNeedSelection(tc)
[folder,r,opts] = fixture(tc,'ambiguous');
nodes = r.nodes; bonds = r.bonds; activeEnds = r.activeEnds; cfg = r.config;
save(fullfile(folder,'old_a.mat'),'nodes','bonds','activeEnds','cfg');
save(fullfile(folder,'old_b.mat'),'nodes','bonds','activeEnds','cfg');
verifyError(tc,@()buildLammpsIngredients(folder,opts),'NetworkGrowth:LammpsBuilderAmbiguous');
opts.sourceFile = 'old_b.mat';
info = buildLammpsIngredients(folder,opts);
verifyTrue(tc,endsWith(info.sourceFile,'old_b.mat'));
end

function bytes = filereadBytes(file)
fid = fopen(file,'rb'); cleanup = onCleanup(@()fclose(fid));
bytes = fread(fid,Inf,'*uint8');
end
