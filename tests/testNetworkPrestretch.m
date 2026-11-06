function tests = testNetworkPrestretch
tests = functiontests(localfunctions);
end

function setupOnce(tc)
tc.TestData.oldPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'SRC'));
tc.TestData.root = tempname; mkdir(tc.TestData.root);
fprintf('Prestretch test artifacts: %s\n',tc.TestData.root);
end

function teardownOnce(tc)
path(tc.TestData.oldPath);
end

function s = snapshot()
s = struct('nodes',[1 1 1;9 1 2;3 4 3], ...
    'bonds',[1 2 2 -1 0;2 3 2 1 0;3 3 3 0 0], ...
    'L',10,'b',2,'nodeCount',3,'bondCount',3, ...
    'usagePercent',100,'phase','final','event',12,'ssaTime',4);
s.completed = struct('kuhnCount',[2;2;3],'contourLength',[4;4;6], ...
    'endToEnd',[2;5;0]);
end

function testLegacyMagnitudesAndContourFallback(tc)
s = snapshot(); s = rmfield(s,{'nodes','bonds','L'});
s.dangling = struct('kuhnCount',100,'contourLength',200);
[lambda,available] = networkDumpPrestretch(s);
verifyTrue(tc,available);
verifyEqual(tc,lambda,[.5;1.25;0]); % self-loop included; no dangling sample
s = rmfield(s,'b');
verifyEqual(tc,networkDumpPrestretch(s),lambda);
s.completed = rmfield(s.completed,'contourLength');
[lambda,available] = networkDumpPrestretch(s);
verifyFalse(tc,available); verifyEmpty(tc,lambda);
end

function testStoredImagesAndNetworkOnly(tc)
s = snapshot(); s = rmfield(s,'completed');
verifyEqual(tc,networkDumpPrestretch(s),[.5;1.25;0]);
shifted = s; shifted.nodes(:,1:2) = shifted.nodes(:,1:2)+[100 -30];
verifyEqual(tc,networkDumpPrestretch(shifted),[.5;1.25;0]);
reversed = s; reversed.bonds(:,1:2) = s.bonds(:,[2 1]);
reversed.bonds(:,4:5) = -s.bonds(:,4:5);
verifyEqual(tc,networkDumpPrestretch(reversed),[.5;1.25;0]);
% The second stored vector is (4,3), not the vector obtained by ignoring its image.
verifyEqual(tc,networkDumpPrestretch(s)*s.b.*s.bonds(:,3),[2;5;0]);
end

function testEmptyInvalidAndMismatchedSamples(tc)
s = struct('b',2,'completed',struct('kuhnCount',[],'endToEnd',[]));
[lambda,available] = networkDumpPrestretch(s);
verifyTrue(tc,available); verifyEmpty(tc,lambda);
s.completed.kuhnCount = [0;2;NaN;1]; s.completed.endToEnd = [1;4;5;Inf];
verifyEqual(tc,networkDumpPrestretch(s),1);
s.completed.endToEnd = [1;2];
[lambda,available] = networkDumpPrestretch(s);
verifyFalse(tc,available); verifyEmpty(tc,lambda);
end

function testSharedAxesAndCroppedNormalization(tc)
middle = snapshot(); late = middle;
late.completed.endToEnd = [1;1;0];
opts = struct('prestretchBins',5);
style = networkDumpHistogramStyle({middle,late},opts);
verifyGreaterThan(tc,style.prestretchXLim(2),1.25);
verifyEqual(tc,numel(style.prestretchEdges),6);
for s = {middle,late}
    fig = plotNetworkDump(s{1},struct('visible','off','histogramStyle',style));
    cleanup = onCleanup(@()close(fig));
    ax = findobj(fig,'Tag','prestretchHistogram'); h = findobj(ax,'Type','histogram');
    verifyEqual(tc,ax.XLim,style.prestretchXLim);
    verifyEqual(tc,ax.YLim,style.prestretchYLim);
    verifyEqual(tc,h.BinEdges,style.prestretchEdges);
    verifyEqual(tc,sum(h.Values),1,'AbsTol',1e-12);
    clear cleanup;
end
fig = plotNetworkDump(middle,struct('visible','off','prestretchXLim',[0 .75], ...
    'prestretchYLim',[0 .8],'prestretchBins',3));
cleanup = onCleanup(@()close(fig));
ax = findobj(fig,'Tag','prestretchHistogram'); h = findobj(ax,'Type','histogram');
verifyEqual(tc,sum(h.Values),2/3,'AbsTol',1e-12);
verifyEqual(tc,ax.XLim,[0 .75]); verifyEqual(tc,ax.YLim,[0 .8]);
end

function testPostprocessExistingDumpsAndComparison(tc)
runDir = makeDumps(tc,'legacy');
% Simulate an older distribution-only run with no stored vector data.
files = dir(fullfile(runDir,'dumps','snapshot_*.mat'));
for k = 1:numel(files)
    file = fullfile(files(k).folder,files(k).name);
    loaded = load(file,'snapshot'); snapshot = loaded.snapshot;
    snapshot.completed = rmfield(snapshot.completed,'endToEndVector');
    save(file,'snapshot');
end
before = dir(fullfile(runDir,'dumps','*.mat'));
opts = struct('makeVideo',false,'showComparisonFigure',false, ...
    'checkpointUsagePercent',[0 50 100],'comparisonUsagePercent',[0 50 100], ...
    'prestretchBins',10);
report = postprocessNetworkDumps(runDir,opts);
verifyTrue(tc,isfile(report.finalFile)); verifyTrue(tc,isfile(report.comparisonFile));
verifyEqual(tc,numel(report.checkpointFiles),3);
verifyEqual(tc,numel(report.histogramStyle.prestretchEdges),11);
index = readtable(fullfile(report.outputDir,'comparison_index.csv'));
verifyEqual(tc,index.prestretchSampleCount,[0;1;2]);
verifyTrue(tc,all(index.prestretchAvailable));
[fig,selection] = plotNetworkDumps(runDir,struct('visible','off', ...
    'usagePercent',[0 50 100],'population','dangling','prestretchEdges',[0 .5 1 1.5]));
cleanup = onCleanup(@()close(fig));
verifyEqual(tc,selection.prestretchSampleCount,[0;1;2]);
ax = findobj(fig,'Tag','prestretchComparison');
curves = findobj(ax,'Type','stair');
verifyEqual(tc,numel(curves),2); % initial empty state gets only a legend entry
for k = 1:numel(curves)
    verifyEqual(tc,sum(curves(k).YData(1:end-1)),1,'AbsTol',1e-12);
end
after = dir(fullfile(runDir,'dumps','*.mat'));
verifyEqual(tc,{after.name},{before.name});
verifyEqual(tc,[after.bytes],[before.bytes]);
verifyEqual(tc,[after.datenum],[before.datenum]);
fprintf('Prestretch plots: %s\n',report.outputDir);
end

function runDir = makeDumps(tc,label)
cfg = defaultNetworkAConfig(); cfg.L = 10; cfg.b = 2; cfg.fA = 4;
cfg.saveDumps = true; cfg.dumpSaveNetwork = false;
runDir = fullfile(tc.TestData.root,label); cfg.outputDir = runDir;
state = initNetworkDumps(cfg,cfg,3,20,1);
nodes = [1 1 0]; bonds = zeros(0,5); ends = [1 1 1 0];
state = writeNetworkDump(state,nodes,bonds,ends,0,0,2,20,'initial');
nodes = [1 1 1;9 1 1]; bonds = [1 2 2 -1 0];
state = writeNetworkDump(state,nodes,bonds,ends,3,1,1,18,'event');
nodes = [nodes;3 4 1]; bonds = [bonds;2 3 2 1 0];
state = writeNetworkDump(state,nodes,bonds,ends,7,2,0,16,'event');
nodes(3,2) = 3;
writeNetworkDump(state,nodes,bonds,ends,7,2,0,16,'final','no_allowed_events');
end
