function tests = testNetworkDumpVisualization
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.oldPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'SRC'));
testCase.TestData.root = tempname;
mkdir(testCase.TestData.root);
fprintf('Offline visualization test artifacts: %s\n',testCase.TestData.root);
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function runDir = makeFixture(testCase,label,networkPayload)
cfg = defaultNetworkAConfig(); cfg.L = 10; cfg.b = 1; cfg.fA = 2;
cfg.Ntotal = 13; cfg.xA = 3/13; cfg.saveDumps = true;
cfg.dumpSaveNetwork = networkPayload;
runDir = fullfile(testCase.TestData.root,label); cfg.outputDir = runDir;
state = initNetworkDumps(cfg,cfg,3,10,1);
nodes = [5 5 0]; bonds = zeros(0,5); ends = repmat([5 5 1 0],2,1);
state = writeNetworkDump(state,nodes,bonds,ends,0,0,2,10,'initial');
nodes = [5 5 1;7 5 1]; bonds = [1 2 2 0 0]; ends = [5 5 1 0;7 5 2 0];
state = writeNetworkDump(state,nodes,bonds,ends,3,1,1,8,'event');
nodes = [5 5 1;7 5 2;7 8 1]; bonds = [bonds;2 3 3 0 0]; ends = [5 5 1 0;7 8 3 0];
state = writeNetworkDump(state,nodes,bonds,ends,7,2.5,0,5,'event');
nodes(3,2) = 7.5; ends(2,2) = 7.5;
writeNetworkDump(state,nodes,bonds,ends,7,2.5,0,5,'final','no_allowed_events');
end

function testOfflineMovieAndPlots(testCase)
runDir = makeFixture(testCase,'full_payload',true);
before = dir(fullfile(runDir,'dumps','*.mat'));
options = struct('videoFrameEvery',2,'videoFPS',4,'videoHoldLastSeconds',.5, ...
    'heatGrid',30,'showComparisonFigure',false,'comparisonUsagePercent',[0 50 100]);
report = postprocessNetworkDumps(runDir,options);
verifyTrue(testCase,isfile(report.videoFile));
verifyTrue(testCase,isfile(report.comparisonFile));
verifyTrue(testCase,isfile(report.finalFile));
verifyEqual(testCase,numel(report.checkpointFiles),3); % shared states for crossed targets
verifyTrue(testCase,all(isfile(report.checkpointFiles)));
index = readtable(fullfile(report.outputDir,'video_frames.csv'),'TextType','string');
verifyEqual(testCase,index.file,["snapshot_000001.mat";"snapshot_000003.mat";"snapshot_000004.mat"]);
verifyEqual(testCase,index.event,[0;7;7]);
verifyEqual(testCase,index.ssaTime,[0;2.5;2.5]);
verifyEqual(testCase,index.firstVideoFrame,[1;2;3]);
verifyEqual(testCase,index.lastVideoFrame,[1;2;5]);
video = VideoReader(report.videoFile); count = 0;
while hasFrame(video)
    rgb = readFrame(video); count = count+1;
    if count == 1, firstFrame = rgb; end
end
verifyEqual(testCase,count,5);
verifyEqual(testCase,mod([size(rgb,1) size(rgb,2)],2),[0 0]);
verifyNotEqual(testCase,firstFrame,rgb);
% The reservoir panel must contain actual red/blue curves, not just empty axes.
right = rgb(:,ceil(size(rgb,2)/2):end,:);
red = right(:,:,1)>150 & right(:,:,2)<100 & right(:,:,3)<100;
blue = right(:,:,3)>150 & right(:,:,1)<100 & right(:,:,2)<100;
verifyGreaterThan(testCase,nnz(red),200);
verifyGreaterThan(testCase,nnz(blue),200);
imwrite(rgb,fullfile(report.outputDir,'video_last_frame_test.png'));
after = dir(fullfile(runDir,'dumps','*.mat'));
verifyEqual(testCase,{after.name},{before.name});
verifyEqual(testCase,[after.bytes],[before.bytes]);
verifyEqual(testCase,[after.datenum],[before.datenum]);
fprintf('Rendered video and plots: %s\n',report.outputDir);
end

function testFilterAndDisabledExports(testCase)
runDir = makeFixture(testCase,'filter',true);
[s,~,index] = loadNetworkDumps(runDir,[5 50 75 100]);
verifyEqual(testCase,numel(s),2);
verifyEqual(testCase,index.file,["snapshot_000002.mat";"snapshot_000003.mat"]);
verifyFalse(testCase,any(index.phase=="final"));
options = struct('makeVideo',false,'saveCheckpointPlots',false, ...
    'saveComparison',false,'saveFinalPlot',false,'showComparisonFigure',false);
report = postprocessNetworkDumps(fullfile(runDir,'dumps'),options);
verifyEmpty(testCase,report.videoFile);
verifyEmpty(testCase,report.checkpointFiles);
verifyEmpty(testCase,report.finalFile);
verifyEmpty(testCase,report.comparisonFile);
verifyTrue(testCase,isfile(fullfile(report.outputDir,'visualization_report.mat')));
end

function testDistributionsWithoutNetwork(testCase)
runDir = makeFixture(testCase,'distributions_only',false);
options = struct('makeVideo',true,'saveCheckpointPlots',false, ...
    'saveFinalPlot',false,'showComparisonFigure',false,'comparisonUsagePercent',[0 50 100]);
warningState = warning('off','NetworkGrowth:NoNetworkVideo');
cleanup = onCleanup(@()warning(warningState));
report = postprocessNetworkDumps(runDir,options);
verifyEmpty(testCase,report.videoFile);
verifyTrue(testCase,isfile(report.comparisonFile));
end

function testSharedHistogramRangesAndPeakCounts(testCase)
% A transient large extension and a concentrated count peak must both fit,
% even though neither occurs in the last state.
early = struct('completed',struct('kuhnCount',[],'endToEnd',[]));
middle = struct('completed',struct('kuhnCount',[2;2;2;8],'endToEnd',[1;1;1;12]));
late = struct('completed',struct('kuhnCount',[2;3;4;8],'endToEnd',[1;2;3;4]));
style = networkDumpHistogramStyle({early,middle,late});
verifyEqual(testCase,style.kuhnXLim,[0 9]);
verifyGreaterThan(testCase,style.endToEndXLim(2),12);
verifyEqual(testCase,style.histogramYLim,[0 .83]);
verifyEqual(testCase,diff(style.nEdges),ones(1,numel(style.nEdges)-1));
manual = networkDumpHistogramStyle({early,middle,late}, ...
    struct('kuhnXLim',[0 20],'endToEndXLim',[0 15],'histogramYLim',[0 50],'endToEndBins',15));
verifyEqual(testCase,manual.kuhnXLim,[0 20]);
verifyEqual(testCase,manual.endToEndXLim,[0 15]);
verifyEqual(testCase,manual.histogramYLim,[0 50]);
verifyEqual(testCase,numel(manual.rEdges),16);
end

function testSnapshotAxesStayFixedAndRdfIsAbsent(testCase)
runDir = makeFixture(testCase,'fixed_axes',true);
options = struct('makeVideo',false,'saveComparison',false, ...
    'checkpointUsagePercent',[0 100],'showComparisonFigure',false);
report = postprocessNetworkDumps(runDir,options);
style = report.histogramStyle;
[snapshots,~,index] = loadNetworkDumps(runDir);
for k = 1:numel(snapshots)
    % File input must inherit run-wide limits, even for an empty initial state.
    fig = plotNetworkDump(fullfile(runDir,'dumps',index.file(k)),struct('visible','off'));
    cleanup = onCleanup(@()close(fig));
    verifyEqual(testCase,numel(findobj(fig,'Type','axes')),5);
    xyAxis = findobj(fig,'Tag','displacementHistogram');
    verifyEqual(testCase,xyAxis.XLim,style.displacementLim);
    verifyEqual(testCase,xyAxis.YLim,style.displacementLim);
    verifyEqual(testCase,xyAxis.CLim,style.displacementCLim);
    nAxis = findobj(fig,'Tag','kuhnHistogram');
    rAxis = findobj(fig,'Tag','weightFractionHistogram');
    verifyEqual(testCase,nAxis.XLim,style.kuhnXLim);
    verifyEqual(testCase,rAxis.XLim,style.weightXLim);
    verifyEqual(testCase,nAxis.YLim,style.kuhnYLim);
    verifyEqual(testCase,rAxis.YLim,style.weightYLim);
    nHist = findobj(nAxis,'Type','histogram');
    rHist = findobj(rAxis,'Type','histogram');
    verifyEqual(testCase,nHist.BinEdges,style.nEdges);
    verifyEqual(testCase,rHist.BinEdges,style.weightEdges);
    verifyEqual(testCase,nHist.Normalization,'probability');
    verifyEqual(testCase,rHist.Normalization,'count'); % precomputed fractions, no second normalization
    if ~isempty(snapshots{k}.completed.kuhnCount)
        verifyEqual(testCase,sum(nHist.Values),1,'AbsTol',1e-12);
        verifyEqual(testCase,sum(rHist.Values),1,'AbsTol',1e-12);
    end
    if k == 1 || k == numel(snapshots)
        exportgraphics(fig,fullfile(report.outputDir,sprintf('fixed_axes_state_%d_test.png',k)),'Resolution',120);
    end
    clear cleanup;
end
verifyTrue(testCase,isfile(report.finalFile));
verifyEqual(testCase,numel(report.checkpointFiles),2);
fprintf('Fixed-axis plot artifacts: %s\n',report.outputDir);
end

function testNumberAndWeightFractionsUseDifferentDenominators(testCase)
runDir = makeFixture(testCase,'probability',true);
states = loadNetworkDumps(runDir); s = states{end};
s.completed.kuhnCount = [1;1;1;4];
s.completed.endToEnd = [1;2;2;3;3;3];
fig = plotNetworkDump(s,struct('visible','off'));
cleanup = onCleanup(@()close(fig));
hn = findobj(findobj(fig,'Tag','kuhnHistogram'),'Type','histogram');
hr = findobj(findobj(fig,'Tag','weightFractionHistogram'),'Type','histogram');
verifyEqual(testCase,hn.Values,histcounts(s.completed.kuhnCount,hn.BinEdges)/4,'AbsTol',1e-12);
verifyEqual(testCase,hr.Values(hr.BinEdges(1:end-1)==.5),3/7,'AbsTol',1e-12);
verifyEqual(testCase,hr.Values(hr.BinEdges(1:end-1)==3.5),4/7,'AbsTol',1e-12);
verifyEqual(testCase,sum(hn.Values),1,'AbsTol',1e-12);
verifyEqual(testCase,sum(hr.Values),1,'AbsTol',1e-12);
% Cropping must not renormalize away out-of-range samples.
fig2 = plotNetworkDump(s,struct('visible','off','weightXLim',[0 2]));
cleanup2 = onCleanup(@()close(fig2));
hr2 = findobj(findobj(fig2,'Tag','weightFractionHistogram'),'Type','histogram');
verifyEqual(testCase,sum(hr2.Values),3/7,'AbsTol',1e-12);
end

function testWeightFractionsGroupedBinsAndEmptyPopulation(testCase)
[w,total] = chainWeightFractions([1;2;4;4],[.5 2.5 5.5]);
verifyEqual(testCase,total,11);
verifyEqual(testCase,w,[3/11 8/11],'AbsTol',1e-12);
verifyEqual(testCase,chainWeightFractions([1;2;4;4],[.5 2.5]),3/11,'AbsTol',1e-12);
verifyEqual(testCase,chainWeightFractions([4;4],[.5 2.5]),0);
[w,total] = chainWeightFractions([],[-.5 .5 1.5]);
verifyEqual(testCase,w,[0 0]); verifyEqual(testCase,total,0);
[w,total] = chainWeightFractions([0;0],[-.5 .5 1.5]);
verifyEqual(testCase,w,[0 0]); verifyEqual(testCase,total,0);
s = struct('completed',struct('kuhnCount',[1;1;1;3],'endToEnd',ones(4,1)));
style = networkDumpHistogramStyle(s);
verifyEqual(testCase,style.kuhnYLim,[0 .83]);
verifyEqual(testCase,style.weightYLim,[0 .55]);
end

function testPeriodicDisplacementHeatmapAndSavedVectors(testCase)
nodes = [1 2 2;9 5 1;4 1 1];
bonds = [1 2 3 -1 0;1 3 4 0 1];
legacy = struct('nodes',nodes,'bonds',bonds,'L',10);
expected = [-2 3;3 9];
[dv,available] = networkDumpDisplacements(legacy);
verifyTrue(testCase,available); verifyEqual(testCase,dv,expected);
translated = legacy; translated.nodes(:,1:2) = translated.nodes(:,1:2)+[100 -40];
verifyEqual(testCase,networkDumpDisplacements(translated),expected);
reversed = legacy; reversed.bonds(:,1:2) = bonds(:,[2 1]); reversed.bonds(:,4:5) = -bonds(:,4:5);
verifyEqual(testCase,networkDumpDisplacements(reversed),-expected);

cfg = defaultNetworkAConfig(); cfg.L = 10; cfg.saveDumps = true; cfg.dumpSaveNetwork = false;
cfg.outputDir = fullfile(testCase.TestData.root,'vector_payload');
state = initNetworkDumps(cfg,cfg,3,20,1);
writeNetworkDump(state,nodes,bonds,zeros(0,4),1,1,0,13,'event');
states = loadNetworkDumps(cfg.outputDir); s = states{1};
verifyFalse(testCase,isfield(s,'nodes'));
verifyEqual(testCase,s.completed.endToEndVector,expected);
verifyEqual(testCase,s.completed.endToEnd,sqrt(sum(expected.^2,2)));
fig = plotNetworkDump(s,struct('visible','off','displacementBins',10));
cleanup = onCleanup(@()close(fig));
ax = findobj(fig,'Tag','displacementHistogram');
im = findobj(ax,'Type','image');
style = networkDumpHistogramStyle(s,struct('displacementBins',10));
prob = histcounts2(expected(:,1),expected(:,2),style.displacementEdges,style.displacementEdges)/2;
verifyEqual(testCase,im.CData,prob');
verifyEqual(testCase,sum(im.CData,'all'),1);
verifyEqual(testCase,im.AlphaData,double(prob'>0));
verifyEqual(testCase,ax.Color,[1 1 1]);
verifyEqual(testCase,ax.YDir,'normal');
fig2 = plotNetworkDump(s,struct('visible','off','displacementLimit',4));
cleanup2 = onCleanup(@()close(fig2));
im2 = findobj(findobj(fig2,'Tag','displacementHistogram'),'Type','image');
verifyEqual(testCase,sum(im2.CData,'all'),.5);
oldDistributionOnly = s; oldDistributionOnly.completed = rmfield(s.completed,'endToEndVector');
[~,available] = networkDumpDisplacements(oldDistributionOnly);
verifyFalse(testCase,available);
end

function testRobustDisplacementColorCapPreservesProbability(testCase)
% A lone early chain must not hide a later broad population. Only colors clip.
[x,y] = meshgrid(-4.5:1:4.5);
early = struct('event',1,'ssaTime',1,'phase','event','usagePercent',5, ...
    'completed',struct('kuhnCount',4,'endToEnd',sqrt(.5),'endToEndVector',[.5 .5]));
late = early; late.completed.kuhnCount = 4*ones(100,1);
late.completed.endToEndVector = [x(:) y(:)];
late.completed.endToEnd = hypot(x(:),y(:));
opts = struct('displacementBins',10,'displacementLimit',5);
style = networkDumpHistogramStyle({early,late},opts);
verifyEqual(testCase,style.displacementCLim,[0 .01],'AbsTol',1e-12);
opts.displacementColorPercentile = 100;
full = networkDumpHistogramStyle({early,late},opts);
verifyEqual(testCase,full.displacementCLim,[0 1]);
opts.displacementCLim = [0 .02];
manual = networkDumpHistogramStyle({early,late},opts);
verifyEqual(testCase,manual.displacementCLim,[0 .02]);
for s={early,late}
    fig = plotNetworkDump(s{1},struct('visible','off','histogramStyle',style));
    cleanup = onCleanup(@()close(fig));
    ax = findobj(fig,'Tag','displacementHistogram'); im = findobj(ax,'Type','image');
    verifyEqual(testCase,ax.CLim,[0 .01],'AbsTol',1e-12);
    verifyEqual(testCase,sum(im.CData,'all'),1,'AbsTol',1e-12);
    if numel(s{1}.completed.kuhnCount)==1
        verifyEqual(testCase,max(im.CData,[],'all'),1); % never clamp stored probabilities
        verifyEqual(testCase,nnz(im.AlphaData),1); % all empty bins transparent
        cb = findobj(fig,'Type','colorbar');
        verifyTrue(testCase,any(contains(string(cb.TickLabels),'\geq')));
    end
    clear cleanup;
end
empty = early; empty.completed.endToEndVector=zeros(0,2);
empty.completed.kuhnCount=[]; empty.completed.endToEnd=[];
emptyStyle=networkDumpHistogramStyle(empty);
verifyEqual(testCase,emptyStyle.displacementCLim,[0 .001]);
end
