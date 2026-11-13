function tests=testNetworkConnectivity
tests=functiontests(localfunctions);
end
function setupOnce(tc)
tc.TestData.oldPath=path;addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'SRC'));
tc.TestData.root=tempname;mkdir(tc.TestData.root);
end
function teardownOnce(tc)
path(tc.TestData.oldPath);
end
function [x,B]=triangle()
x=[1 1 0;3 1 0;2 3 0];B=[1 2 2 0 0;2 3 2 0 0;3 1 2 0 0];
end
function testTriangleAndSharedHinge(tc)
[x,B]=triangle();r=networkConnectivityMetrics(x,B,10);
verifyEqual(tc,r.N_rigid,1);verifyTrue(tc,r.wholeNetworkRigid);
verifyEqual(tc,r.N_path_y,0);verifyEqual(tc,r.rigidity.clusterWindingRank,0);
x=[x;4 4 0;5 2 0];B=[B;3 4 2 0 0;4 5 2 0 0;5 3 2 0 0];
r=networkConnectivityMetrics(x,B,10);
verifyEqual(tc,r.N_rigid,2);verifyFalse(tc,r.wholeNetworkRigid);
verifyEqual(tc,r.rigidity.multiClusterJunctionIds,3);
verifyEqual(tc,sort(r.rigidity.clusterSizes),[3;3]);
end
function testOneClusterIsNotWholeAndBarsNotCounted(tc)
[x,B]=triangle();x=[x;4 4 0];B=[B;3 4 2 0 0];
r=networkConnectivityMetrics(x,B,10);
verifyEqual(tc,r.N_rigid,1);verifyFalse(tc,r.wholeNetworkRigid);
verifyEqual(tc,r.rigidity.uncoveredJunctionIds,4);
verifyEqual(tc,r.rigidity.standaloneRigidBarCount,1);
r=networkConnectivityMetrics(x(1:2,:),[1 2 2 0 0;1 1 3 0 0],10);
verifyEqual(tc,r.N_rigid,0);verifyEqual(tc,r.selfLoopCount,1);
verifyEqual(tc,r.rigidity.standaloneRigidBarCount,1);
end
function testEmptySeed(tc)
r=networkConnectivityMetrics([1 1 0],zeros(0,5),10);
verifyEqual(tc,[r.N_path_x r.N_path_y r.N_rigid],[0 0 0]);
verifyFalse(tc,r.wholeNetworkRigid);
end
function testPeriodicFlowGaugeContourAndDoubleWinding(tc)
% Two circuits share one junction, but no bond. One winds twice vertically.
x=[1 1;2 2;3 1;4 2;5 1];
B=[1 2 2 0 0;2 3 2 0 0;3 1 2 0 1;1 4 3 0 0;4 5 3 0 0;5 1 3 0 2];
p=periodicNetworkPaths(x,B);verifyEqual(tc,p{2}.capacity,3);
verifyEqual(tc,p{2}.positiveWindingCycles,2);
verifyEqual(tc,p{2}.highlightedKuhn,15);
verifyEqual(tc,p{2}.primalDualGap,0,'AbsTol',1e-9);
% Rewrapping vertex coordinates and compensating edge images is a gauge change.
g=[1 -2;0 1;-1 0;2 2;0 -1];Bg=B;Bg(:,4:5)=B(:,4:5)+g(B(:,1),:)-g(B(:,2),:);
pg=periodicNetworkPaths(x+10*g,Bg);
verifyEqual(tc,pg{2}.capacity,p{2}.capacity);verifyEqual(tc,pg{2}.highlightedKuhn,p{2}.highlightedKuhn);
% A long parallel route shares the winding edge: optimizer must use short route.
B=[1 2 1 0 0;2 3 1 0 0;3 1 1 0 1;1 4 10 0 0;4 3 10 0 0];
p=periodicNetworkPaths(x(1:4,:),B);
verifyEqual(tc,p{2}.capacity,1);verifyEqual(tc,p{2}.highlightedKuhn,3);
end
function testFinalSelectionBatchBothModesAndSnapshots(tc)
batch=fullfile(tc.TestData.root,'batch');mkdir(batch);
[nodes,bonds]=triangle();activeEnds=zeros(0,4);
for mode={'growing','preformed'}
    dir=fullfile(batch,mode{1});mkdir(dir);mkdir(fullfile(dir,'dumps'));
    metadata=struct('config',struct('L',10,'b',1.6,'fA',4,'growthMode',mode{1}));
    save(fullfile(dir,'dumps','metadata.mat'),'metadata');
    snapshot=struct('nodes',nodes,'bonds',bonds,'activeEnds',activeEnds,'L',10,'phase','final');
    save(fullfile(dir,'dumps','snapshot_000002.mat'),'snapshot');
    snapshot.phase='usage';snapshot.bonds=zeros(0,5);snapshot.usagePercent=100;
    save(fullfile(dir,'dumps','snapshot_000003.mat'),'snapshot');
    [s,source]=loadFinalNetworkState(fullfile(dir,'dumps'));
    verifyEqual(tc,s.bonds,bonds);verifyEqual(tc,source.phase,'final');
end
opts=struct('outputRoot',tc.TestData.root,'saveSnapshots',true,'snapshotResolution',90,'continueOnError',false);
r=postprocessNetworkConnectivity(batch,opts);
verifyEqual(tc,r.completed,2);verifyEqual(tc,r.failed,0);verifyEqual(tc,r.summary.N_rigid,[1;1]);
verifyEqual(tc,sort(r.summary.growthMode),["growing";"preformed"]);
verifyTrue(tc,isfile(fullfile(r.summary.outputDirectory(1),'paths_y.png')));
verifyTrue(tc,isfile(fullfile(r.summary.outputDirectory(1),'rigid_clusters.png')));
for name={'paths_x.png','paths_y.png','rigid_clusters.png'}
    img=imfinfo(fullfile(r.summary.outputDirectory(1),name{1}));
    verifyEqual(tc,[img.Width img.Height],[360 360]);
end
verifyEqual(tc,height(readtable(fullfile(r.summary.outputDirectory(1),'rigid_cluster_nodes.csv'))),3);
end
function testIncompleteIsNotZeroAndInvalidGeometry(tc)
dir=fullfile(tc.TestData.root,'incomplete');mkdir(dir);mkdir(fullfile(dir,'dumps'));
[nodes,bonds]=triangle();snapshot=struct('nodes',nodes,'bonds',bonds,'activeEnds',zeros(0,4),'L',10,'phase','usage');
save(fullfile(dir,'dumps','snapshot_000001.mat'),'snapshot');
verifyError(tc,@()loadFinalNetworkState(dir),'NetworkGrowth:NoFinalState');
warningState=warning('off','NetworkGrowth:ConnectivityRunFailed');clean=onCleanup(@()warning(warningState));
r=postprocessNetworkConnectivity(dir,struct('outputRoot',tc.TestData.root,'saveSnapshots',false));
verifyEqual(tc,r.failed,1);verifyTrue(tc,isnan(r.summary.N_rigid));
verifyError(tc,@()networkConnectivityMetrics(nodes,[1 1 2 0 1],10),'NetworkGrowth:WindingSelfLoop');
end
