function result = networkConnectivityMetrics(nodes,bonds,L,options)
% Fixed-cell bar-rigidity clusters and unit-strand periodic winding capacity.
% Counts geometry/topology, NOT stiffness. Works for both growth chemistries.
% Bonds: [node1 node2 N_Kuhn imageX imageY], in the SAVED image convention.
if nargin<4, options=struct(); end
options=mergeConfig(struct('rankTolerance',1e-9,'pairTolerance',1e-8, ...
    'imageRange',2,'checkRankStability',true),options);
validateattributes(nodes,{'numeric'},{'2d','real','finite','nonempty'});
assert(size(nodes,2)>=2,'NetworkGrowth:ConnectivityNodes','Nodes need x and y coordinates.');
validateattributes(L,{'numeric'},{'scalar','real','finite','positive'});
if isempty(bonds),bonds=zeros(0,5);end
validateattributes(bonds,{'numeric'},{'2d','real','finite','ncols',5});
validateattributes(bonds(:,1:2),{'numeric'},{'integer','>=',1,'<=',size(nodes,1)});
validateattributes(bonds(:,3),{'numeric'},{'integer','positive'});
validateattributes(bonds(:,4:5),{'numeric'},{'integer'});
validateattributes(options.rankTolerance,{'numeric'},{'scalar','finite','positive','<',1});
validateattributes(options.pairTolerance,{'numeric'},{'scalar','finite','positive','<',1});
validateattributes(options.imageRange,{'numeric'},{'scalar','integer','>=',1});
loopRows=find(bonds(:,1)==bonds(:,2));physicalRows=find(bonds(:,1)~=bonds(:,2));
assert(~any(bonds(loopRows,4:5),'all'),'NetworkGrowth:WindingSelfLoop', ...
    'A self-edge with nonzero winding is not a primary self-loop; this analysis cannot discard it.');
B=bonds(physicalRows,:);
if ~isempty(B)
    v=nodes(B(:,2),1:2)-nodes(B(:,1),1:2)+L*B(:,4:5);
    assert(all(vecnorm(v,2,2)>1e-12*L),'NetworkGrowth:ZeroBar', ...
        'A zero-length inter-junction strand has undefined bar direction.');
end
paths=periodicNetworkPaths(nodes,B);
for a=1:2
    paths{a}.sourceBondRows=physicalRows;
    for k=1:numel(paths{a}.cycles)
        paths{a}.cycles{k}.sourceBondRows=physicalRows(paths{a}.cycles{k}.edgeIndices);
    end
end
rigidity=networkRigidClusters(nodes,B,L,options);
for k=1:numel(rigidity.clusters)
    ids=rigidity.clusters{k};rows=find(ismember(B(:,1),ids)&ismember(B(:,2),ids));
    rigidity.clusterSourceBondRows{k,1}=physicalRows(rows);
    rigidity.clusterWindingRank(k,1)=networkWindingRank(size(nodes,1),B(rows,:));
end
rigidity.wholeNetworkRigid=rigidity.rigidClusterCount==1 && ...
    numel(rigidity.clusters{1})==size(nodes,1);
result=struct('schemaVersion',1,'junctionCount',size(nodes,1),'physicalBondCount',size(B,1), ...
    'selfLoopCount',numel(loopRows),'physicalBondRows',physicalRows,'selfLoopRows',loopRows, ...
    'N_path_x',paths{1}.capacity,'N_path_y',paths{2}.capacity, ...
    'N_rigid',rigidity.rigidClusterCount,'wholeNetworkRigid',rigidity.wholeNetworkRigid, ...
    'pathsX',paths{1},'pathsY',paths{2},'rigidity',rigidity,'options',options);
end
