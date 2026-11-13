function out=networkRigidClusters(nodes,B,L,options)
% Port of the certified figure analysis, now entirely in MATLAB.
% Fixed cell, one-cell-periodic first-order bar constraints. No force model.
V=size(nodes,1);E=size(B,1);x=nodes(:,1:2);p=B(:,1);q=B(:,2);
rv=x(q,:)-x(p,:)+B(:,4:5)*L;unit=rv./vecnorm(rv,2,2);
R=sparse(repmat((1:E)',4,1),[p;q;p+V;q+V], ...
    [-unit(:,1);unit(:,1);-unit(:,2);unit(:,2)],E,2*V);
if E==0
    values=[];Z=eye(2*V);rankR=0;counts=zeros(1,4);
else
    [~,S,Z]=svd(full(R));d=min(size(S));values=diag(S(1:d,1:d));
    rankR=nnz(values>options.rankTolerance*max(values));
    counts=arrayfun(@(t)nnz(values>t*max(values)),[1e-7,1e-8,1e-9,1e-10]);
end
if options.checkRankStability && any(counts~=rankR)
    error('NetworkGrowth:RigidityRankSensitive','Rigidity rank changes over 1e-7..1e-10; no cluster claim is made.');
end
internal=2*V-2-rankR;Q=Z(:,rankR+1:end);
tx=[ones(V,1);zeros(V,1)]/sqrt(V);ty=[zeros(V,1);ones(V,1)]/sqrt(V);
Q=Q-tx*(tx'*Q)-ty*(ty'*Q);[Q,D]=svd(Q,'econ');Q=Q(:,diag(D)>1e-6);
assert(size(Q,2)==internal && norm(R*Q,'fro')<1e-7, ...
    'NetworkGrowth:RigidityNullspace','Could not certify the periodic mechanism basis.');
clusters={};bars={};certificates=struct('size',{},'rigidMotionResidual',{},'internalRank',{});
if internal==0 && V>=3
    clusters={(1:V)'};
elseif E>0
    wrapped=mod(x,L);range=-options.imageRange:options.imageRange;
    [sy,sx]=ndgrid(range,range);shifts=[sx(:) sy(:)];
    adj=false(V);imageX=zeros(V);imageY=zeros(V);
    for i=1:V
        for j=i+1:V
            dv=wrapped(j,:)-wrapped(i,:)+L*shifts;
            du=[Q(j,:)-Q(i,:);Q(j+V,:)-Q(i+V,:)];motionNorm=norm(du,'fro');
            residual=vecnorm(dv*du,2,2)./max(vecnorm(dv,2,2),1e-12);
            accepted=find(residual<options.pairTolerance*motionNorm | motionNorm<1e-10);
            if ~isempty(accepted)
                [~,choice]=min(vecnorm(dv(accepted,:),2,2));im=shifts(accepted(choice),:);
                adj(i,j)=true;adj(j,i)=true;
                imageX(i,j)=im(1);imageX(j,i)=-im(1);imageY(i,j)=im(2);imageY(j,i)=-im(2);
            end
        end
    end
    assert(all(adj(sub2ind([V V],p,q))),'NetworkGrowth:RigidityImageRange', ...
        'A saved bar is missing from the closure. Increase imageRange or inspect tolerances.');
    cliques=maximalCliques(adj);
    [~,order]=sort(cellfun(@numel,cliques),'descend');cliques=cliques(order);
    actual=sparse([p;q],[q;p],true,V,V);
    for ci=1:numel(cliques)
        c=sort(cliques{ci}(:));n=numel(c);
        if n<3
            if n==2 && actual(c(1),c(2)),bars{end+1,1}=c;end
            continue;
        end
        anchor=c(1);coords=wrapped(c,:)+L*[imageX(anchor,c)' imageY(anchor,c)'];
        centered=coords-mean(coords,1);
        assert(rank(centered,1e-7)==2,'NetworkGrowth:DegenerateRigidCluster', ...
            'A candidate rigid cluster is collinear; cannot certify a multi-junction rigid body.');
        design=[ones(n,1) zeros(n,1) -centered(:,2);zeros(n,1) ones(n,1) centered(:,1)];
        target=[Q(c,:);Q(c+V,:)];residual=norm(design*(design\target)-target,'fro');
        assert(residual<2e-5,'NetworkGrowth:RigidImages', ...
            'Candidate cluster images do not define a single first-order rigid body.');
        inside=ismember(p,c)&ismember(q,c);localR=full(R(inside,[c;c+V]));
        sv=svd(localR);localRank=nnz(sv>1e-8*max(sv));
        assert(localRank>=2*n-3,'NetworkGrowth:ExternalRigidConstraints', ...
            'A candidate subset relies on external constraints; refusing an uncertified cluster count.');
        clusters{end+1,1}=c;
        certificates(end+1)=struct('size',n,'rigidMotionResidual',residual,'internalRank',localRank);
    end
end
membership=zeros(V,1);
for i=1:numel(clusters),membership(clusters{i})=membership(clusters{i})+1;end
uncovered=find(membership==0);uncoveredSizes=[];
if ~isempty(uncovered)
    g=graph(p,q,[],V);components=conncomp(subgraph(g,uncovered));
    uncoveredSizes=accumarray(components(:),1);
end
sizes=cellfun(@numel,clusters);largest=max([0;sizes(:)]);
out=struct('rigidClusterCount',numel(clusters),'clusters',{clusters}, ...
    'clusterSizes',sizes(:),'clusterSourceBondRows',{{}},'clusterWindingRank',zeros(0,1), ...
    'largestRigidClusterFraction',largest/V,'coveredJunctionCount',nnz(membership), ...
    'junctionClusterMultiplicity',membership,'multiClusterJunctionIds',find(membership>1), ...
    'uncoveredJunctionIds',uncovered,'uncoveredRegionSizes',uncoveredSizes, ...
    'uncoveredRegionCount',numel(uncoveredSizes),'standaloneRigidBars',{bars}, ...
    'standaloneRigidBarCount',numel(bars),'internalInfinitesimalBarMechanisms',internal, ...
    'statesOfSelfStress',E-rankR,'rank',rankR,'rankToleranceSweep',[1e-7 1e-8 1e-9 1e-10;counts], ...
    'nullspaceResidual',norm(R*Q,'fro'),'clusterCertificates',certificates, ...
    'rigidClusterCertificatesPassed',true,'wholeNetworkRigid',false);
end

function cliques=maximalCliques(adj)
% Iterative Bron-Kerbosch with pivoting; avoids MATLAB's recursion-depth limit.
cliques={};stack={frame([],1:size(adj,1),[])};
while ~isempty(stack)
    cur=stack{end};
    if isempty(cur.P)&&isempty(cur.X)
        cliques{end+1}=cur.chosen;stack(end)=[];continue;
    end
    if isempty(cur.candidates),stack(end)=[];continue;end
    v=cur.candidates(1);cur.candidates(1)=[];
    child=frame([cur.chosen v],cur.P(adj(v,cur.P)),cur.X(adj(v,cur.X)));
    cur.P(cur.P==v)=[];cur.X=sort([cur.X v]);stack{end}=cur;stack{end+1}=child;
end
    function f=frame(chosen,P,X)
        candidates=[];
        if ~isempty(P)||~isempty(X)
            both=unique([P X]);[~,pick]=max(sum(adj(both,P),2));
            candidates=P(~adj(both(pick),P));
        end
        f=struct('chosen',chosen,'P',P,'X',X,'candidates',candidates);
    end
end
