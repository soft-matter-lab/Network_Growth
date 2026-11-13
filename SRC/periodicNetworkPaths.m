function out=periodicNetworkPaths(nodes,bonds)
% Maximum unit-capacity periodic circulation along x and y.
% Vertices obey flow conservation; every physical strand has |flow|<=1.
% The objective is net integer winding per cell, NOT a tensile force.
V=size(nodes,1);E=size(bonds,1);u=bonds(:,1);v=bonds(:,2);
if E==0
    out=cell(1,2);
    for axis=1:2
        out{axis}=struct('axis',axis,'capacity',0,'flow',zeros(0,1), ...
            'highlightFlow',zeros(0,1),'cycles',{{}},'highlightedStrands',0, ...
            'highlightedKuhn',0,'positiveWindingCycles',0,'dualPotential',zeros(V,1), ...
            'dualResidual',zeros(0,1),'dualObjective',0,'solverIterations',0, ...
            'primalResidual',0,'primalDualGap',0);
    end
    return;
end
if exist('linprog','file')~=2 || ~license('test','Optimization_Toolbox')
    error('NetworkGrowth:ConnectivityToolbox','Periodic path counts require MATLAB Optimization Toolbox (linprog).');
end
D=sparse([(1:E)';(1:E)'],[u;v],[-ones(E,1);ones(E,1)],E,V);
opts=optimoptions('linprog','Display','none','Algorithm','dual-simplex-highs');
out=cell(1,2);
for axis=1:2
    q=bonds(:,3+axis);
    [j,fun,flag,solver]=linprog(-q,[],[],D',zeros(V,1),-ones(E,1),ones(E,1),opts);
    assert(flag>0 && max(abs(j-round(j)))<1e-7);j=round(j);capacity=round(-fun);
    assert(max(abs(D'*j))<1e-8 && max(abs(j))<=1 && abs(q'*j-capacity)<1e-7);
    % Among maximum-capacity integer circulations choose minimum total contour.
    % The directed incidence constraints remain totally unimodular, so the
    % LP has an integral vertex optimum. A dominant reward fixes capacity first.
    cost=bonds(:,3);reward=2*sum(cost)+1;
    [pm,~,costFlag]=linprog([cost-reward*q;cost+reward*q],[],[],[D',-D'],zeros(V,1), ...
        zeros(2*E,1),ones(2*E,1),opts);
    assert(costFlag>0 && max(abs(pm-round(pm)))<1e-6);pm=round(pm);
    j=pm(1:E)-pm(E+1:end);
    assert(max(abs(D'*j))<1e-8 && abs(q'*j-capacity)<1e-7 && ~any(pm(1:E)&pm(E+1:end)));
    % Independent dual certificate: min sum |q+D*phi|, one gauge pinned.
    Dr=D(:,2:end);Aineq=[Dr,-speye(E);-Dr,-speye(E)];bineq=[-q;q];
    [dual,dval,dflag]=linprog([zeros(V-1,1);ones(E,1)],Aineq,bineq,[],[], ...
        [-inf(V-1,1);zeros(E,1)],[],opts);
    assert(dflag>0 && abs(dval-capacity)<1e-7);
    potential=[0;dual(1:V-1)];residual=q+D*potential;
    assert(abs(sum(abs(residual))-capacity)<1e-7);
    % Separate a conservative integer flow into edge-disjoint directed cycles.
    remaining=j;cycles={};highlight=zeros(E,1);total=0;
    while any(remaining)
        seed=find(remaining,1);if remaining(seed)>0,start=u(seed);else,start=v(seed);end
        visited=zeros(V,1);walkVertices=zeros(0,1);walkEdges=zeros(0,1);walkSigns=zeros(0,1);cur=start;
        while visited(cur)==0
            visited(cur)=numel(walkVertices)+1;walkVertices(end+1,1)=cur;
            e=find((remaining>0 & u==cur)|(remaining<0 & v==cur),1);assert(~isempty(e));
            direction=sign(remaining(e));walkEdges(end+1,1)=e;walkSigns(end+1,1)=direction;
            if direction>0,cur=v(e);else,cur=u(e);end
        end
        first=visited(cur);ee=walkEdges(first:end);ss=walkSigns(first:end);
        winding=sum(bonds(ee,4:5).*ss,1);remaining(ee)=remaining(ee)-ss;
        assert(winding(axis)>=0,'Negative objective cycle in an optimum');
        if winding(axis)>0
            cycles{end+1}=struct('edgeIndices',ee,'edgeDirections',ss,'winding',winding); %#ok<AGROW>
            highlight(ee)=ss;total=total+winding(axis);
        end
    end
    assert(total==capacity && max(abs(D'*highlight))<1e-8);
    out{axis}=struct('axis',axis,'capacity',capacity,'flow',j,'highlightFlow',highlight,'cycles',{cycles}, ...
        'highlightedStrands',nnz(highlight),'highlightedKuhn',sum(cost.*abs(highlight)), ...
        'positiveWindingCycles',numel(cycles),'dualPotential',potential, ...
        'dualResidual',residual,'dualObjective',dval,'solverIterations',solver.iterations, ...
        'primalResidual',max(abs(D'*j)),'primalDualGap',abs(dval-capacity));
end
end
