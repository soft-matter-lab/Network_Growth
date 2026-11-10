function tests = testNetworkColors
tests = functiontests(localfunctions);
end
function setupOnce(tc)
tc.TestData.oldPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'SRC'));
tc.TestData.root = tempname; mkdir(tc.TestData.root);
fprintf('Network color artifacts: %s\n',tc.TestData.root);
end
function teardownOnce(tc)
path(tc.TestData.oldPath);
end
function s = snapshot()
s = struct('L',10,'b',1,'fA',4,'growthMode','preformed','phase','final', ...
    'event',6,'ssaTime',1,'usagePercent',100,'nodeCount',5,'bondCount',6);
s.nodes = [9 1 1;5 8 2;8 6 3;2 6 2;1 1 4];
s.bonds = [1 5 10 1 0;2 5 30 0 0;3 5 60 0 0;4 5 10 0 0;2 3 30 0 0;3 4 60 0 0];
s.activeEnds = zeros(0,4);
% Input order deliberately differs from the ascending-length color order.
s.preformed = struct('chainLengths',[60 10 30]);
v = s.nodes(s.bonds(:,2),1:2)+s.bonds(:,4:5)*s.L-s.nodes(s.bonds(:,1),1:2);
s.completed = struct('kuhnCount',s.bonds(:,3),'endToEnd',vecnorm(v,2,2));
end
function testFixedDegreeScaleAndBlackToggle(tc)
s = snapshot(); s.fA = 5; % highest observed degree is 4: must not become black
fig = plotNetworkDump(s,struct('visible','off','showDanglingClouds',false));
cleanup = onCleanup(@()close(fig));
nodes = findall(fig,'Tag','networkNodes');
verifyEqual(tc,nodes.CData,repmat([.8;.6;.4;.6;.2],1,3),'AbsTol',1e-14);
ax = findobj(fig,'Tag','network'); verifyEqual(tc,ax.CLim,[1 5]);
verifyEqual(tc,ax.Colormap([1 end],:),[.8 .8 .8;0 0 0]);
fig2 = plotNetworkDump(s,struct('visible','off','showNodeConnectivity',false,'colorPreformedBonds',false));
cleanup2 = onCleanup(@()close(fig2));
nodes2 = findall(fig2,'Tag','networkNodes'); verifyEqual(tc,nodes2.CData,zeros(5,3));
links = findobj(fig2,'Tag','networkBonds'); verifyEqual(tc,numel(links),1);
verifyEqual(tc,links.Color,[0 0 0]);
cb = findobj(fig2,'Type','colorbar');
verifyFalse(tc,any(arrayfun(@(h)strcmp(h.Label.String,'junction degree'),cb)));
end
function testLengthGroupsPeriodicPiecesAndLiveConsistency(tc)
s = snapshot();
fig = plotNetworkDump(s,struct('visible','off','showDanglingClouds',false, ...
    'outputFile',fullfile(tc.TestData.root,'trimodal_colors.png')));
cleanup = onCleanup(@()close(fig));
expected = [0 0 0;1 0 0;128/255 0 32/255];
lengths = [10 30 60];
for k = 1:3
    line = findobj(fig,'Tag','networkBonds','DisplayName',sprintf('N = %d',lengths(k)));
    verifyEqual(tc,line.Color,expected(k,:));
    if k==1
        % The N=10 strand crossing x=10/0 retains one color on both pieces.
        verifyTrue(tc,all(ismember([0 10],line.XData)));
    end
end
cfg = defaultNetworkBConfig(); cfg.L = s.L; cfg.b = s.b; cfg.fA = s.fA;
cfg.chainLengths = [60 10 30]; cfg.heatGrid = 20;
viz = initVisualization(cfg,0,0,[0 0 0],'off');
cleanup2 = onCleanup(@()close(viz.fig));
updateLivePlot(viz,s.nodes,s.bonds,s.activeEnds,cfg,0,0,[0 0 0]);
staticNodes = findall(fig,'Tag','networkNodes');
verifyEqual(tc,viz.hNodes.CData,staticNodes.CData);
verifyEqual(tc,viz.hNodes.CData(end,:),[0 0 0]);
for k = 1:3, verifyEqual(tc,viz.hChain(k).Color,expected(k,:)); end
drawMotion(s.nodes,s.bonds,s.L,s.activeEnds,viz,s.fA,s.b);
verifyEqual(tc,viz.hNodes.CData,staticNodes.CData); % relaxation redraw keeps truecolor
% Early state contains only the longest species; it must stay burgundy.
early = s.bonds(s.bonds(:,3)==60,:);
updateLivePlot(viz,s.nodes,early,s.activeEnds,cfg,0,0,[0 0 0]);
verifyEmpty(tc,viz.hChain(1).XData); verifyEmpty(tc,viz.hChain(2).XData);
verifyNotEmpty(tc,viz.hChain(3).XData); verifyEqual(tc,viz.hChain(3).Color,expected(3,:));
end
function testBimodalSingleAndLegacyStyles(tc)
s = snapshot(); s.preformed.chainLengths = [30 10];
style = networkColorStyle(s);
verifyEqual(tc,style.chainLengths,[10 30]); verifyEqual(tc,style.bondColors,[0 0 0;1 0 0]);
s.preformed.chainLengths = 30;
verifyEqual(tc,networkColorStyle(s).bondColors,[0 0 0]);
% Legacy Case A must never acquire species colors from its observed lengths.
s = rmfield(s,{'preformed','growthMode'});
verifyEmpty(tc,networkColorStyle(s).chainLengths);
% Self-loop uses two connections; isolated seed remains light gray.
s.bonds = [s.bonds;1 1 10 0 0];
degree = computeDegrees(s.bonds,6);
verifyEqual(tc,networkNodeColors(degree,4,true), ...
    repmat([.8/3;1.6/3;.8/3;1.6/3;0;.8],1,3),'AbsTol',1e-14);
end
