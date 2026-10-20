function results = runNetworkA_cleaned(cfg)
% Run one Case-A network realization with cleaned modular layout.
%
% Core model behavior is preserved from NetworkA_single_v2.m:
% - manuscript-consistent SSA event construction with user-approved choices
% - b*n^0.6 placement for newly created crosslinkers
% - tensile-only bonded-chain spring
% - harmonic pairwise repulsion (bonded pairs optional via cfg.repelBondedPairs)
% - periodic image tracking on bonds

if nargin < 1
    cfg = struct();
end

defaults = defaultNetworkAConfig();
cfg = mergeConfig(defaults, cfg);
cfg.Ntotal = resolveNtotal(cfg);

outputDir = '';
if isfield(cfg, 'outputDir') && ~isempty(cfg.outputDir)
    outputDir = cfg.outputDir;
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
end
dataPath = resolveOutputPath(cfg.dataFile, outputDir);
videoPathBase = stripMp4Suffix(resolveOutputPath(cfg.videoName, outputDir));

Abox = cfg.L^2;

NA0 = round(cfg.xA * cfg.Ntotal);
NB0 = cfg.Ntotal - NA0;

if mod(cfg.fA * NA0, 2) ~= 0
    warning(['fA*NA0 = %d is odd. Full saturation requires an even total number ' ...
             'of crosslinker arms because every completed chain consumes two arms.'], ...
             cfg.fA * NA0);
end

% Constant feed-based mean-field repulsion radius.
nAvgFeed = NB0 / max(NA0, 1);
rEq = 0.5 * cfg.b * sqrt(nAvgFeed);

nA = NA0 - 1;
nB = NB0;
t = 0;

nodes = [cfg.L/2, cfg.L/2, 0];                      % [x,y,degree]
activeEnds = repmat([cfg.L/2, cfg.L/2, 1, 0], cfg.fA, 1);   % [x,y,parent,n]
bonds = zeros(0,5);                                  % [i,j,n,imgX,imgY]
history = zeros(0,10);                               % [event,t,dt,type,nA,nB,nNodes,nBonds,nActive,eqIters]

histT = 0;
histA = nA;
histB = nB;

viz = [];
if cfg.showLivePlot
    viz = initVisualization(cfg, histT, histA, histB);
    updateLivePlot(viz, nodes, bonds, activeEnds, cfg, histT, histA, histB);
end

recording = cfg.makeVideo && cfg.showLivePlot;
if recording
    vid = VideoWriter(videoPathBase, 'MPEG-4');
    vid.FrameRate = cfg.videoFPS;
    open(vid);
    writeVideo(vid, getframe(viz.fig));
end

runtime = tic;

eventID = 0;
while true
    eventID = eventID + 1;

    if isempty(activeEnds)
        break;
    end

    nodes(:,3) = computeDegrees(bonds, size(nodes,1));
    expectedEnds = sum(cfg.fA - nodes(:,3));
    if size(activeEnds,1) ~= expectedEnds
        error('Active-end bookkeeping mismatch: have %d, expected %d.', size(activeEnds,1), expectedEnds);
    end

    M = size(activeEnds,1);
    R = zeros(0,4); % [type,end1,end2,weight]

    nBare = sum(activeEnds(:,4) == 0);
    nGrowing = M - nBare;
    CM = nB / Abox;
    CX = nA / Abox;

    for e = 1:M
        n = activeEnds(e,4);
        if nB > 0
            R(end+1,:) = [1, e, 0, cfg.kChem * CM];
        end
        if nA > 0 && n >= 1
            R(end+1,:) = [2, e, 0, cfg.fA * cfg.kChem * CX];
        end
    end

    if M >= 2
        adjNow = makeAdjacency(bonds, size(nodes,1));
        for e = 1:M-1
            p1 = activeEnds(e,3);
            n1 = activeEnds(e,4);
            for o = e+1:M
                p2 = activeEnds(o,3);
                n2 = activeEnds(o,4);

                if adjNow(p1,p2)
                    continue;
                end
                if p1 == p2 && ~cfg.allowSelfLoops
                    continue;
                end
                if n1 == 0 && n2 == 0
                    continue;
                end

                nTot = n1 + n2;
                if p1 == p2
                    d2 = 0;
                else
                    dv = nodes(p2,1:2) - nodes(p1,1:2);
                    dv = dv - cfg.L * round(dv / cfg.L);
                    d2 = sum(dv.^2);
                end

                wEndEnd = cfg.kChem / (pi * nTot * cfg.b^2) * exp(-d2 / (nTot * cfg.b^2));
                if wEndEnd > 1e-14
                    R(end+1,:) = [3, e, o, wEndEnd];
                end
            end
        end
    end

    if isempty(R)
        degNow = computeDegrees(bonds, size(nodes,1));
        fprintf('\nSSA JAMMED at t = %.6g\n', t);
        fprintf('Free crosslinkers = %d, free B2 = %d, active ends = %d\n', nA, nB, M);
        fprintf('Bare crosslinker arms (n=0) = %d, growing B-chain ends (n>0) = %d\n', nBare, nGrowing);
        fprintf('Degree counts:');
        for dd = 0:cfg.fA
            fprintf(' d%d=%d', dd, sum(degNow == dd));
        end
        fprintf('\n');
        break;
    end

    w = R(:,4);
    a0 = sum(w);
    dtRxn = -log(max(rand, realmin)) / a0;
    t = t + dtRxn;

    pick = find(cumsum(w) >= rand * sum(w), 1, 'first');
    type = R(pick,1);
    e = R(pick,2);
    o = R(pick,3);

    p = activeEnds(e,3);
    n = activeEnds(e,4);
    newBond = false;

    switch type
        case 1
            activeEnds(e,4) = n + 1;
            nB = nB - 1;

        case 2
            nChain = n;
            r0 = cfg.b * nChain^0.6;
            r = r0 * sqrt(-log(max(rand, realmin)));
            theta = 2 * pi * rand;
            newPos = nodes(p,1:2) + r * [cos(theta), sin(theta)];

            q = size(nodes,1) + 1;
            nodes(q,:) = [newPos, 1];
            bonds(end+1,:) = [p, q, nChain, 0, 0];
            nodes(p,3) = nodes(p,3) + 1;

            newPosWrap = wrapPoint(newPos, cfg.L);
            activeEnds(e,:) = [newPosWrap, q, 0];
            if cfg.fA > 2
                activeEnds(end+1:end+(cfg.fA-2),:) = repmat([newPosWrap, q, 0], cfg.fA-2, 1);
            end

            nA = nA - 1;
            newBond = true;

        case 3
            q = activeEnds(o,3);
            no = activeEnds(o,4);

            if p == q
                img = [0, 0];
            else
                img = -round((nodes(q,1:2) - nodes(p,1:2)) / cfg.L);
            end
            bonds(end+1,:) = [p, q, n + no, img];
            nodes(p,3) = nodes(p,3) + 1;
            nodes(q,3) = nodes(q,3) + 1;

            rm = sort([e, o], 'descend');
            activeEnds(rm(1),:) = [];
            activeEnds(rm(2),:) = [];
            newBond = true;
    end

    nodes(:,3) = computeDegrees(bonds, size(nodes,1));
    validateNetworkState(nodes, bonds, activeEnds, cfg.fA);

    if cfg.showLivePlot
        drawMotion(nodes, bonds, cfg.L, activeEnds, viz, cfg.fA, cfg.b);
        drawnow;
    end

    eqIters = 0;
    rmsSpeed = 0;
    if newBond
        [nodes, rmsSpeed, eqIters] = equilibrateNetwork( ...
            nodes, bonds, activeEnds, cfg, rEq, viz);

        if eqIters >= cfg.maxEq
            warning('Event %d: equilibration hit step cap (%d), RMS|v|=%.3e not yet below tolV.', ...
                eventID, cfg.maxEq, rmsSpeed);
        end

        if ~isempty(activeEnds)
            parent = activeEnds(:,3);
            activeEnds(:,1:2) = wrapPoint(nodes(parent,1:2), cfg.L);
        end
        validateNetworkState(nodes, bonds, activeEnds, cfg.fA);
    end

    history(end+1,:) = [eventID, t, dtRxn, type, nA, nB, size(nodes,1), size(bonds,1), size(activeEnds,1), eqIters];
    histT(end+1) = t;
    histA(end+1) = nA;
    histB(end+1) = nB;

    if cfg.showLivePlot && mod(eventID, cfg.plotEvery) == 0
        updateLivePlot(viz, nodes, bonds, activeEnds, cfg, histT, histA, histB);
        drawnow;
        if recording
            writeVideo(vid, getframe(viz.fig));
        end
    end
end

if ~isempty(bonds)
    [nodes, rmsSpeedF, itFinal] = equilibrateNetwork( ...
        nodes, bonds, activeEnds, localFinalCfg(cfg), rEq, viz);
else
    rmsSpeedF = 0;
    itFinal = 0;
end

if ~isempty(activeEnds)
    parent = activeEnds(:,3);
    activeEnds(:,1:2) = wrapPoint(nodes(parent,1:2), cfg.L);
end

nodes(:,3) = computeDegrees(bonds, size(nodes,1));
unsat = find(nodes(:,3) < cfg.fA);
elapsed = toc(runtime);

fprintf('\n=== FINAL STATE ===\n');
fprintf('SSA events: %d   SSA time: %.6g   runtime: %.2f s\n', size(history,1), t, elapsed);
fprintf('Junctions: %d   Chains: %d\n', size(nodes,1), size(bonds,1));
fprintf('Free crosslinkers: %d   Free B2: %d   Active ends: %d\n', nA, nB, size(activeEnds,1));
fprintf('Unsaturated junctions: %d\n', numel(unsat));
if ~isempty(activeEnds)
    fprintf('Remaining bare crosslinker arms: %d   remaining growing chain ends: %d\n', ...
        sum(activeEnds(:,4) == 0), sum(activeEnds(:,4) > 0));
end
fprintf('Final equilibration: %d integration steps, RMS node speed = %.3e\n', itFinal, rmsSpeedF);
if rmsSpeedF >= cfg.tolVfinal
    warning('Final mechanical relaxation did not reach tolVfinal: RMS|v| = %.3e (target %.3e).', ...
        rmsSpeedF, cfg.tolVfinal);
end

[Etot, Ebond, Erep] = networkEnergy(nodes, bonds, cfg.b, cfg.kT, cfg.kRep, rEq, cfg.L, ...
    cfg.repelBondedPairs, cfg.useRepulsionBinning,     cfg.repulsionFullPairThreshold, cfg.kRepAnh, cfg.chainModel);
    fprintf('Bonded-chain free energy (%s): %.4g kT\n', cfg.chainModel, Ebond);
fprintf('Pairwise repulsion energy: %.4g\n', Erep);
fprintf('Total model energy (bonded + repulsion): %.4g\n', Etot);

if isempty(unsat)
    fprintf('All crosslinker junctions reached degree %d through SSA events.\n', cfg.fA);
else
    fprintf('Unsaturated node indices: ');
    fprintf('%d ', unsat);
    fprintf('\n');
end

if cfg.saveData
    save(dataPath, 'nodes', 'bonds', 'activeEnds', 'history', 't', 'nA', 'nB', 'cfg');
end

if cfg.showLivePlot
    updateLivePlot(viz, nodes, bonds, activeEnds, cfg, histT, histA, histB);
    title(viz.ax1, sprintf('Final network: %d junctions, %d chains', size(nodes,1), size(bonds,1)));
    drawnow;
end

if recording
    for extra = 1:cfg.videoFPS*2
        writeVideo(vid, getframe(viz.fig));
    end
    close(vid);
    fprintf('Saved movie: %s.mp4\n', videoPathBase);
end

if ~isempty(bonds)
    plotChainDistributions(nodes, bonds, cfg, outputDir, cfg.showFinalHistogram);
end

results = struct();
results.nodes = nodes;
results.bonds = bonds;
results.activeEnds = activeEnds;
results.history = history;
results.t = t;
results.nA = nA;
results.nB = nB;
results.unsaturatedNodes = unsat;
results.finalRmsSpeed = rmsSpeedF;
results.finalEqIters = itFinal;
results.energy = struct('total', Etot, 'bonded', Ebond, 'repulsion', Erep);
results.rEq = rEq;
results.nAvgFeed = nAvgFeed;
results.config = cfg;
results.runtimeSeconds = elapsed;
results.outputDir = outputDir;
results.dataFile = dataPath;
if recording
    results.videoFile = [videoPathBase '.mp4'];
else
    results.videoFile = '';
end
end

function cfgOut = localFinalCfg(cfgIn)
% Reuse same settings but with final-stage tolerances/caps.
cfgOut = cfgIn;
cfgOut.tolV = cfgIn.tolVfinal;
cfgOut.maxEq = cfgIn.maxEqFinal;
end

function outPath = resolveOutputPath(inPath, outputDir)
if isempty(outputDir) || isAbsolutePath(inPath)
    outPath = inPath;
else
    outPath = fullfile(outputDir, inPath);
end
end

function tf = isAbsolutePath(p)
tf = false;
if isempty(p)
    return;
end
if numel(p) >= 2 && p(2) == ':'
    tf = true;   % Windows drive path, e.g. C:\...
    return;
end
if numel(p) >= 2 && strcmp(p(1:2), '\\')
    tf = true;   % Windows UNC path
    return;
end
if p(1) == '/' || p(1) == '\'
    tf = true;
end
end

function out = stripMp4Suffix(pathNoExt)
if numel(pathNoExt) >= 4 && strcmpi(pathNoExt(end-3:end), '.mp4')
    out = pathNoExt(1:end-4);
else
    out = pathNoExt;
end
end
