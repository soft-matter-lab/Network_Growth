function frameIndex = exportNetworkDumpVideo(runOrDumpDir, videoFile, options)
% Replay stored network snapshots to an MP4, loading one selected frame at a time.
% Equal duration per saved state; no reaction or mechanical-state interpolation.
if nargin < 3, options = struct(); end
options = mergeConfig(struct('videoFrameEvery',1,'videoFPS',6, ...
    'videoHoldLastSeconds',2,'videoPreview',false,'heatGrid',150, ...
    'showDanglingClouds',true,'danglingCloudOpacity',.85),options);
validateattributes(options.videoFrameEvery,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(options.videoFPS,{'numeric'},{'scalar','real','positive','finite'});
validateattributes(options.videoHoldLastSeconds,{'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(options.heatGrid,{'numeric'},{'scalar','integer','>=',2,'finite'});
validateattributes(options.showDanglingClouds,{'logical','numeric'},{'scalar','binary'});
validateattributes(options.danglingCloudOpacity,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
[files,metadata,dumpDir] = listNetworkDumpFiles(runOrDumpDir);
if isempty(files), error('NetworkGrowth:NoSnapshots','There are no completed snapshot files.'); end
chosen = unique([1:options.videoFrameEvery:numel(files),numel(files)]);
cfg = metadata.config; cfg.heatGrid = options.heatGrid;
visible = 'off'; if options.videoPreview, visible = 'on'; end
% Empty plot([],[]) returns no Line handles, so seed both reservoir curves.
viz = initVisualization(cfg,0,metadata.initialFreeCrosslinkers,metadata.initialMonomers,visible);
figureCleanup = onCleanup(@()close(viz.fig));
% Offline replay has a white background and optional blue dangling clouds.
% The live simulation keeps its existing palette and rendering behavior.
viz.showDanglingClouds = options.showDanglingClouds;
set(viz.hHeat,'CData',ones(cfg.heatGrid,cfg.heatGrid,3), ...
    'AlphaData',options.danglingCloudOpacity,'AlphaDataMapping','none','Tag','danglingClouds');
set(viz.hChain,'Color',[.30 .30 .30],'LineWidth',1.4);
set(viz.ax1,'Color','w');
set(viz.fig,'PaperPositionMode','auto','InvertHardcopy','off');
writer = VideoWriter(videoFile,'MPEG-4'); writer.FrameRate = options.videoFPS;
open(writer);
videoCleanup = onCleanup(@()close(writer));
histT = zeros(1,numel(chosen)); histA = histT; histB = histT;
frameIndex = table();
for j = 1:numel(chosen)
    file = files(chosen(j)).name;
    data = load(fullfile(dumpDir,file),'snapshot'); s = data.snapshot;
    if ~isfield(s,'nodes')
        error('NetworkGrowth:NoNetworkVideo','Video requires dumpSaveNetwork = true: %s',file);
    end
    histT(j) = s.ssaTime; histA(j) = s.freeCrosslinkers; histB(j) = s.freeMonomers;
    updateLivePlot(viz,s.nodes,s.bonds,s.activeEnds,cfg,histT(1:j),histA(1:j),histB(1:j));
    title(viz.ax1,sprintf('Saved state | A usage %.3g%% | event %d',s.usagePercent,s.event));
    title(viz.ax2,sprintf('t = %.5g | %s | equal frame intervals',s.ssaTime,s.phase),'Interpreter','none');
    drawnow;
    % Render hidden figures without screen capture or opening a GUI window.
    rgb = print(viz.fig,'-RGBImage','-r100');
    if mod(size(rgb,1),2), rgb(end+1,:,:) = 255; end
    if mod(size(rgb,2),2), rgb(:,end+1,:) = 255; end
    writeVideo(writer,rgb);
    repeats = 0;
    if j == numel(chosen)
        repeats = round(options.videoHoldLastSeconds*options.videoFPS);
        for holdFrame = 1:repeats, writeVideo(writer,rgb); end
    end
    row = table(j,j+repeats,string(file),s.event,s.ssaTime,s.usagePercent,string(s.phase), ...
        'VariableNames',{'firstVideoFrame','lastVideoFrame','file','event','ssaTime','usagePercent','phase'});
    frameIndex = [frameIndex;row]; %#ok<AGROW>
    if j == 1 || mod(j,50) == 0 || j == numel(chosen)
        fprintf('Video: rendered %d/%d saved states.\n',j,numel(chosen));
    end
end
% Close the writer before returning so the MP4 is ready for playback.
clear videoCleanup figureCleanup;
end
