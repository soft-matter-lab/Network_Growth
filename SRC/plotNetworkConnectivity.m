function plotNetworkConnectivity(s,r,folder,options)
% Borderless high-resolution snapshots; all selected paths have equal width.
xy=mod(s.nodes(:,1:2),s.L);B=s.bonds(r.physicalBondRows,:);
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 4 4]);
set(fig,'PaperUnits','inches','PaperSize',[4 4],'PaperPosition',[0 0 4 4], ...
    'PaperPositionMode','manual','InvertHardcopy','off');
cleanup=onCleanup(@()close(fig));ax=axes(fig,'Position',[.02 .02 .96 .96]);
for mode=1:3
    cla(ax);hold(ax,'on');axis(ax,'equal');axis(ax,[0 s.L 0 s.L]);axis(ax,'off');
    [x,y]=buildSkeleton(s.nodes,B,s.L);plot(ax,x,y,'Color',[.85 .86 .87],'LineWidth',.45);
    if mode<3
        if mode==1,p=r.pathsX;name='paths_x';else,p=r.pathsY;name='paths_y';end
        [x,y]=buildSkeleton(s.nodes,B(p.highlightFlow~=0,:),s.L);
        plot(ax,x,y,'Color',[.72 .07 .13],'LineWidth',1);
        endpoints=s.bonds(:,1:2);degree=accumarray(endpoints(:),1,[size(xy,1) 1]);
        fA=max([2;degree]);if isfield(s.config,'fA'),fA=s.config.fA;end
        colors=networkNodeColors(degree,fA,options.showNodeConnectivity);
        scatter(ax,xy(:,1),xy(:,2),7,colors,'filled','MarkerEdgeColor','none');
        loopNodes=unique(s.bonds(r.selfLoopRows,1));
        if ~isempty(loopNodes)
            scatter(ax,xy(loopNodes,1),xy(loopNodes,2),20,'MarkerEdgeColor',[.9 .48 .1],'LineWidth',.6);
        end
    else
        name='rigid_clusters';rows=unique(vertcat(r.rigidity.clusterSourceBondRows{:}));
        [x,y]=buildSkeleton(s.nodes,s.bonds(rows,:),s.L);
        plot(ax,x,y,'Color',[.13 .40 .58],'LineWidth',.85);
        covered=r.rigidity.junctionClusterMultiplicity>0;
        scatter(ax,xy(~covered,1),xy(~covered,2),4,[.75 .77 .78],'filled');
        scatter(ax,xy(covered,1),xy(covered,2),5,[.13 .40 .58],'filled');
    end
    print(fig,fullfile(folder,[name '.png']),'-dpng',sprintf('-r%d',options.snapshotResolution));
end
end
