function [fig,result] = compareOriginalFloryFits(averageData,options)
% Compare an original, normalized Flory fit against the earlier conditional fit.
% Uses saved equal-replicate histograms; no simulation or histogram rescaling.
if nargin<2, options=struct(); end
options=mergeConfig(struct('minCount',4,'anchorCount',4,'anchorWeight',10, ...
    'bootstrapSamples',2000,'bootstrapSeed',1729,'kuhnYLim',[0 .5], ...
    'showSpread',true,'showConditionalFit',true,'visible','on','outputFile',''),options);
if ~isempty(options.kuhnYLim)
    validateattributes(options.kuhnYLim,{'numeric'},{'vector','numel',2,'real','finite','nonnegative','increasing'});
end
edges=averageData.edges{1}; counts=(edges(1:end-1)+edges(2:end))/2;
if any(abs(diff(edges)-1)>1e-9) || any(abs(counts-round(counts))>1e-9)
    error('NetworkGrowth:OriginalFloryBins', ...
        'This comparison needs individual integer-count bins. Re-run average_pkg_cleaned with nBinWidth=1.');
end
counts=round(counts); targets=averageData.targets; T=numel(targets);
result=struct('targets',targets,'fit',{{}},'fitTable',table(),'options',options, ...
    'sourceWeighting',averageData.weighting,'outputFile',options.outputFile);
curves=cell(1,T); conditional=cell(1,T); upper=[];
for t=1:T
    rows=averageData.perReplicate{1}(:,:,t); rows=rows(all(isfinite(rows),2),:);
    f=fitOriginalFloryNumberDistribution(rows,counts,options);
    result.fit{t}=f;
    row=table(targets(t),f.p,f.pCI95(1),f.pCI95(2),f.unanchoredP, ...
        f.anchorCount,f.anchorWeight,f.anchorObserved,f.anchorPredicted, ...
        f.anchorResidual,f.anchorMaximum,f.exactAnchorFeasible, ...
        f.tailSSE,f.unanchoredTailSSE,f.weightedSSE, ...
        f.observedTailMass,f.modelTailMass,f.replicateCount, ...
        'VariableNames',{'targetPercent','p','pCI95Low','pCI95High','unanchoredP', ...
        'anchorCount','anchorWeight','anchorObserved','anchorPredicted', ...
        'anchorResidual','anchorMaximum','exactAnchorFeasible', ...
        'tailSSE','unanchoredTailSSE','weightedSSE','observedTailMass','modelTailMass','replicateCount'});
    result.fitTable=[result.fitTable;row]; %#ok<AGROW>
    xx=unique([linspace(1,edges(end),1000),counts(counts>=1),f.anchorCount]);
    yy=(1-f.p)*f.p.^(xx-1);
    if f.p==0, xx=1; yy=1; end
    curves{t}=[xx;yy]; upper=[upper;yy(:)]; %#ok<AGROW>
    if options.showConditionalFit && isfield(averageData,'floryFit') && ...
            numel(averageData.floryFit)>=t && ~isempty(averageData.floryFit{t}) && ...
            (~isfield(averageData.floryFit{t},'modelStyle') || strcmp(averageData.floryFit{t}.modelStyle,'shifted'))
        previous=averageData.floryFit{t};
        [cx,cy]=floryDisplayCurve(previous,edges);
        keep=cx>=previous.minCount;
        conditional{t}=[cx(keep);cy(keep)];
        upper=[upper;cy(keep)']; %#ok<AGROW>
    end
end
meanValues=averageData.mean{1}; sdValues=averageData.sd{1};
spread=sdValues; spread(~isfinite(spread))=0;
if options.showSpread, upper=[upper;reshape(min(1,meanValues+spread),[],1)];
else, upper=[upper;meanValues(:)]; end
sharedYLim=options.kuhnYLim;
if isempty(sharedYLim), sharedYLim=[0 max(.05,1.1*max(upper))]; end
fig=figure('Color','w','Visible',options.visible,'Name','Original Flory comparison', ...
    'Position',[100 100 min(1800,180+420*T) 560]);
layout=tiledlayout(fig,1,T,'TileSpacing','compact','Padding','compact');
histX=reshape([edges(1:end-1);edges(2:end)],1,[]);
legendHandles=gobjects(0); legendAxes=[];
for t=1:T
    ax=nexttile(layout); hold(ax,'on'); f=result.fit{t};
    set(ax,'FontName','Arial','FontSize',11,'LineWidth',.8,'Box','off', ...
        'TickDir','out','XGrid','off','YGrid','on','GridAlpha',.1,'Layer','top');
    m=repelem(meanValues(t,:),2); sd=repelem(sdValues(t,:),2);
    hSpread=gobjects(0); hPrevious=gobjects(0);
    if options.showSpread && f.replicateCount>=2
        hSpread=fill(ax,[histX fliplr(histX)],[max(0,m-sd) fliplr(min(1,m+sd))],[0 .38 .68], ...
            'FaceAlpha',.1,'EdgeColor','none','DisplayName','+/- 1 replicate SD');
    end
    hMean=plot(ax,histX,m,'Color',[0 .38 .68],'LineWidth',2,'DisplayName','Replicate mean');
    if ~isempty(conditional{t})
        curve=conditional{t};
        hPrevious=plot(ax,curve(1,:),curve(2,:),'--','Color',[.6 .6 .6],'LineWidth',1.6, ...
            'DisplayName','Shifted Flory fit','Tag','conditionalComparisonCurve');
    end
    curve=curves{t};
    hOriginal=plot(ax,curve(1,:),curve(2,:),'k--','LineWidth',2.2, ...
        'DisplayName','Original Flory (soft anchor)','Tag','originalFloryCurve');
    if f.p==0, set(hOriginal,'Marker','o'); end
    plot(ax,f.anchorCount,f.anchorObserved,'o','MarkerSize',7,'LineWidth',1.5, ...
        'MarkerEdgeColor',[0 .38 .68],'MarkerFaceColor','w', ...
        'HandleVisibility','off','Tag','originalFloryAnchor');
    title(ax,sprintf('%g%% usage',targets(t)),'FontSize',14,'FontWeight','bold');
    subtitle(ax,{sprintf('Original p = %.4f',f.p), ...
        sprintf('n = %g: data %.3f | fit %.3f',f.anchorCount,f.anchorObserved,f.anchorPredicted)},'FontSize',11);
    xlim(ax,[edges(1) edges(end)]); ylim(ax,sharedYLim);
    handles=[hMean hSpread hOriginal hPrevious];
    if numel(handles)>numel(legendHandles), legendHandles=handles; legendAxes=ax; end
end
title(layout,sprintf('Original Flory distribution | n >= %g fit, anchor weight %g', ...
    options.minCount,options.anchorWeight),'FontName','Arial','FontSize',16,'FontWeight','bold');
xlabel(layout,'n (Kuhn segments)','FontName','Arial','FontSize',12);
ylabel(layout,'Chain number fraction','FontName','Arial','FontSize',12);
lgd=legend(legendAxes,legendHandles,'Orientation','horizontal','Box','off','FontSize',10);
lgd.Layout.Tile='south';
if ~isempty(options.outputFile), exportgraphics(fig,options.outputFile,'Resolution',180); end
end
