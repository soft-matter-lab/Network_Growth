function rgb = danglingCloudColors(activeEnds,nodes,L,b,Xg,Yg)
% Visualization-only dangling-chain envelopes: white tails to dark-blue centers.
% Reuse the live cloud widths/periodic geometry; omit completed self-loop clouds.
% Truecolor keeps junction-degree and distribution colormaps independent.
H = buildHeatmap(activeEnds,zeros(0,5),nodes,L,b,Xg,Yg);
stops = [1 1 1; .80 .91 .98; .32 .60 .83; .025 .15 .38];
colors = interp1([0 .2 .6 1],stops,min(1,max(0,H(:))));
rgb = reshape(colors,[size(H),3]);
end
