function updateSphereNodes(h, xy, baseColors, markerArea, L)
% A batched, shaded 2D bead mesh: no 3D geometry or lighting enters the model.
% Truecolor vertices keep lighting independent of the axes/heatmap colormap.
persistent unitXY triangles diffuse highlight cachedResolution
% Small network markers need few triangles; enlarged/exported bead samples
% get a smoother outline and highlight without slowing ordinary live frames.
sectors = max(24,min(96,ceil(sqrt(markerArea))));
rings = max(6,min(24,ceil(sqrt(markerArea)/4)));
if isempty(unitXY) || ~isequal(cachedResolution,[sectors rings])
    cachedResolution = [sectors rings];
    angle = (0:sectors-1)'*(2*pi/sectors);
    radius = sin((1:rings)*pi/(2*rings));
    unitXY = [0 0;reshape(cos(angle)*radius,[],1),reshape(sin(angle)*radius,[],1)];
    triangles = zeros(sectors*(2*rings-1),3); row = 0;
    for j = 1:sectors
        next = mod(j,sectors)+1;
        row = row+1; triangles(row,:) = [1,1+j,1+next];
    end
    for ring = 2:rings
        inner = 1+(ring-2)*sectors; outer = inner+sectors;
        for j = 1:sectors
            next = mod(j,sectors)+1;
            row = row+1; triangles(row,:) = [inner+j,outer+j,outer+next];
            row = row+1; triangles(row,:) = [inner+j,outer+next,inner+next];
        end
    end
    z = sqrt(max(0,1-sum(unitXY.^2,2)));
    normal = [unitXY,z];
    light = [-.65 .75 1]; light = light/norm(light);
    halfway = light+[0 0 1]; halfway = halfway/norm(halfway);
    diffuse = (.25+.75*max(0,normal*light')).*(.45+.55*z.^.4);
    highlight = .26*max(0,normal*halfway').^18;
end
n = size(xy,1); verticesPerNode = size(unitXY,1);
if n==0
    set(h,'Vertices',zeros(0,2),'Faces',zeros(0,3),'FaceVertexCData',zeros(0,3));
    return;
end
% Match the former scatter marker's diameter in points on a square network.
ax = ancestor(h,'axes'); pixels = getpixelposition(ax,true);
pointsToPixels = get(groot,'ScreenPixelsPerInch')/72;
beadRadius = .5*sqrt(markerArea)*pointsToPixels*L/max(1,min(pixels(3:4)));
vertices = repmat(unitXY*beadRadius,n,1)+repelem(xy,verticesPerNode,1);
% A small neutral reflectance makes even black beads visibly rounded.
material = .12+.88*baseColors;
colors = repelem(material,verticesPerNode,1).*repmat(diffuse,n,3)+repmat(highlight,n,3);
colors = min(1,max(0,colors));
faces = repmat(triangles,n,1)+repelem((0:n-1)'*verticesPerNode,size(triangles,1),3);
set(h,'Vertices',vertices,'Faces',faces,'FaceVertexCData',colors);
end
