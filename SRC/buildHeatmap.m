function H = buildHeatmap(activeEnds, bonds, nodes, L, b, Xg, Yg)
% Visualization-only Gaussian cloud field.

H = zeros(size(Xg));
fadeRadiusFactor = 1.25;
edgeLevel = 0.05;
fadeLambda = -log(edgeLevel) / (fadeRadiusFactor^2);

for e = 1:size(activeEnds,1)
    p = activeEnds(e,3);
    n = activeEnds(e,4);
    if n <= 0
        continue;
    end

    rRMS = b * sqrt(n);
    s2 = rRMS^2 / max(fadeLambda, 1e-12);
    nodeW = wrapPoint(nodes(p,1:2), L);
    x0 = nodeW(1);
    y0 = nodeW(2);

    dx = Xg - x0;
    dx = dx - L * round(dx / L);
    dy = Yg - y0;
    dy = dy - L * round(dy / L);
    r2 = dx.^2 + dy.^2;
    H = max(H, exp(-r2 / s2));
end

if ~isempty(bonds)
    loopRows = find(bonds(:,1) == bonds(:,2) & bonds(:,3) > 0);
    for kk = 1:numel(loopRows)
        s = loopRows(kk);
        p = bonds(s,1);
        n = bonds(s,3);

        rRMS = b * sqrt(n);
        s2 = rRMS^2 / max(fadeLambda, 1e-12);
        nodeW = wrapPoint(nodes(p,1:2), L);
        x0 = nodeW(1);
        y0 = nodeW(2);

        dx = Xg - x0;
        dx = dx - L * round(dx / L);
        dy = Yg - y0;
        dy = dy - L * round(dy / L);
        r2 = dx.^2 + dy.^2;
        H = max(H, exp(-r2 / s2));
    end
end
end

