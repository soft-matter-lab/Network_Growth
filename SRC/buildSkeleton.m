function [skX, skY] = buildSkeleton(nodes, bonds, L)
% Build drawable line segments for bonds with periodic boundary splitting.

skX = [];
skY = [];
for s = 1:size(bonds,1)
    i = bonds(s,1);
    j = bonds(s,2);
    if i == j
        continue;
    end

    p0 = wrapPoint(nodes(i,1:2), L);
    dv = bondVector(nodes, bonds(s,:), L);
    startPt = p0;
    rem = dv;
    piece = 0;

    while norm(rem) > 1e-12 && piece < 20
        piece = piece + 1;
        endPt = startPt + rem;
        if all(endPt >= -1e-12) && all(endPt <= L + 1e-12)
            endPt = min(max(endPt,0),L);
            skX = [skX, startPt(1), endPt(1), NaN];
            skY = [skY, startPt(2), endPt(2), NaN];
            break;
        end

        tcand = [];
        if rem(1) > 0 && endPt(1) > L
            tcand(end+1) = (L - startPt(1)) / rem(1);
        elseif rem(1) < 0 && endPt(1) < 0
            tcand(end+1) = (0 - startPt(1)) / rem(1);
        end
        if rem(2) > 0 && endPt(2) > L
            tcand(end+1) = (L - startPt(2)) / rem(2);
        elseif rem(2) < 0 && endPt(2) < 0
            tcand(end+1) = (0 - startPt(2)) / rem(2);
        end

        tcand = tcand(tcand > 1e-12 & tcand <= 1 + 1e-12);
        if isempty(tcand)
            skX = [skX, startPt(1), endPt(1), NaN];
            skY = [skY, startPt(2), endPt(2), NaN];
            break;
        end

        tcross = min(tcand);
        crossPt = startPt + tcross * rem;
        crossPt = min(max(crossPt,0),L);
        skX = [skX, startPt(1), crossPt(1), NaN];
        skY = [skY, startPt(2), crossPt(2), NaN];

        rem = (1 - tcross) * rem;
        startPt = crossPt;

        if abs(startPt(1)) < 1e-10 && rem(1) < 0
            startPt(1) = L;
        elseif abs(startPt(1)-L) < 1e-10 && rem(1) > 0
            startPt(1) = 0;
        end
        if abs(startPt(2)) < 1e-10 && rem(2) < 0
            startPt(2) = L;
        elseif abs(startPt(2)-L) < 1e-10 && rem(2) > 0
            startPt(2) = 0;
        end
    end
end
end

