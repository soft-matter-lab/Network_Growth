function deg = computeDegrees(bonds, N)
deg = zeros(N,1);
if isempty(bonds)
    return;
end
deg = accumarray([bonds(:,1); bonds(:,2)], 1, [N,1]);
end

