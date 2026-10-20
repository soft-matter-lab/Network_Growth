function adj = makeAdjacency(bonds, N)
adj = false(N,N);
if isempty(bonds)
    return;
end
ind1 = sub2ind([N,N], bonds(:,1), bonds(:,2));
ind2 = sub2ind([N,N], bonds(:,2), bonds(:,1));
adj([ind1; ind2]) = true;
end

