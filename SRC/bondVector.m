function dv = bondVector(nodes, bond, L)
p = bond(1);
q = bond(2);
img = bond(4:5);
dv = nodes(q,1:2) + img * L - nodes(p,1:2);
end

