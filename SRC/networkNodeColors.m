function rgb = networkNodeColors(degree, fA, enabled)
% Degree one (and isolated seeds) is light gray; degree fA is black.
rgb = zeros(numel(degree),3);
if enabled
    shade = .8 * (fA-min(fA,max(1,degree(:))))/(fA-1);
    rgb = repmat(shade,1,3);
end
end
