function rdf = networkRdf(nodes, L, edges)
% 2D junction-junction RDF in a square periodic box, including bonded pairs.
% Unordered pair counts / [N*(N-1)/2 * annulus area / L^2].
% Restrict r <= L/2 so circular shells fit the minimum-image square.
validateattributes(L, {'numeric'}, {'scalar','real','finite','positive'});
if nargin < 3 || isempty(edges), edges = linspace(0,L/2,41); end
validateattributes(edges, {'numeric'}, {'vector','real','finite','nonnegative'});
edges = edges(:)';
if numel(edges) < 2 || any(diff(edges) <= 0) || edges(end) > L/2
    error('NetworkGrowth:RdfEdges','RDF edges must increase and lie in [0,L/2].');
end
counts = zeros(1,numel(edges)-1);
N = size(nodes,1);
for i = 1:N-1
    dv = nodes(i+1:end,1:2) - nodes(i,1:2);
    dv = dv - L*round(dv/L);
    counts = counts + histcounts(sqrt(sum(dv.^2,2)),edges);
end
expected = N*(N-1)/2 * pi*diff(edges.^2)/L^2;
g = nan(size(counts));
if N >= 2, g = counts ./ expected; end
rdf = struct('edges',edges,'r',(edges(1:end-1)+edges(2:end))/2, ...
    'counts',counts,'g',g,'expectedIdealCounts',expected,'nodeCount',N);
end
