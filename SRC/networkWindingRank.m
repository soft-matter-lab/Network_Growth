function value=networkWindingRank(nNodes,bonds)
% Rank of the integer cycle windings of an induced periodic graph.
adj=cell(nNodes,1);
for e=1:size(bonds,1)
    i=bonds(e,1);j=bonds(e,2);q=bonds(e,4:5);
    adj{i}(end+1,:)=[j q];adj{j}(end+1,:)=[i -q];
end
seen=false(nNodes,1);offset=zeros(nNodes,2);cycles=zeros(0,2);
for root=1:nNodes
    if seen(root)||isempty(adj{root}),continue;end
    seen(root)=true;stack=root;
    while ~isempty(stack)
        i=stack(end);stack(end)=[];
        for a=1:size(adj{i},1)
            j=adj{i}(a,1);trial=offset(i,:)+adj{i}(a,2:3);
            if ~seen(j)
                seen(j)=true;offset(j,:)=trial;stack(end+1)=j;
            else
                q=trial-offset(j,:);if any(q),cycles(end+1,:)=q;end
            end
        end
    end
end
value=rank(cycles);
end
