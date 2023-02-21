A = [0, 1, 1, 0, 0; 1, 0, 1, 0, 0; 1, 1, 0, 1, 1; 0, 0, 1, 0, 1; 0, 0, 1, 1, 0];
findCliques(A)  
function [findClique] = findCliques(adjacency)
    n = size(adjacency, 2);
    findClique = [];
    current = [];
    nodes = 1:n;
    used = [];
    bron_kerbosch(current, nodes, used)

    function [] = bron_kerbosch(current, nodes, used)
        if isempty(nodes) && isempty(used)
            newfindclique = zeros(1, n);
            newfindclique(current) = 1;
            findClique = [findClique newfindclique.'];
        else
            for u = nodes
            nodes = setxor(nodes,u);
            used = [used u];
            newCurrent = [current u];
            Nu = find(adjacency(u,:));
            newNodes = intersect(nodes,Nu);
            newUsed = intersect(used,Nu);
            bron_kerbosch(newCurrent, newNodes, newUsed);
     
            end
        end
    end  
end
