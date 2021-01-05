function Gv = voronoi_graph(xNodes, yNodes)
    % Returns graph structure corresponding to the voronoi graph
    % with nodes at positions (xNodes, yNodes).

    % check input sizes
    if size(xNodes, 1) == 1
        xNodes = xNodes';
    elseif (size(xNodes, 1) ~= 1) && (size(xNodes, 2) ~= 1)
        exit('Invalid array xNodes.')
    end
    if size(yNodes, 1) == 1
        yNodes = yNodes';
    elseif (size(yNodes, 1) ~= 1) && (size(yNodes, 2) ~= 1)
        exit('Invalid array yNodes.')
    end 
    
    % Get vertices and indeces (v,c) of voronoi diagram
    % corresponding to xNodes, yNodes.
    P = [xNodes, yNodes];
    numNodes = size(P,1);
    [v,c] = voronoin(P);

    % store nodes connected to vertices in connectedNodes
    connectedNodes = cell(1, size(v,1));
    for n=1:length(c)
        connectedVertices = c{n};
        for vj=connectedVertices
            if ~ismember(n, connectedNodes{vj})
                connectedNodes{vj} = [connectedNodes{vj}, n];
            end
        end
    end
    % remove first entry (corresponding to the vertex at infinity)
    connectedNodes(1) = [];

    % Initialize graph
    Gv = graph();
    Gv = addnode(Gv, numNodes);

    % Fill graph: nodes that share a vertex are connected
    % Loop over verteces
    for i=1:length(connectedNodes)
        % Get array with nodes connected to this vertex
        nodeArray = connectedNodes{i};
        while any(nodeArray)
            % Pick first node in this array and
            % loop over the remaining nodes.
            n1 = nodeArray(1);
            for j = 2:length(nodeArray)
                if ~(findedge(Gv, n1, nodeArray(j)))
                    Gv = addedge(Gv, n1, nodeArray(j));
                end
            end
            % remove first node from the array.
            % continue unitil there are no more nodes left.
            nodeArray(1) = [];
        end
    end
end
