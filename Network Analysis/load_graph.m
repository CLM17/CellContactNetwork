function [G, xNodes, yNodes] = load_graph(well_folder)

    %----------------------------------------------------------------------
    % This function reads the .csv files created by fiji and converts the 
    % information into a ML graph structure.
    
    % INPUT PARAMETERS
    %-----------------
    
    % well_folder: string
    % Path to well folder.
    
    % OUTPUT PARAMETERS
    %------------------
    
    % G: graph
    % Graph object of cellular contact matrix.
    
    % xNodes: Nx1 array, where N is the number of nodes in the graph.
    % X positions of all nodes in graph.
    
    % yNodes: Nx1 array, where N is the number of nodes in the graph.
    % Y positions of all nodes in graph.
    
    %----------------------------------------------------------------------
    
    % Read csv files with edges & node positions
    edges = csvread(fullfile(well_folder,'edges.csv'));
    nuclei_com = csvread(fullfile(well_folder,'nuclei_com.csv'));
    numNodes = size(nuclei_com,1);
    
    % Node positions
    xNodes = nuclei_com(:,2);
    yNodes = nuclei_com(:,3);

    % Initialize graph
    G = graph();
    G = addnode(G, numNodes);

    % Add all edges in a for-loop
    for i = 1:size(edges,1)
        node1 = int64(edges(i,1));
        node2 = int64(edges(i,2));
        G = addedge(G, node1, node2);
    end
end