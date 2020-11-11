function allData = update_all_data(allData, well, well_folder, T, scale)

    %----------------------------------------------------------------------
    % This function updates the allData datastructure with the graph
    % corresponding to 'well'.
    % Node positions are converted from pixel values to measures in um.
    
    % INPUT PARAMETERS
    %-----------------
    
    % allData: struct
    % Fields are the wells. Each well field contains relevant data for that 
    % well, including the graph and node positions.
    
    % well: string
    % Well name, e.g. 'B02'.
    
    % well_folder: string
    % path to well folder.
    
    % T: table
    % table with well locations. Read from  'Well locations.xlsx'.
    
    % scale: double
    % Scale of one pixel (in um/pixel).
    
    % OUTPUT PARAMETERS
    %------------------
    
    % allData: struct
    % The function adds a field for the well. It adds the centers
    % coordinates of the well and its diameter (xc, yc, diameter). Also, it
    % adds the graph (G) and the node positions (xNodes, yNodes).
    %----------------------------------------------------------------------
    
    
    row = find( strcmp(T.well, well) );
    
    allData.(well).diameter = T.diameter(row) * scale;
    allData.(well).xc = T.xc(row) * scale;
    allData.(well).yc = T.yc(row) * scale;

    [G, xNodes, yNodes] = load_graph(well_folder);
    allData.(well).G = G;
    allData.(well).xNodes = xNodes * scale;
    allData.(well).yNodes = yNodes * scale;
    
end