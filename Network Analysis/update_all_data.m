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
    
    allData.(well).diameter = str2double( T.diameter(row) ) * scale;
    xc = str2double( T.xc(row) );
    yc = str2double( T.yc(row) );

    allData.(well).xc = xc * scale;
    allData.(well).yc = yc * scale;
    
    [G, xNodes, yNodes] = load_graph(well_folder);
    
    % Read cell measurements table
    [Measurements, measurementNames, cellValues] = read_measurements_table(well_folder, scale);
    % Remove the nodes that are not part of a cell
    for i = 1:numnodes(G)
         if ~any( ismember(cellValues, i) )
            G = rmnode(G, i);
            xNodes(i) = [];
            yNodes(i) = [];
         end
    end

    allData.(well).G = G;
    allData.(well).xNodes = xNodes * scale;
    allData.(well).yNodes = yNodes * scale;
    allData.(well).area = Measurements.area;
    allData.(well).circularity = Measurements.circularity;
    allData.(well).longness = Measurements.longness;
    allData.measurementNames = measurementNames;

end