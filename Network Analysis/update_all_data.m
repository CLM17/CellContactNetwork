function allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier)

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
    
    if ~isfloat(T.diameter(row))
        diameter = str2double( T.diameter(row) );
    else
        diameter = T.diameter(row);
    end 
    if ~isfloat(T.xc(row))
        xc = str2double( T.xc(row) );
    else
        xc = T.xc(row);
    end
    if ~isfloat(T.yc(row))
        yc = str2double( T.yc(row) );
    else
        yc = T.yc(row);
    end 
    
    allData.(experiment).(magnification).(well).diameter = diameter * scale;
    allData.(experiment).(magnification).(well).xc = xc * scale;
    allData.(experiment).(magnification).(well).yc = yc * scale;
    
    [G, xNodes, yNodes] = load_graph(well_folder, network_specifier);
    
    % Read cell measurements table
    [Measurements, measurementNames, cellValues] = read_measurements_table(well_folder, scale, network_specifier);
    area = Measurements.area;
    circularity = Measurements.circularity;
    longness = Measurements.longness;
    % Remove the nodes that are not part of a cell
    c = 0;
    removeNodes = [];
    for i = 1:numnodes(G)
         if ~any( ismember(cellValues, i) )
             c = c + 1;
             removeNodes(c) = i;
         end
    end
    G = rmnode(G, removeNodes);
    xNodes(removeNodes) = [];
    yNodes(removeNodes) = [];
    
    if scale ~= 1
        cutOffArea = 8e4;
    else
        cutOffArea = Inf;
    end
    
    % remove node with degree 16
    highDegree = find(degree(G) == 16);
    highArea = find(area > cutOffArea);
    excluded = union(highDegree, highArea);
    
    G = rmnode(G, excluded);
    xNodes(excluded) = [];
    yNodes(excluded) = [];
    area(excluded) = [];
    circularity(excluded) = [];
    longness(excluded) = [];
    
    allData.(experiment).(magnification).(well).G = G;
    allData.(experiment).(magnification).(well).xNodes = xNodes * scale;
    allData.(experiment).(magnification).(well).yNodes = yNodes * scale;
    allData.(experiment).(magnification).(well).area = area;
    allData.(experiment).(magnification).(well).circularity = circularity;
    allData.(experiment).(magnification).(well).longness = longness;
    allData.measurementNames = measurementNames;

end