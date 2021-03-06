function simulate(x, experiment, magnification, well, network_specifier,...
                  fieldSizeString, keepPositionsString, description)
    
    root = fullfile('..','..','Experiments', experiment, magnification);
    cutoffDistance = 250; % cutoff in micron
    
    if strcmp(keepPositionsString, 'true')
        keepPositions = true;
    elseif strcmp(keepPositionsString, 'false')
        keepPositions = false;
    else
        error('Invalid choice for keepDistancesString')
    end
    
    fieldSize = str2double(fieldSizeString);
    nr = str2double(x);

    %% ----------------------------START CODE------------------------------

    % generate random stream
    myStream = RandStream('mlfg6331_64', 'Seed', nr);
    
    % Load image and graph if this wasn't done already
    
    well_folder = fullfile(root, well);
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);

    allData = struct;
    scale = calculate_scale(magnification, fieldSize);
    allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier);
    allData = add_distributions_to_all_data(allData, experiment, magnification, well);
    
    disp('All data loaded.')
    
    % well locations
    xc = allData.(experiment).(magnification).(well).xc;
    yc = allData.(experiment).(magnification).(well).yc;

    % Get data of current well
    G = allData.(experiment).(magnification).(well).G;
    xNodes = allData.(experiment).(magnification).(well).xNodes - xc; % set center of well to x=0
    yNodes = yc - allData.(experiment).(magnification).(well).yNodes; % set center of well to y=0

    % pdf
    pd_dist = allData.(experiment).(magnification).(well).pd_dist;
    pd_edgeLength = allData.(experiment).(magnification).(well).pd_edgeLength;
    density = allData.(experiment).(magnification).(well).density;

    %% Generate random node positions with the same radial and angular 
    % distributions as the experimental result.
    
    nNodes = numnodes(G);
    
    if keepPositions
        
        xSim = xNodes;
        ySim = yNodes;
        
    else
        
        % Experimental radial and angular distributions
        r = sqrt( xNodes.^2 + yNodes.^2 );
        theta = atan2(yNodes, xNodes);

        % Fit a kernel and truncate
        pdf_r = fitdist(r,'kernel','Width',5);
        pdf_r = truncate(pdf_r, 0, max(r));
        pdf_theta = fitdist(theta,'kernel','Width',5);
        pdf_theta = truncate(pdf_theta, -pi, pi);

        % Generate random simulations with the same distribution
        rSim = random(pdf_r, myStream, nNodes, 1);
        theta_sim = random(pdf_theta, myStream, nNodes, 1);

        xSim = rSim .* cos(theta_sim);
        ySim = rSim .* sin(theta_sim);
        
    end

    % Calculate simulated distances
    [X,Y] = meshgrid(xSim, ySim);
    distances = sqrt( (X-X').^2 + (Y-Y').^2);

    % Extract the distances smaller than the cutoff distance, with the
    % corresponding indeces.
    smallDistanceBoolean = (distances > 0 & distances < cutoffDistance);
    smallDistances = distances(smallDistanceBoolean);
    [iArray,jArray,~] = find(sparse(smallDistanceBoolean));

    % Find the values of the probability density functions of the distance and
    % edge distributions at these small distances. This takes a few minutes.
    pDist = pdf(pd_dist, smallDistances);
    pEdgelength = pdf(pd_edgeLength, smallDistances);

    % Apply Bayes
    pConnect = density * pEdgelength ./ pDist;

    % Fill an array with probabilities of making an edge.
    % It is 0 for all distances larger than the cutoff distance.
    edgeProbabilities = zeros(nNodes, nNodes);
    for n = 1:length(smallDistances)
        i = iArray(n);
        j = jArray(n);
        edgeProbabilities(i,j) = pConnect(n);
    end

    % An edge is made if the edge probability is larger than some random
    % number.
    edges = edgeProbabilities > tril(rand(myStream, nNodes));

    % Make symmetric and make diagonal zeros.
    edges = triu(edges',1) + tril(edges);
    edges = edges - diag( diag(edges) );

    GSim = graph(edges);
    
    fname = fullfile(well_folder, [well,'_simulatedGraph_',description, '_', x, '.mat']);
    save(fname, 'G', 'xNodes', 'yNodes', 'GSim', 'xSim', 'ySim', 'pDist', 'pEdgelength', 'pConnect', 'smallDistances', 'description')
    disp('Output saved.')
end

%% ---------------------------HELPER FUNCTIONS-----------------------------

function allData = add_distributions_to_all_data(allData, experiment, magnification, well)

    G = allData.(experiment).(magnification).(well).G;
    xNodes = allData.(experiment).(magnification).(well).xNodes;
    yNodes = allData.(experiment).(magnification).(well).yNodes;
    [distances, existingEdgeLength] = distances_distributions(G, xNodes, yNodes);
    [pd_dist, pd_edgeLength, density] = find_probability_distributions(G, distances, existingEdgeLength);
    allData.(experiment).(magnification).(well).pd_dist = pd_dist;
    allData.(experiment).(magnification).(well).pd_edgeLength = pd_edgeLength;
    allData.(experiment).(magnification).(well).density = density;
    
end

function density = calculate_density(G)
    nNodes = numnodes(G);
    nEdges = numedges(G);
    density = 2 * nEdges / (nNodes * (nNodes - 1));
end

function [distances, existingEdgeLength] = distances_distributions(G, xNodes, yNodes)
    % distance matrix
    [X,Y] = meshgrid(xNodes,yNodes);
    full_matrix = full( adjacency(G) );
    distances = sqrt( (X-X').^2 + (Y-Y').^2);
    edgeLength = distances .* full_matrix;
    existingEdgeLength = edgeLength(edgeLength > 0);
end

function [pd_dist, pd_edgeLength, density] = find_probability_distributions(G, distances, existingEdgeLength)
    
    pd_dist = fitdist(distances(distances > 0),'kernel','Width',5);
    pd_dist = truncate(pd_dist, 0, inf);
    pd_edgeLength = fitdist(existingEdgeLength, 'kernel','Width',5);
    pd_edgeLength = truncate(pd_edgeLength, 0, inf);

    density = calculate_density(G);
    disp(['Density is ', num2str(density)]);
end

