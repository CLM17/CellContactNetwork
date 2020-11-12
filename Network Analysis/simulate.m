function simulate(experiment, magnification, well, fieldSize, keepPositionsString, description)
    
    root = fullfile('..','..','Experiments', experiment, magnification);
    cutoffDistance = 250; % cutoff in micron
    
    if strcmp(keepPositionsString, 'true')
        keepPositions = true;
    elseif strcmp(keepPositionsString, 'false')
        keepPositions = false;
    else
        error('Invalid choice for keepDistancesString')
    end
    
    fieldSize = str2double(fieldSize);
    
    if strcmp(magnification, '10x')
        scale = 873.96 / fieldSize; % 873.96 is field size for 10x in CX7
    elseif strcmp(magnification, '20x')
        scale = 441.41 / fieldSize; % 441.41 is field size for 10x in CX7
    else
        error('Invalid choice for magnification')
    end

    %% ------------------------------START CODE--------------------------------

    % Load image and graph if this wasn't done already
    if ~exist('T','var')
        well_folder = fullfile(root, well);
        xlsfileName = fullfile(root, 'Well locations.xlsx');
        T = readtable(xlsfileName);
    end

    if ~exist('allData','var')
        allData = struct;
    end

    if ~isfield(allData, well)
        allData = update_all_data(allData, well, well_folder, T, scale);
        allData = add_distributions_to_all_data(allData, well);
    end
    disp('All data loaded.')

    % well locations
    xc = allData.(well).xc;
    yc = allData.(well).yc;
    diameter = allData.(well).diameter;
    
    % Get data of current well
    G = allData.(well).G;
    xNodes = allData.(well).xNodes - xc; % set center of well to x=0
    yNodes = yc - allData.(well).yNodes; % set center of well to y=0

    % pdf
    pd_dist = allData.(well).pd_dist;
    pd_edgeLength = allData.(well).pd_edgeLength;
    density = allData.(well).density;

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
        pdf_r = fitdist(r,'kernel');
        pdf_r = truncate(pdf_r, 0, inf);
        pdf_theta = fitdist(theta,'kernel');
        pdf_theta = truncate(pdf_theta, -pi, pi);

        % Generate random simulations with the same distribution
        rSim = random(pdf_r, nNodes, 1);
        theta_sim = random(pdf_theta, nNodes, 1);

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
    edges = edgeProbabilities > tril(rand(nNodes));

    % Make symmetric and make diagonal zeros.
    edges = triu(edges',1) + tril(edges);
    edges = edges - diag( diag(edges) );

    GSim = graph(edges);
    
    fname = fullfile(well_folder, [well,'_simulatedGraph_',description,'.mat']);
    save(fname, 'G', 'xNodes', 'yNodes', 'GSim', 'xSim', 'ySim', 'pDist', 'pEdgelength', 'pConnect', 'smallDistances', 'description')
    disp('Output saved.')
end

%% ------------------------------FUNCTIONS---------------------------------

function allData = add_distributions_to_all_data(allData, well)

    G = allData.(well).G;
    xNodes = allData.(well).xNodes;
    yNodes = allData.(well).yNodes;
    [distances, existingEdgeLength] = distances_distributions(G, xNodes, yNodes);
    [pd_dist, pd_edgeLength, density] = find_probability_distributions(G, distances, existingEdgeLength);
    allData.(well).pd_dist = pd_dist;
    allData.(well).pd_edgeLength = pd_edgeLength;
    allData.(well).density = density;
    
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
    
    pd_dist = fitdist(distances(distances > 0),'kernel');
    pd_dist = truncate(pd_dist, 0, inf);
    pd_edgeLength = fitdist(existingEdgeLength, 'kernel');
    pd_edgeLength = truncate(pd_edgeLength, 0, inf);

    density = calculate_density(G);
    disp(['Density is ', num2str(density)]);
end

