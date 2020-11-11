function simulate(well, experiment)
    root = '/Users/lukasvdh/Documents/BEP/Experiments/WKS024/10x';
    well = 'B02';               % well name
    experiment = 'WKS024';
    cutoffDistance = 250; % cutoff in micron

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
        allData = update_all_data(allData, well, well_folder, T);
    end

    % Get data of current well
    G = allData.(well).G;
    xNodes = allData.(well).xNodes;
    yNodes = allData.(well).yNodes;

    % well locations
    xc = allData.(well).xc;
    yc = allData.(well).yc;
    diameter = allData.(well).diameter;

    % pdf
    pd_dist = allData.(well).pd_dist;
    pd_edgeLength = allData.(well).pd_edgeLength;
    density = allData.(well).density;

    %% Generate random node positions with the same radial and angular 
    % distributions as the experimental result.

    % Experimental radial and angular distributions
    r = sqrt( (xNodes - xc).^2 + (yNodes - yc).^2 );
    theta = atan2((yc - yNodes), (xNodes - xc));

    % Fit a kernel and truncate
    pdf_r = fitdist(r,'kernel');
    pdf_r = truncate(pdf_r, 0, inf);
    pdf_theta = fitdist(theta,'kernel');
    pdf_theta = truncate(pdf_theta, -pi, pi);

    % Generate random simulations with the same distribution
    nNodes = numnodes(G);
    rSim = random(pdf_r, nNodes, 1);
    theta_sim = random(pdf_theta, nNodes, 1);

    xSim = rSim .* cos(theta_sim);
    ySim = rSim .* sin(theta_sim);

    % Calculate simulated distances
    [X,Y] = meshgrid(xSim,ySim);
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
    
end

%% Calulcate centralities
Centralities = struct();
centralityNames = {'degree', 'betweenness', 'closeness', 'pagerank', 'eigenvector'};
nC = length(centralityNames);

for i = 1:nC
    cName = centralityNames{i};
    Centralities.(well).('G').(cName) = centrality(G, cName);
    Centralities.(well).('GSim').(cName) = centrality(GSim, cName);
end

%% Plot distributions
figure()
plot(smallDistances, pDist, '.')
hold on
plot(smallDistances, pEdgelength, '.')
plot(smallDistances, pConnect, '.')
hold off
legend('P(d_{ij})', 'P(d_{ij}|E_{ij})', 'P(E_{ij}|d_{ij})')

%% Plot Centralities

figure()
for i = 1:nC
    cName = centralityNames{i};
    cObserved = Centralities.(well).('G').(cName);
    cSimulated = Centralities.(well).('GSim').(cName);
    
    subplot(3,nC,i)
    p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
    p1.NodeCData = cObserved;
    colormap jet
    title(cName)
    
    subplot(3,nC,i + nC)
    p2 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
    p2.NodeCData = cSimulated;
    colormap jet
    
    subplot(3,nC,i + 2*nC)
    histogram(cObserved)
    hold on
    histogram(cSimulated)
    hold off
end

subplot(3,nC,1)
ylabel('Observed')

subplot(3,nC,nC + 1)
ylabel('Simulated')

subplot(3,nC,2*nC + 1)
ylabel('Histogram')

% figure()
% % Empty graph
% subplot(2,6,1)
% p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
% title('Observed graph')
% colormap jet
% colorbar
% 
% % degree centrality (hub identifyer)
% subplot(2,6,2)
% p2 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
% p2.NodeCData = dc;
% colormap jet
% title('Degree centrality')
% 
% % closeness centrality (nodes that communicate quickly with others)
% subplot(2,6,3)
% p3 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
% p3.NodeCData = cc;
% colormap jet
% title('Closeness centrality')
% 
% % betweeness centrality (nodes that often lie on on paths between other nodes)
% subplot(2,6,4)
% p4 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
% p4.NodeCData = bc;
% colormap jet
% title('Betweenness centrality')
% 
% % pagerank centrality (average time spent at each node during a random walk)
% subplot(2,6,5)
% p5 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
% p5.NodeCData = prc;
% colormap jet
% title('Pagerank centrality')
% 
% % eigenvector centrality (nodes connected to important neighbors)
% subplot(2,6,6)
% p6 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
% p6.NodeCData = evc;
% colormap jet
% title('Eigenvector centrality')
% 
% % Empty graph simulated
% subplot(2,6,7)
% p1 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
% title('Simulated graph')
% colormap jet
% colorbar
% 
% % degree centrality (hub identifyer)
% subplot(2,6,8)
% p2 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
% p2.NodeCData = dc_sim;
% colormap jet
% title('Degree centrality')
% 
% % closeness centrality (nodes that communicate quickly with others)
% subplot(2,6,9)
% p3 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
% p3.NodeCData = cc_sim;
% colormap jet
% title('Closeness centrality')
% 
% % betweeness centrality (nodes that often lie on on paths between other nodes)
% subplot(2,6,10)
% p4 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
% p4.NodeCData = bc_sim;
% colormap jet
% title('Betweenness centrality')
% 
% % pagerank centrality (average time spent at each node during a random walk)
% subplot(2,6,11)
% p5 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
% p5.NodeCData = prc_sim;
% colormap jet
% title('Pagerank centrality')
% 
% % eigenvector centrality (nodes connected to important neighbors)
% subplot(2,6,12)
% p6 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
% p6.NodeCData = evc_sim;
% colormap jet
% title('Eigenvector centrality')

% figure()
% % Plot histograms
% subplot(2,3,1)
% histogram(r)
% xlabel('Distance to center')
% hold on
% 
% subplot(2,3,2)
% histogram(dc)
% xlabel('Degree centrality')
% hold on
% 
% subplot(2,3,3)
% histogram(cc)
% xlabel('Closeness centrality')
% hold on
% 
% subplot(2,3,4)
% histogram(bc)
% xlabel('Betweenness centrality')
% hold on
% 
% subplot(2,3,5)
% histogram(prc)
% xlabel('Pagerank centrality')
% hold on
% 
% subplot(2,3,6)
% histogram(evc)
% xlabel('Eigenvector centrality')
% hold on
% 
% % Plot histograms
% subplot(2,3,1)
% histogram(rSim)
% xlabel('Distance to center')
% 
% subplot(2,3,2)
% histogram(dc_sim)
% xlabel('Degree centrality')
% 
% subplot(2,3,3)
% histogram(cc_sim)
% xlabel('Closeness centrality')
% 
% subplot(2,3,4)
% histogram(bc_sim)
% xlabel('Betweenness centrality')
% 
% subplot(2,3,5)
% histogram(prc_sim)
% xlabel('Pagerank centrality')
% 
% subplot(2,3,6)
% histogram(evc_sim)
% xlabel('Eigenvector centrality')
% 
% %% Plot
% figure()
% plot(x,pdf_dist / sum(pdf_dist))
% hold on
% plot(x,pdf_edgeLength / sum(pdf_edgeLength))
% hold on
% plot(x,pdf_connect / sum(pdf_connect))
% legend('P(d_{ij})','P(d_{ij}|connect)','P(connect|d_{ij})')
% 
% figure()
% histogram(existingEdgeLength,'Normalization','probability','FaceColor','k')
% xlabel('Edge length (pixels)')
% ylabel('Probability')
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'Color','w','Units','inches','Position',[1 1 6 4])
% figName = fullfile('Figures',[experiment, '_', well,'_existingEdgeLengths.png']);
% saveas(gcf, figName)
% 
% figure()
% histogram(distances(distances>0),'Normalization','probability','FaceColor','k')
% xlabel('Edge length (pixels)')
% ylabel('Probability')
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'Color','w','Units','inches','Position',[1 1 6 4])
% figName = fullfile('Figures',[experiment, '_', well,'_allEdgeLengths.png']);
% saveas(gcf, figName)


%% ------------------------------FUNCTIONS---------------------------------

function allData = update_all_data(allData, well, well_folder, T)

    [G, xNodes, yNodes] = load_graph(well_folder);
    allData.(well).G = G;
    allData.(well).xNodes = xNodes;
    allData.(well).yNodes = yNodes;
    
    row = find( strcmp(T.well, well) );
    allData.(well).xc = T.xc(row);
    allData.(well).yc = T.yc(row);
    allData.(well).diameter = T.diameter(row);
    
    [distances, existingEdgeLength] = distances_distributions(G, xNodes, yNodes);
    [pd_dist, pd_edgeLength, density] = find_probability_distributions(G, distances, existingEdgeLength);
    allData.(well).pd_dist = pd_dist;
    allData.(well).pd_edgeLength = pd_edgeLength;
    allData.(well).density = density;
    
end

function [G, xNodes, yNodes] = load_graph(well_folder)

    % This function reads the .csv files created by fiji and converts the 
    % information into a ML graph structure.
    
    % Read csv files with edges & node positions
    tic
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