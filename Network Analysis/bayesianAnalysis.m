
well = 'D02';               % well name
experiment = 'WKS024';
magnification = 'M20';
fieldSize = 1104;
cutoffDistance = 250; % cutoff in micron
network_specifier = '_ml';

%% ------------------------------START CODE--------------------------------

% Load image and graph if this wasn't done already
root = fullfile('..','..','Experiments',experiment,magnification);
%root = 'M:\tnw\bn\dm\Shared\Kasper\PhD\MinimalCellCultures\Experiments_and_Analysis\Experiments\WKS024\10x';
well_folder = fullfile(root, well);

% Initialise allData
if ~exist('allData','var')
    allData.(experiment).(magnification) = struct;
end
if ~isfield(allData, experiment)
    allData.(experiment) = struct;
else
    if ~isfield(allData.(experiment), magnification)
        allData.(experiment).(magnification) = struct;
    end
end

% Read csv 
if ~exist('T','var')
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
end

scale = calculate_scale(magnification, fieldSize);
if ~isfield(allData.(experiment).(magnification), well)
    disp(['Loading the graph of well ', well, '...'])
    allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier);
end

% Get data of current well
G = allData.(experiment).(magnification).(well).G;
xNodes = allData.(experiment).(magnification).(well).xNodes;
yNodes = allData.(experiment).(magnification).(well).yNodes;

% well locations
xc = allData.(experiment).(magnification).(well).xc;
yc = allData.(experiment).(magnification).(well).yc;
diameter = allData.(experiment).(magnification).(well).diameter;

% pdf
[distances, existingEdgeLength] = distances_distributions(G, xNodes, yNodes);
[pd_dist, pd_edgeLength, density] = find_probability_distributions(G, distances, existingEdgeLength);


% Specify distance axis you want to plot
dAxis = 0:cutoffDistance;

 % Find the values of the probability density functions of the distance and
% edge distributions at these small distances. This takes a few minutes.
pDist = pdf(pd_dist, dAxis);
pEdgeLength = pdf(pd_edgeLength, dAxis);

% Apply Bayes
pConnect = density * pEdgeLength ./ pDist;

%% Plot distributions
figure()
plot(dAxis, pConnect)
hold on
plot(dAxis, pEdgeLength)
plot(dAxis, pDist)
hold off
%set(gca,'YScale', 'log')
legend('P(d_{ij})', 'P(d_{ij}|E_{ij})', 'P(E_{ij}|d_{ij})')

save('bayesianAnalysis.mat', 'well', 'magnification', 'experiment', 'dAxis',...
     'pDist', 'pEdgeLength', 'pConnect','pd_dist','pd_edgeLength', 'density')
 
 %%
figure()
[hDistances, eDistances] = histcounts(distances);
[hExistingEdgeLength, eExistingEdgeLength] = histcounts(existingEdgeLength);
bar(bin_centers(eDistances), hDistances / max(hDistances), ...
    'barWidth', 1, 'LineStyle', 'None', 'FaceAlpha', 0.7)
hold on
bar(bin_centers(eExistingEdgeLength), hExistingEdgeLength / max(hExistingEdgeLength), ...
    'barWidth', 1, 'LineStyle', 'None', 'FaceAlpha', 0.7)
set(gca, 'XScale', 'log')

%% Plot Centralities
% 
% figure()
% for i = 1:nC
%     cName = centralityNames{i};
%     cObserved = Centralities.(well).('G').(cName);
%     cSimulated = Centralities.(well).('GSim').(cName);
%     
%     subplot(3,nC,i)
%     p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
%     p1.NodeCData = cObserved;
%     colormap jet
%     title(cName)
%     
%     subplot(3,nC,i + nC)
%     p2 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
%     p2.NodeCData = cSimulated;
%     colormap jet
%     
%     subplot(3,nC,i + 2*nC)
%     histogram(cObserved)
%     hold on
%     histogram(cSimulated)
%     hold off
% end
% 
% subplot(3,nC,1)
% ylabel('Observed')
% 
% subplot(3,nC,nC + 1)
% ylabel('Simulated')
% 
% subplot(3,nC,2*nC + 1)
% ylabel('Histogram')



%% ------------------------------FUNCTIONS---------------------------------


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
    %pd_dist = truncate(pd_dist, min(existingEdgeLength), inf);
    pd_dist = truncate(pd_dist, 0, inf);
    pd_edgeLength = fitdist(existingEdgeLength, 'kernel','Width',5);
    %pd_edgeLength = truncate(pd_edgeLength, min(distances(distances > 0)), inf);
    pd_edgeLength = truncate(pd_edgeLength, 0, inf);

    density = calculate_density(G);
    disp(['Density is ', num2str(density)]);
    
end

function centers = bin_centers(edges)
    centers = (edges(2:end) + edges(1:end-1)) / 2;
end