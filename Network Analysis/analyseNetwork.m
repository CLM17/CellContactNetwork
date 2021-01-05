close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters

%   // experiment: string
%   Experiment name.

%   // well: string
%   Name of well (example: 'B02')

%   // fieldSize: double
%   Size of 1 field (is used to calculate the scale in um/pixel).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experiment = 'WKS024';
magnification = '20x';
well = 'B02';
fieldSize = 1104;
network_specifier = '_ml';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate scale in um/pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale = calculate_scale(magnification, fieldSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load xlsx file with well locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root = fullfile('../../Experiments', experiment, magnification);
well_folder = fullfile(root, well);

% Load image and graph if this wasn't done already
if ~exist('T','var')
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
end

if ~exist('allData','var')
    allData = struct;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis for 1 well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data for this well (if it wasn't done already)
if ~isfield(allData, well)
    allData = update_all_data(allData, well, well_folder, T, scale, network_specifier);
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

%% Centrality measures

% Calculate centralities
NetworkProperties = struct();
NetworkCentralities = struct();
propertyNames = {'rDensity', 'theta'};
centralityNames = {'degree', 'betweenness', 'closeness', 'pagerank', 'eigenvector'};

% Calculate the cell density as a function of radius
nBins = 20;
r = sqrt( xNodes.^2 + yNodes.^2 );
binEdges = linspace(0, max(r), nBins+1);
binArea = zeros(1, nBins);
rBinned = zeros(1, nBins);
for i = 1:nBins
    rBinned(i) = sum(r >= binEdges(i) & r<binEdges(i+1));
    binArea(i) = pi * ( binEdges(i+1)^2 - binEdges(i)^2 ); 
end

binCenters = (binEdges(2:end) + binEdges(1:end-1)) / 2;
NetworkProperties.(well).('rDensity') = rBinned ./ binArea;  
NetworkProperties.(well).('theta') = atan2(yNodes, xNodes);

% Calculate network centralities
nC = length(centralityNames);
for i = 1:nC
    cName = centralityNames{i};
    c = calculate_normalized_centrality(G, cName);
    NetworkCentralities.(well).(cName) = c;
end

% Plot biophysical properties (R, theta)
figure()
subplot(1,2,1)
bar(binCenters, NetworkProperties.(well).('rDensity'),'FaceColor','k');
xlabel('Radial distance to well center (um)')
ylabel('Density (# cells / um^2)')

subplot(1,2,2)
histogram(NetworkProperties.(well).('theta'), nBins, 'Normalization', 'probability','FaceColor','k');
xlabel('Angular position (radians)')
ylabel('Frequency')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 4])
figName = fullfile('Figures/AnalyzeWells/',[experiment, '_',magnification, '_', well,'_properties.png']);
saveas(gcf, figName)

% Plot centrality measures
figure()
for i = 1:nC

    cName = centralityNames{i};
    c = NetworkCentralities.(well).(cName);
    meanC = mean(c);
    medianC = median(c);

    subplot(2,nC,i)
    p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
    p1.NodeCData = c;
    colormap jet
    title(cName)

    subplot(2,nC,i+nC)
    histogram(c)
    xlabel(cName)
    hold on
    xline(meanC, 'r', 'LineWidth', 1.5);
    %xline(medianC, 'g', 'LineWidth', 2);
    if i==nC
        legend('Data', 'Mean')
    end

end

% save
set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 18 7.5])
figName = fullfile('Figures/AnalyzeWells/',[experiment, '_',magnification, '_', well,'_centralities.png']);
saveas(gcf, figName)
%disp('Centralities saved')

% Calculate distributions of distances
[X,Y] = meshgrid(xNodes, yNodes);
distances = sqrt( (X-X').^2 + (Y-Y').^2);
edgeLength = distances .* adjacency(G);
existingEdgeLength = edgeLength(edgeLength > 0);

% Plot distributions of distances
figure()
subplot(1,2,1)
histogram(existingEdgeLength,'Normalization','probability','FaceColor','k')
xlabel('Edge length (um)')
ylabel('Probability')

subplot(1,2,2)
histogram(distances(distances>0),'Normalization','probability','FaceColor','k')
xlabel('Distances (um)')
ylabel('Probability')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 4])
figName = fullfile('Figures/AnalyzeWells/',[experiment, '_',magnification, '_', well,'_distanceDistributions.png']);
saveas(gcf, figName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function density = calculate_density(G)
    nNodes = numnodes(G);
    nEdges = numedges(G);
    density = 2 * nEdges / (nNodes * (nNodes - 1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = shortest_path_lengths(G)
    d = distances(G);
    maxD = max(d(d~=inf));
    d(d==inf) = maxD + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
