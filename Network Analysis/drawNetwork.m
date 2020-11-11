close all

% This script plots the HeLa cell network on an RGB microcsopy image.
% It saves the image as <well>_network.tif, where <well> is the well name.
% You can find the resulting image in the well's folder.
  

%% --------------------------SPECIFY PARAMETERS----------------------------

root = '../Experiments/WKS024/20x';
well = 'B03';               % well name
nodeSize = 6;               % node size
nodeColor = 'w';            % color of nodes ('w'=white, 'k'=black, 'g'=green, etc)
lineWidth = 4;              % thickness of edges (= lines)
edgeColor = 'w';            % color of edges ('w'=white, 'k'=black, 'g'=green, etc)
edgeTransparency = 0.5;     % transparency of edges (0=fully transparent, 1=not transparent)

%% ------------------------------START CODE--------------------------------

% Load image and graph if this wasn't done already
well_folder = fullfile(root, well);

if ~exist('allData','var')
    allData = struct;
end
if ~isfield(allData, well)
    allData = load_data_from_well(allData, well, well_folder);
end

% Get image and graph of current well
fused = allData.(well).fused;
cmap = allData.(well).cmap;
G = allData.(well).G;
xNodes = allData.(well).xNodes;
yNodes = allData.(well).yNodes;

N = size(fused, 1);
M = size(fused, 2);

% Draw graph as network & save the file
figure()
imshow(fused,cmap)

hold on

plot(G, 'XData', xNodes, 'YData', yNodes,...
    'MarkerSize',nodeSize,...
    'NodeLabel',{},...
    'NodeColor',nodeColor,...
    'LineWidth',lineWidth,...
    'EdgeColor',edgeColor,...
    'EdgeAlpha',edgeTransparency);

saveas(gcf, fullfile(well_folder,[well,'_network.tif']))
disp(['Image saved in ', well_folder,' directory.'])

%% ------------------------------FUNCTIONS---------------------------------

function [G, xNodes, yNodes] = load_graph(well_folder)

    % This function reads the .csv files created by fiji and converts the 
    % information into a ML graph structure.
    
    % Load image
    tic

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

function allData = load_data_from_well(allData, well, well_folder)

    % This function loads the image and graph corresponding to a well
    % and stores them in the allData datastructure.
    
    disp(['Loading the image of well ',well,'...'])
    [fused, cmap] = imread(fullfile(well_folder,[well,'_fused_RGB.tif']));
    allData.(well).fused = fused;
    allData.(well).cmap = cmap;
    disp(['Loading the graph of well ',well,'...'])
    [G, xNodes, yNodes] = load_graph(well_folder);
    allData.(well).G = G;
    allData.(well).xNodes = xNodes;
    allData.(well).yNodes = yNodes;
end