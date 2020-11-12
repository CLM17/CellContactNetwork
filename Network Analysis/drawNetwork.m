close all

% This script plots the HeLa cell network on an RGB microcsopy image.
% It saves the image as <well>_network.tif, where <well> is the well name.
% You can find the resulting image in the well's folder.

%% --------------------------SPECIFY PARAMETERS----------------------------

root = '../../Experiments/WKS024/10x';
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
if ~exist('T','var')
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
end

if ~isfield(allData, well)
    scale = 1;
    disp(['Loading the graph of well ', well, '...'])
    allData = update_all_data(allData, well, well_folder, T, scale);
    disp(['Loading the image for well ', well, '...'])
    [fused, cmap] = imread(fullfile(well_folder,[well,'_fused_RGB.tif']));
    allData.(well).fused = fused;
    allData.(well).cmap = cmap;
    G = allData.(well).G;
    n = numnodes(G);
    allData.(well).betweenness = centrality(G, 'betweenness');% / ((n-1)*(n-2));
    
end

% Get image and graph of current well
fused = allData.(well).fused;
cmap = allData.(well).cmap;
G = allData.(well).G;
xNodes = allData.(well).xNodes;
yNodes = allData.(well).yNodes;
n = numnodes(G);

N = size(fused, 1);
M = size(fused, 2);

% Draw graph as network & save the file
figure()
imshow(fused,cmap)

hold on

rgb = vals2colormap(allData.(well).betweenness);

p = plot(G, 'XData', xNodes, 'YData', yNodes, ...
    'NodeCData', allData.(well).betweenness, ...
    'NodeColor', rgb,...
    'MarkerSize',nodeSize,...
    'NodeLabel',{},...
    'LineWidth',lineWidth,...
    'EdgeColor',edgeColor,...
    'EdgeAlpha',edgeTransparency);

saveas(gcf, fullfile(well_folder,[well,'_network.tif']))
disp(['Image saved in ', well_folder,' directory.'])

%% ------------------------------FUNCTIONS---------------------------------

% function [G, xNodes, yNodes] = load_graph(well_folder)
% 
%     % This function reads the .csv files created by fiji and converts the 
%     % information into a ML graph structure.
%     
%     % Load image
%     tic
% 
%     % Read csv files with edges & node positions
%     tic
%     edges = csvread(fullfile(well_folder,'edges.csv'));
%     nuclei_com = csvread(fullfile(well_folder,'nuclei_com.csv'));
%     numNodes = size(nuclei_com,1);
%     
%     % Node positions
%     xNodes = nuclei_com(:,2);
%     yNodes = nuclei_com(:,3);
% 
%     % Initialize graph
%     G = graph();
%     G = addnode(G, numNodes);
% 
%     % Add all edges in a for-loop
%     for i = 1:size(edges,1)
%         node1 = int64(edges(i,1));
%         node2 = int64(edges(i,2));
%         G = addedge(G, node1, node2);
%     end
% end
% 
% function allData = load_data_from_well(allData, well, well_folder)
% 
%     % This function loads the image and graph corresponding to a well
%     % and stores them in the allData datastructure.
%     
%     disp(['Loading the image of well ',well,'...'])
%     [fused, cmap] = imread(fullfile(well_folder,[well,'_fused_RGB.tif']));
%     allData.(well).fused = fused;
%     allData.(well).cmap = cmap;
%     disp(['Loading the graph of well ',well,'...'])
%     [G, xNodes, yNodes] = load_graph(well_folder);
%     allData.(well).G = G;
%     allData.(well).xNodes = xNodes;
%     allData.(well).yNodes = yNodes;
% end

function rgb = vals2colormap(vals, colormap, crange)

% Take in a vector of N values and return and return a Nx3 matrix of RGB

% values associated with a given colormap

%

% rgb = AFQ_vals2colormap(vals, [colormap = 'jet'], [crange])

%

% Inputs:

% vals     = A vector of values to map to a colormap or a cell array of

%            vectors of values

% colormap = A matlab colormap. Examples: colormap = 'autumn';

%            colormap = 'jet'; colormap = 'hot';

% crange   = The values to map to the minimum and maximum of the colormap.

%            Defualts to the full range of values in vals.

%

% Outputs:

% rgb      = Nx3 matrix of rgb values mapping each value in vals to the

%            corresponding rgb colors.  If vals is a cell array then rgb

%            will be a cell array of the same length

%

% Example:

% vals = rand(1,100);

% rgb = AFQ_vals2colormap(vals, 'hot');

%

% Copyright Jason D. Yeatman, June 2012



if ~exist('colormap','var') || isempty(colormap)

    colormap = 'jet';

end



%

if ~iscell(vals)

    if ~exist('crange','var') || isempty(crange)

        crange = [min(vals) max(vals)];

    end

    % Generate the colormap

    cmap = eval([colormap '(256)']);

    % Normalize the values to be between 1 and 256

    vals(vals < crange(1)) = crange(1);

    vals(vals > crange(2)) = crange(2);

    valsN = round(((vals - crange(1)) ./ diff(crange)) .* 255)+1;

    % Convert any nans to ones

    valsN(isnan(valsN)) = 1;

    % Convert the normalized values to the RGB values of the colormap

    rgb = cmap(valsN, :);

elseif iscell(vals)

    if ~exist('crange','var') || isempty(crange)

        crange = [min(vertcat(vals{:})) max(vertcat(vals{:}))];

    end

    % Generate the colormap

    cmap = eval([colormap '(256)']);

    for ii = 1:length(vals)

        % Normalize the values to be between 1 and 256 for cell ii

        valsN = vals{ii};

        valsN(valsN < crange(1)) = crange(1);

        valsN(valsN > crange(2)) = crange(2);

        valsN = round(((valsN - crange(1)) ./ diff(crange)) .* 255)+1;

        % Convert any nans to ones

        valsN(isnan(valsN)) = 1;

        % Convert the normalized values to the RGB values of the colormap

        rgb{ii} = cmap(valsN, :);

    end

end

return
end