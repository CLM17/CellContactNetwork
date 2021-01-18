close all

% This script plots the HeLa cell network on an RGB microcsopy image.
% It saves the image as <well>_network.tif, where <well> is the well name.
% You can find the resulting image in the well's folder.

%% --------------------------SPECIFY PARAMETERS----------------------------

experiment = 'WKS024';
magnification = 'M20';
well = 'D02';               % well name
network_specifier = '_ml';
nodeData = 'none';   % change to centrality name ('betweenness', 'closeness', etc) or cell measurement ('area', 'circularity', 'longness') or 'none'.
nodeSize = 1;               % node size
nodeColor = 'k';            % color of nodes ('w'=white, 'k'=black, 'g'=green, etc)
lineWidth = 0.5;              % thickness of edges (= lines)
edgeColor = 'k';            % color of edges ('w'=white, 'k'=black, 'g'=green, etc)
edgeTransparency = 0.5;     % transparency of edges (0=fully transparent, 1=not transparent)
qualityCheck = false;

%% ------------------------------START CODE--------------------------------

green = [108,186,116] / 255;
pink = [144, 108, 186] / 255;
yellow = [250,255,164] / 255;
close all

colorHela = green;
colorCos = pink;

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

scale = 1; % Do not convert pixels to um
if ~isfield(allData.(experiment).(magnification), well)
    disp(['Loading the graph of well ', well, '...'])
    allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier);
end

% Get image and graph of current well
measurementNames = allData.measurementNames;
G = allData.(experiment).(magnification).(well).G;
xc = allData.(experiment).(magnification).(well).xc;
yc = allData.(experiment).(magnification).(well).yc;
xNodes = allData.(experiment).(magnification).(well).xNodes - xc; % set center of well to x=0
yNodes = yc - allData.(experiment).(magnification).(well).yNodes; % set center of well to y=0n = numnodes(G);

% Calculate centrality if the user asked for it
centralityNames = {'degree', 'betweenness', 'closeness', 'pagerank', 'eigenvector'};
if (ismember(nodeData, centralityNames) && ...
        ~isfield(allData.(experiment).(magnification).(well), nodeData))
    data = centrality(G, nodeData);
    allData.(experiment).(magnification).(well).(nodeData) = data;
end

% Draw graph as network & save the file
figure()

% p = plot(G, 'XData', xNodes, 'YData', yNodes, ...
%     'NodeColor', nodeColor,...
%     'MarkerSize',nodeSize,...
%     'NodeLabel',{},...
%     'LineWidth',lineWidth,...
%     'EdgeColor',edgeColor,...
%     'EdgeAlpha',edgeTransparency);

plot(xNodes, yNodes, '.', 'Color', colorHela, 'MarkerSize', 5)

if ismember(nodeData, centralityNames) || ismember(nodeData, measurementNames)
    rgb = vals2colormap(allData.(experiment).(magnification).(well).(nodeData));
    p.NodeColor = rgb;
end

% bc = calculate_normalized_centrality(G, 'betweenness');
% colors = ones(n, 3);
% colors(bc > 0.1,:) = repmat([0,1,0], sum(bc>0.1), 1);
% p.NodeColor = colors;

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 6 6 ])
saveas(gcf, fullfile('Figures', 'FullNetworks',[experiment, '_', well,'_positions_only.pdf']))

%% Quality check

if qualityCheck
    nImgs = 5;
    imgSize = 2 * 512;
    % set seed
    rng(1);
    xPos = randi(N - imgSize, [1, nImgs]);
    yPos = randi(M - imgSize, [1, nImgs]);

    for i = 1:nImgs
        x = xPos(i);
        y = yPos(i);
        img = fused(y:y+imgSize, x:x+imgSize, :);

        nodesInImg = (xNodes >= x & xNodes <= x+imgSize) & (yNodes >= y & yNodes <= y+imgSize);
        subG = subgraph(G, nodesInImg);

        disp(['Img ', num2str(i),': N nodes = ', num2str(numnodes(subG)), ', N edges = ', num2str(numedges(subG)),'.'])

        figure()
        imshow(img)
        hold on

        p = plot(subG, 'XData', xNodes(nodesInImg)-x, 'YData', yNodes(nodesInImg)-y, ...
            'NodeColor', nodeColor,...
            'MarkerSize',nodeSize,...
            'NodeLabel',{},...
            'LineWidth',lineWidth,...
            'EdgeColor',edgeColor,...
            'EdgeAlpha',edgeTransparency);

        saveas(gcf, fullfile('Figures', 'Quality assessment',['img',num2str(i),'_',network_specifier,'.tif']))
    end
end

%%
bc = calculate_normalized_centrality(G, 'betweenness');


%% ------------------------------FUNCTIONS---------------------------------

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