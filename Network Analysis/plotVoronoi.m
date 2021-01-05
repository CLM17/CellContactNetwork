close all

%% --------------------------SPECIFY PARAMETERS----------------------------

experiment = 'JJ005';
magnification = '20x';
well = 'C02';               % well name
network_specifier = '';
nodeSize = 7;               % node size
nodeColor = 'w';            % color of nodes ('w'=white, 'k'=black, 'g'=green, etc)
lineWidth = 3;              % thickness of edges (= lines)
edgeColor = 'w';            % color of edges ('w'=white, 'k'=black, 'g'=green, etc)
edgeTransparency = 0.5;     % transparency of edges (0=fully transparent, 1=not transparent)
qualityCheck = false;

%% ------------------------------START CODE--------------------------------

% Load image and graph if this wasn't done already
root = fullfile('..','..','Experiments',experiment,magnification);
%root = 'M:\tnw\bn\dm\Shared\Kasper\PhD\MinimalCellCultures\Experiments_and_Analysis\Experiments\WKS024\10x';
well_folder = fullfile(root, well);

if ~exist('allData','var')
    allData = struct;
end
if ~exist('T','var')
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
end

scale = 1; % Do not convert pixels to um
if ~isfield(allData, well)
    disp(['Loading the graph of well ', well, '...'])
    allData = update_all_data(allData, well, well_folder, T, scale, network_specifier);
    disp(['Loading the image for well ', well, '...'])
    [fused, cmap] = imread(fullfile(well_folder,[well,'_fused_RGB.tif']));
    allData.(well).fused = fused;
    allData.(well).cmap = cmap;
end

% Get image and graph of current well
measurementNames = allData.measurementNames;
fused = allData.(well).fused;
cmap = allData.(well).cmap;
G = allData.(well).G;
xNodes = allData.(well).xNodes;
yNodes = allData.(well).yNodes;

%% Get voronoi graph. Takes a few mins.
Gv = voronoi_graph(xNodes, yNodes);
%% Plot voronoi graph
figure()
p = plot(Gv, 'XData', P(:,2), 'YData', P(:,1));
p.NodeCData = centrality(Gv, 'betweenness');
figure()
p = plot(G, 'XData', P(:,2), 'YData', P(:,1));
p.NodeCData = centrality(G, 'betweenness');

%% Calculate centrality
Centralities = struct();
centralityNames = {'degree', 'closeness', 'betweenness'};

for i = 1:length(centralityNames)
    cName = centralityNames{i};
    Centralities.(well).('G').(cName) = calculate_normalized_centrality(G, cName);
    Centralities.(well).('Gv').(cName) = calculate_normalized_centrality(Gv, cName);    
end

%%
figure()
for i = 1:length(centralityNames) 
    cName = centralityNames{i}; 
    c = Centralities.(well).('G').(cName);
    cv = Centralities.(well).('Gv').(cName);
    h = kstest2(c', cv');
    
    subplot(3,3,i)
    p = plot(G, 'XData', P(:,2), 'YData', P(:,1));
    p.NodeCData = c;
    title([cName, ', h=',num2str(h)])
    ylabel('observed')
    
    subplot(3,3,i+3)
    p = plot(Gv, 'XData', P(:,2), 'YData', P(:,1));
    p.NodeCData = cv;
    ylabel('voronoi')
    
    subplot(3,3,i+6)
    histogram(c)
    hold on
    histogram(cv)
    hold off
    legend('observed', 'voronoi')
end


