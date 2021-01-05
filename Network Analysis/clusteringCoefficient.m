close all

experiment = 'WKS024';
magnification = '20x';
well = 'D02';
fieldSize = 1104;
network_specifier = '_ml';

%% ------------------------------START CODE--------------------------------

root = fullfile('..','..','Experiments', experiment, magnification);
well_folder = fullfile(root, well);

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
    scale = calculate_scale(magnification, fieldSize);
    allData = update_all_data(allData, well, well_folder, T, scale, network_specifier);
end
disp('All data loaded.')

% Get data of current well
measurementNames = allData.measurementNames;
G = allData.(well).G;
xNodes = allData.(well).xNodes;
yNodes = allData.(well).yNodes;

% well locations
xc = allData.(well).xc;
yc = allData.(well).yc;
diameter = allData.(well).diameter;
area = allData.(well).area;

%%

addpath('BrainConnectivity');
A = adjacency(G);
tic
C = clustering_coef_bu(A);
disp(toc)

%%
figure()
p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',5);

p1.NodeCData = C;
colormap jet
colorbar

%%
flagDirected = 0;
r = assortativity(A,flagDirected);
