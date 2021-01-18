close all
clear all


orange = [243,146,0] / 255;
blue = [30,144,255] / 255;

well = 'D02';               % well name
experiment = 'WKS024';
magnification = 'M20';
fieldSize = 1104;
cutoffDistance = 250; % cutoff in micron
network_specifier = '_ml';

%% ----------------------------- Load data --------------------------------

load('bayesianAnalysis.mat');
magnification = 'M20';

root = fullfile('..','..','Experiments',experiment,magnification);
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

scale = calculate_scale(magnification, fieldSize); % Do not convert pixels to um
disp(scale)
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

%% -------------------------------- Plot ----------------------------------

figure()
%subplot(2,1,1)
[hDistances, eDistances] = histcounts(distances);
[hExistingEdgeLength, eExistingEdgeLength] = histcounts(existingEdgeLength);
%bar(bin_centers(eDistances), hDistances / max(hDistances), ...
%    'barWidth', 1, 'LineStyle', 'None', 'FaceAlpha', 0.7)
centersDist = bin_centers(eDistances);
centersEdge = bin_centers(eExistingEdgeLength);

fill([centersDist,fliplr(centersDist)],... 
     [hDistances / max(hDistances), zeros(1,length(hDistances))], ...
     blue, 'FaceAlpha', 0.5,  'edgeColor', blue)

hold on
% bar(bin_centers(eExistingEdgeLength), hExistingEdgeLength / max(hExistingEdgeLength), ...
%     'barWidth', 1, 'LineStyle', 'None', 'FaceAlpha', 0.7)

fill([centersEdge,fliplr(centersEdge)],... 
     [hExistingEdgeLength / max(hExistingEdgeLength), zeros(1,length(hExistingEdgeLength))], ...
     orange, 'FaceAlpha', 0.5,  'edgeColor', orange)

set(gca, 'XScale', 'log')
ylim([0,1.4])
xlim([min(existingEdgeLength), 16000])

ylabel('Frequency / max frequency')
xlabel('Distance (\mum)')
leg = legend('Distances between all cells', 'Distances between connected cells', 'Location', 'NorthWest');

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 6 2.5])
figName = fullfile('Figures/Bayes',[well,'observedDistributions.pdf']);
saveas(gcf, figName)

%%

figure()
%subplot(2,1,2)
plot(pConnect, '-', 'LineWidth', 1.5, 'Color', orange)
hold on
%maxDist = dAxis(pConnect == max(pConnect));
%line([maxDist, maxDist], [0,max(pConnect)],'LineStyle', '--', 'Color', 'k','LineWidth',2);
%line([0, maxDist], [max(pConnect),max(pConnect)],'LineStyle', '--', 'Color', 'k','LineWidth',2);

xlim([0,250])
ylim([0,1.1])
xlabel('Distance, d_{ij} (\mum)')
ylabel('Probability of connection')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 6 2.5])
figName = fullfile('Figures/Bayes',[well,'inferredDistribution.pdf']);
saveas(gcf, figName)

function centers = bin_centers(edges)
    centers = (edges(2:end) + edges(1:end-1)) / 2;
end

function [distances, existingEdgeLength] = distances_distributions(G, xNodes, yNodes)
    % distance matrix
    [X,Y] = meshgrid(xNodes,yNodes);
    full_matrix = full( adjacency(G) );
    distances = sqrt( (X-X').^2 + (Y-Y').^2);
    edgeLength = distances .* full_matrix;
    existingEdgeLength = edgeLength(edgeLength > 0);
end