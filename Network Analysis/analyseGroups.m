close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters

%   // experiment: string
%   Experiment name.

%   // well: string
%   Name of well (example: 'B02'), or 'all' if you want to analyze all
%   wells of the experiment (in groups).

%   // plotting: string
%   If well = 'all', define here how you want the output to be plotted.
%   There are two options:
%   (1) plotting = 'barplot' --> the output will be a barplots with one bar per
%       group.
%   (2) plotting = 'scatterplot' --> the output will be scatter plots
%       with one datapoint per group.

%   // fieldSize: double
%   Size of 1 field (is used to calculate the scale in um/pixel).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experiment = 'WKS024';
magnification = '10x';
plotting = 'scatterplot';
fieldSize = 1104;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define groups

%   // groups: structure array.
%   Each entry in the datastructure corresponds to one group.
%   Fields:
%   (1) groups(i).Description: string with a short description of the
%       group.
%   (2) groups(i).Wells: cell array with all wells in the group.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(strcmp(plotting, 'barplot') || strcmp(plotting, 'scatterplot'))
    error('Invalid plotting method. Choose from barplot or scatterplot.')
end

groups = struct();
if strcmp(experiment, 'WKS024')
    
    root = '../../Experiments/WKS024/10x';
        
    groups(1).('Description') = '1';
    groups(1).('Wells') = {'B02', 'C02', 'D02'};
    groups(2).('Description') = '2';
    groups(2).('Wells') = {'B03', 'C03', 'D03'};
    groups(3).('Description') = '4';
    groups(3).('Wells') = {'B04', 'C04', 'D04'};
    groups(4).('Description') = '8';
    groups(4).('Wells') = {'B05', 'C05', 'D05'};
    groups(5).('Description') = '16';
    groups(5).('Wells') = {'B06', 'C06', 'D06'};
    groups(6).('Description') = '32';
    groups(6).('Wells') = {'B07', 'C07', 'D07'};
    
elseif strcmp(experiment, 'WKS023')
    
    root = '../../Experiments/WKS023/2020-09-09';
                 
    groups(1).('Description') = 'Coating';
    groups(1).('Wells') = {'B03', 'B04', 'B05'};
    groups(2).('Description') = 'No coating';
    groups(2).('Wells') = {'D03', 'D04', 'D05'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate scale in um/pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale = calculate_scale(magnification, fieldSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load xlsx file with well locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load image and graph if this wasn't done already
if ~exist('T','var')
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
end

if ~exist('allData','var')
    allData = struct;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis for well group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nGroups = length(groups);
descriptionArray = cell(1,nGroups);

for i = 1:nGroups
    disp(['Group ',num2str(i)])
    wellArray = groups(i).('Wells');
    descriptionArray{i} = groups(i).('Description');

    numCells = zeros(1, length(wellArray));
    meanDegree = zeros(1, length(wellArray));
    meanCloseness = zeros(1, length(wellArray));
    meanBetweenness = zeros(1, length(wellArray));
    meanR = zeros(1, length(wellArray));
    density = zeros(1, length(wellArray));
    meanD = zeros(1, length(wellArray));

    for j = 1:length(wellArray)
        well = wellArray{j};
        well_folder = fullfile(root, well);

        % Load data for this well (if it wasn't done already)
        if ~isfield(allData, well)
            allData = update_all_data(allData, well, well_folder, T, scale);
        end
        disp(['Data loaded for well ', well])

       % well locations

        xc = allData.(well).xc;
        yc = allData.(well).yc;
        diameter = allData.(well).diameter;


        % Get data of current well
        G = allData.(well).G;
        xNodes = allData.(well).xNodes - xc; % set center of well to x=0
        yNodes = yc - allData.(well).yNodes; % set center of well to y=0

        numNodes = numnodes(G);

        r = sqrt( (xNodes).^2 + (yNodes).^2 );

        % calculate density
        density(j) = calculate_density(G);
        % shortest paths
        distances = shortest_path_lengths(G);
        meanD(j) = sum(distances,'all') * 2 / (numNodes*(numNodes-1));

        % Calculate centralities
        dc = calculate_normalized_centrality(G, 'degree');
        bc = calculate_normalized_centrality(G, 'betweenness');
        cc = calculate_normalized_centrality(G, 'closeness');

        % Store parameters
        numCells(j) = numNodes;
        meanDegree(j) = mean(dc);
        meanCloseness(j) = mean(cc);
        meanBetweenness(j) = mean(bc);
        meanR(j) = mean(r);

    end

    groups(i).('NumCells') = numCells;
    groups(i).('MeanDegree') = meanDegree;
    groups(i).('MeanCloseness') = meanCloseness;
    groups(i).('MeanBetweenness') = meanBetweenness;
    groups(i).('MeanR') = meanR;
    groups(i).('Density') = density;
    groups(i).('AvgPathLength') = meanD;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define plot measures and labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colors = {'b','r','b','r','b','r'};

measures = {'NumCells', 'MeanDegree', 'MeanCloseness', 'MeanBetweenness', 'MeanR', 'Density', 'AvgPathLength'};
%measures = {'NumCells', 'MeanDegree'};
ylabels = {'Number of cells', 'Mean degree', 'Mean Closeness', 'Mean betweenness', 'Mean radial distance to center', 'Graph density', 'Average path length'};
%ylabels = {'Number of cells', 'Mean number of neighbours'};
onLogScale = {'NumCells', 'MeanBetweenness'};
nMeasures = length(measures);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot barplots if plotting = 'barplots'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(plotting, 'barplot')

    figure()
    for i = 1:nGroups
        wellArray = groups(i).('Wells');

        for j = 1:nMeasures

            measure = measures{j};
            data = groups(i).(measure);

            subplot(1,nMeasures,j)
            bar(i, mean(data), colors{i}, 'FaceAlpha', 0.3, 'EdgeColor', 'None')
            hold on
            x = zeros(1, length(wellArray)) + i;
            scatter(x, data, ['o',colors{i}],'MarkerFaceColor', colors{i}, 'MarkerFaceAlpha', 0.5)
            hold on
            errorbar(i, mean(data), std(data), colors{i})
        end
    end

    for j = 1:nMeasures
        subplot(1,nMeasures,j)
        ylabel(ylabels{j})
        xticks(1:nGroups)
        xticklabels(descriptionArray)
        xtickangle(45)
    end

    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Color','w','Units','inches','Position',[1 1 12 4])
    figName = fullfile('Figures/AnalyzeGroups/',[experiment, '_', magnification, '_barplotGroups.pdf']);
    saveas(gcf, figName)        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot measures vs number of cells if plotting = 'numberOfCells'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(plotting, 'scatterplot')

    figure()
    colors = lines(nGroups);

    for i = 1:nGroups
        wellArray = groups(i).('Wells');
        numCells = groups(i).NumCells;
        meanNumCells = mean(numCells);
        stdNumCells = std(numCells);
        c = colors(i,:);

        for j = 2:nMeasures
            subplot(2,floor(nMeasures/2),j-1)

            measure = measures{j};
            data = groups(i).(measure);

            errorbar(meanNumCells, mean(data), std(data), std(data),...
                std(numCells), std(numCells), 's','Color',c, 'MarkerFaceColor',c);
            hold on
        end
    end

    for j = 2:nMeasures
        measure = measures{j};
        subplot(2,floor(nMeasures/2),j-1)
        set(gca, 'XScale', 'log')
        if any(strcmp(onLogScale, measure))
            set(gca, 'YScale', 'log')
        end
        ylabel(ylabels{j})
        xlabel('Number of cells')
    end

    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Color','w','Units','inches','Position',[1 1 12 6])
    figName = fullfile('Figures/AnalyzeGroups/',[experiment, '_', magnification '_scatterplotGroups.png']);
    saveas(gcf, figName)

    %% plot one figure with most important measure
    figure()

    measure = 'MeanDegree';
    c = 'r';
    maxY = 0;
    maxX = 0;

    for i = 1:nGroups
        wellArray = groups(i).('Wells');
        numCells = groups(i).NumCells;
        meanNumCells = mean(numCells);
        stdNumCells = std(numCells);

        data = groups(i).(measure);

        errorbar(meanNumCells, mean(data), std(data), std(data),...
            std(numCells), std(numCells), 's','Color',c, 'MarkerFaceColor',c);
        hold on

        if mean(data) > maxY
            maxY = mean(data);
        end
        if meanNumCells > maxX
            maxX = meanNumCells;
        end
    end

    set(gca, 'XScale', 'log');
    xlabel('Number of cells')
    ylabel('Mean number of neighbours')
    ylim([0,maxY + 1])
    xlim([0,maxX + 1e3])

    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Color','w','Units','inches','Position',[1 1 6 4])
    figName = fullfile('Figures/AnalyzeGroups/',[experiment, '_', magnification '_meanDegreeGroups.png']);
    saveas(gcf, figName)

end % if plotting = 'scatterplot'

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

function nc = calculate_normalized_centrality(G, cName)
    
     c = centrality(G, cName);
     N = numnodes(G);
     
     if strcmp(cName, 'closeness')
         nc = c * (N - 1);
     elseif strcmp(cName, 'betweenness')
         nc = 2*c / ( (N-1)*(N-2) );
     else
         nc = c;
     end

end