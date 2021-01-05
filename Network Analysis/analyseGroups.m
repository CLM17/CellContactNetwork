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
magnification = 'M10';
plotting = 'scatterplot';
fieldSize = 1104;
network_specifier = '';

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
    
    root = '../../Experiments/WKS024/M10';
        
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
    confluency = zeros(1, length(wellArray));
    meanDegree = zeros(1, length(wellArray));
    meanCloseness = zeros(1, length(wellArray));
    meanBetweenness = zeros(1, length(wellArray));
    meanR = zeros(1, length(wellArray));
    density = zeros(1, length(wellArray));
    meanD = zeros(1, length(wellArray));
    connected = zeros(1, length(wellArray));

    for j = 1:length(wellArray)
        well = wellArray{j};
        well_folder = fullfile(root, well);

        % Load data for this well (if it wasn't done already)
        allData = update_all_data(allData, experiment, magnification,...
                          well, well_folder, T, scale, network_specifier);

        disp(['Data loaded for well ', well])

       % well locations

        xc = allData.(experiment).(magnification).(well).xc;
        yc = allData.(experiment).(magnification).(well).yc;
        diameter = allData.(experiment).(magnification).(well).diameter;

        % Get data of current well
        G = allData.(experiment).(magnification).(well).G;
        xNodes = allData.(experiment).(magnification).(well).xNodes - xc; % set center of well to x=0
        yNodes = yc - allData.(experiment).(magnification).(well).yNodes; % set center of well to y=0
        area = allData.(experiment).(magnification).(well).area;

        numNodes = numnodes(G);
        conComp = conncomp(G);
        sizeLargestGraph = sum( conComp == mode(conComp) );
        connected(j) = sizeLargestGraph / numnodes(G);

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
        confluency(j) = 4*sum(area) / (pi * diameter^2);
        meanDegree(j) = mean(dc);
        meanCloseness(j) = mean(cc);
        meanBetweenness(j) = mean(bc);
        meanR(j) = mean(r);

    end

    groups(i).('NumCells') = numCells;
    groups(i).('Confluency') = confluency;
    groups(i).('MeanDegree') = meanDegree;
    groups(i).('MeanCloseness') = meanCloseness;
    groups(i).('MeanBetweenness') = meanBetweenness;
    groups(i).('MeanR') = meanR;
    groups(i).('Density') = density;
    groups(i).('AvgPathLength') = meanD;
   	groups(i).('Connected') = connected;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define plot measures and labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colors = {'b','r','b','r','b','r'};

%measures = {'Confluency','NumCells', 'MeanDegree', 'MeanCloseness', 'MeanBetweenness', 'MeanR', 'Density', 'AvgPathLength'};
measures = {'Connected', 'MeanDegree', 'MeanCloseness', 'MeanBetweenness', 'AvgPathLength'};
%ylabels = {'Confluency', 'Number of cells', 'Mean degree', 'Mean Closeness', 'Mean betweenness', 'Mean radial distance to center', 'Graph density', 'Average path length'};
ylabels = {'Size of main graph / |V|', 'Mean degree', 'Mean closeness', 'Mean betweenness', 'Average path length'};
onLogScale = {'MeanBetweenness'};
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
    figName = fullfile('Figures/Densities/',[experiment, '_', magnification, '_barplotGroups.pdf']);
    saveas(gcf, figName)        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot measures vs number of cells if plotting = 'numberOfCells'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(plotting, 'scatterplot')

    figure(1)
    figure(2)
    colors = lines(nGroups);

    for i = 1:nGroups
        wellArray = groups(i).('Wells');
        connected = groups(i).('Connected');
        confluency = groups(i).Confluency;
        meanConfluency = mean(confluency);
        stdConfluency = std(confluency);
        c = colors(i,:);

        for j = 1:nMeasures
            figure(1)
            subplot(2,ceil(nMeasures/2),j)

            measure = measures{j};
            data = groups(i).(measure);

            errorbar(meanConfluency, mean(data), std(data), std(data),...
                std(confluency), std(confluency), 's','Color','k', 'MarkerFaceColor','k');
            hold on
            
        end
    end

    for j = 1:nMeasures
        measure = measures{j};
        figure(1)
        subplot(2,ceil(nMeasures/2),j)
        set(gca, 'XScale', 'log')
        if any(strcmp(onLogScale, measure))
            set(gca, 'YScale', 'log')
        end
        ylabel(ylabels{j})
        
        if j==nMeasures
            xlabel('Confluency')
        end
    end

    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Color','w','Units','inches','Position',[1 1 6 4])
    figName = fullfile('Figures/Densities/',[experiment, '_', magnification '_scatterplotGroups.png']);
    saveas(gcf, figName)

    %% plot one figure with most important measure
    figure()

    measure = 'AvgPathLength';
    c = 'r';
    maxY = 0;
    maxX = 0;

    for i = 1:nGroups
        wellArray = groups(i).('Wells');
        numcells = groups(i).('NumCells');
        confluency = groups(i).Confluency;
        meanNumCells = mean(numcells);
        meanConfluency = mean(confluency);
        stdConfluency = std(confluency);

        data = groups(i).(measure);

        errorbar(log(meanNumCells), mean(data), std(data), std(data),...
            log(std(numcells)), log(std(numcells)), 's','Color',c, 'MarkerFaceColor',c);
        hold on

        if mean(data) > maxY
            maxY = mean(data);
        end

    end

    %set(gca, 'XScale', 'log');
    
    xlabel('Num cells')
    ylabel('AvgPathLength')
    ylim([0,maxY + 1])

    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Color','w','Units','inches','Position',[1 1 6 4])
    figName = fullfile('Figures/Densities/',[experiment, '_', magnification '_meanDegreeGroups.png']);
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
    %maxD = max(d(d~=inf));
    %d(d==inf) = maxD + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
