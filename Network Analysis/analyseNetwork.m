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

% Hello hello I am making a change!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experiment = 'WKS023';
magnification = '10x';
well = 'all';
plotting = 'barplot';
disp(experiment)

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
    
    root = '../Experiments/WKS024/10x';
    
    %wellArray = {'B02', 'B03', 'B04', 'B05', 'B06', 'B07'};
    
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
    
    root = '../Experiments/WKS023/2020-09-09';
    
    %wellArray = {'B03', 'B04', 'B05', ...
    %             'D03', 'D04', 'D05'};
             
    groups(1).('Description') = 'Coating';
    groups(1).('Wells') = {'B03', 'B04', 'B05'};
    groups(2).('Description') = 'No coating';
    groups(2).('Wells') = {'D03', 'D04', 'D05'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load xlsx file with well locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load image and graph if this wasn't done already
if ~exist('T','var')
    well_folder = fullfile(root, well);
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
end

if ~exist('allData','var')
    allData = struct;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis for 1 well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(well, 'all')
    
    % Load data for this well (if it wasn't done already)
    if ~isfield(allData, well)
        allData = update_all_data(allData, well, well_folder, T);
    end
    disp('All data loaded.')

    % Get data of current well
    G = allData.(well).G;
    xNodes = allData.(well).xNodes;
    yNodes = allData.(well).yNodes;

    % well locations
    xc = allData.(well).xc;
    yc = allData.(well).yc;
    diameter = allData.(well).diameter;
    
    %% Centrality measures
    
    % Calculate centralities
    Centralities = struct();
    centralityNames = {'R','degree', 'betweenness', 'closeness', 'pagerank', 'eigenvector'};
    nC = length(centralityNames);
    
    Centralities.(well).('R') = sqrt( (xNodes - xc).^2 + (yNodes - yc).^2 );
    
    for i = 2:nC
        cName = centralityNames{i};
        c = centrality(G, cName);
        Centralities.(well).(cName) = c;
    end
    
    % Plot centralities
    figure()
    for i = 1:nC
        
        cName = centralityNames{i};
        c = Centralities.(well).(cName);
        
        subplot(2,nC,i)
        p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
        p1.NodeCData = c;
        colormap jet
        title(cName)
        
        % If cName = R, then normalize the distance per ring area
        if strcmp(cName, 'R')
            [centers, values, dA] = normalize_radial_distance_per_ring_area(c);
            
            subplot(2,nC,i+nC)
            bar(centers, values ./ dA)
            xlabel('R (normalized by ring area)')  
            continue
        end
        
        subplot(2,nC,i+nC)
        histogram(c)
        xlabel(cName)
        
    end
    
    % save
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Color','w','Units','inches','Position',[1 1 18 7.5])
    figName = fullfile('Figures/AnalyzeWells/',[experiment, '_',magnification, '_', well,'_centralities.png']);
    saveas(gcf, figName)
    disp('Centralities saved')

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis for well group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(well, 'all')

    nGroups = length(groups);
    descriptionArray = cell(1,nGroups);

    for i = 1:nGroups
        disp(i)
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
                allData = update_all_data(allData, well, well_folder, T);
            end
            disp(['Data loaded for well ', well])

            % Get data of current well
            G = allData.(well).G;
            xNodes = allData.(well).xNodes;
            yNodes = allData.(well).yNodes;

            % well locations
            xc = allData.(well).xc;
            yc = allData.(well).yc;
            diameter = allData.(well).diameter;

            numNodes = numnodes(G);

            r = sqrt( (xNodes - xc).^2 + (yNodes - yc).^2 );

            % calculate density
            density(j) = calculate_density(G);
            % shortest paths
            distances = shortest_path_lengths(G);
            meanD(j) = sum(distances,'all') * 2 / (numNodes*(numNodes-1));

            % Calculate centralities
            dc = centrality(G, 'degree');
            bc = centrality(G, 'betweenness');
            cc = centrality(G, 'closeness');

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
        figName = fullfile('Figures/AnalyzeGroups/',[experiment, '_', magnification '_scatterplotGroups.pdf']);
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
        figName = fullfile('Figures/AnalyzeGroups/',[experiment, '_', magnification '_meanDegreeGroups.pdf']);
        saveas(gcf, figName)
        
    end % if plotting = 'scatterplot'
    
end % if wells = 'all'

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

function allData = update_all_data(allData, well, well_folder, T)

    row = find( strcmp(T.well, well) );
    diameter = T.diameter(row);
    scale = 6300 / diameter; % (um / pixel)
    
    allData.(well).diameter = diameter * scale;
    allData.(well).xc = T.xc(row) * scale;
    allData.(well).yc = T.yc(row) * scale;

    [G, xNodes, yNodes] = load_graph(well_folder);
    allData.(well).G = G;
    allData.(well).xNodes = xNodes * scale;
    allData.(well).yNodes = yNodes * scale;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G, xNodes, yNodes] = load_graph(well_folder)

    % This function reads the .csv files created by fiji and converts the 
    % information into a ML graph structure.
    
    % Read csv files with edges & node positions
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [centers, values, dA] = normalize_radial_distance_per_ring_area(c)

    [values, edges] = histcounts(c);
    dr = edges(2) - edges(1);
    dA = pi*dr^2 + 2*pi*dr*edges(2:end);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    
end
