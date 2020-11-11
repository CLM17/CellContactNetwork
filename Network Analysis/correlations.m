experiment = 'WKS024';
magnification = '10x';
well = 'D03';

%% ------------------------------START CODE--------------------------------

root = fullfile('..','Experiments', experiment, magnification);
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

%%
%measurementNames = Area	Mode	BX	BY	Width	Height	Major	Minor	Angle	Circ.	AR	Round	Solidity
%cellMeasurements = csvread(fullfile(well_folder,'cell_measurements.csv'), 1);
cellMeasurementsTable = readtable(fullfile(well_folder,'cell_measurements.csv'));

measurementNames = {'area', 'circularity', 'longness'};
nM = length(measurementNames);

cellValues = cellMeasurementsTable.Mode;
[~,ind] = sort(cellValues);

Measurements = struct();
Measurements.(well).('area') = cellMeasurementsTable.Area(ind);
Measurements.(well).('circularity') = cellMeasurementsTable.Circ_(ind);
Measurements.(well).('longness') = cellMeasurementsTable.Major(ind) ./ cellMeasurementsTable.Minor(ind);

%% Calulcate centralities
Centralities = struct();
centralityNames = {'degree', 'betweenness', 'closeness', 'pagerank', 'eigenvector'};
nC = length(centralityNames);

%% Vector dimensions don't agree: please have a look at this later!!
for i = 1:nC
    cName = centralityNames{i};
    c = centrality(G, cName);
    for j = 1:length(c)
        if ~any( ismember(cellValues, j) )
            c(j) = [];
        end
    end
        
    Centralities.(well).(cName) = c;
end

%% Plot cell measures versus centralities
figure()
area = Measurements.(well).area;
circularity = Measurements.(well).circularity;
longness = Measurements.(well).longness;

count = 0;
for i = 1:nM
    
    mName = measurementNames{i};
    m = Measurements.(well).(mName);
    
    for j = 1:nC
        count = count + 1;
        
        cName = centralityNames{j};
        c = Centralities.(well).(cName);

        % fit linear model c = b1 + b2*m
        mdl = fitlm(m, c);
        b1 = mdl.Coefficients.Estimate(1);
        b2 = mdl.Coefficients.Estimate(2);
        R2 = mdl.Rsquared.Ordinary;

        subplot(3, nC, count)
        plot(m, c, '.', m, b1 + b2*m, '-r')
        xlabel(mName)
        ylabel(cName)
        title(['R^2 = ', num2str(R2)])
    end
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 15 9])
figName = fullfile('Figures/Correlations/',[experiment, '_', magnification, '_', well,'_measurementCorrelations.png']);
saveas(gcf, figName)

%% Plot all centralities
figure()
count = 0;
for i = 1:nC
    cName1 = centralityNames{i};
    c1 = Centralities.(well).(cName1);
    for j = 1:nC
        count = count + 1;
        
        cName2 = centralityNames{j};
        c2 = Centralities.(well).(cName2);
        
        % fit linear model c2 = b1 + b2*c1
        mdl = fitlm(c1, c2);
        b1 = mdl.Coefficients.Estimate(1);
        b2 = mdl.Coefficients.Estimate(2);
        R2 = mdl.Rsquared.Ordinary;
        
        subplot(nC, nC, count)
        plot(c1, c2, '.', c1, b1 + b2*c1, '-r')
        xlabel(cName1)
        ylabel(cName2)
        title(['R^2 = ', num2str(R2)])
        
    end
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 15 12])
figName = fullfile('Figures/Correlations/',[experiment, '_', magnification, '_', well,'_centralityCorrelations.png']);
saveas(gcf, figName)

%% ------------------------------FUNCTIONS---------------------------------

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

function [G, xNodes, yNodes] = load_graph(well_folder)

    % This function reads the .csv files created by fiji and converts the 
    % information into a ML graph structure.
    
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