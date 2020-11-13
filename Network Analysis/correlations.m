experiment = 'WKS024';
magnification = '10x';
well = 'B03';
fieldSize = 1104;

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
    allData = update_all_data(allData, well, well_folder, T, scale);
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
cellMeasurementsTable = readtable(fullfile(well_folder,'cell_measurements.csv'));

measurementNames = {'area', 'circularity', 'longness'};
nM = length(measurementNames);

cellValues = cellMeasurementsTable.Mode;
[~,ind] = sort(cellValues);

Measurements = struct();
Measurements.(well).('area') = cellMeasurementsTable.Area(ind) * scale^2;
Measurements.(well).('circularity') = cellMeasurementsTable.Circ_(ind);
Measurements.(well).('longness') = cellMeasurementsTable.Major(ind) ./ cellMeasurementsTable.Minor(ind);

%% Calulcate centralities
Centralities = struct();
centralityNames = {'degree', 'betweenness', 'closeness', 'pagerank', 'eigenvector'};
nC = length(centralityNames);

%% Remove nuclei that do not overlap with a cell (i.e. they have no cell measure)
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
set(gcf,'Color','w','Units','inches','Position',[1 1 12 9])
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
set(gcf,'Color','w','Units','inches','Position',[1 1 12 9])
figName = fullfile('Figures/Correlations/',[experiment, '_', magnification, '_', well,'_centralityCorrelations.png']);
saveas(gcf, figName)
