close all

experiment = 'WKS024';
magnification = 'M20';
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
    allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier);
end
disp('All data loaded.')

% Get data of current well
measurementNames = {'area', 'circularity'};
G = allData.(experiment).(magnification).(well).G;
xNodes = allData.(experiment).(magnification).(well).xNodes;
yNodes = allData.(experiment).(magnification).(well).yNodes;

% well locations
xc = allData.(experiment).(magnification).(well).xc;
yc = allData.(experiment).(magnification).(well).yc;
diameter = allData.(experiment).(magnification).(well).diameter;

%% Calulcate centralities
Centralities = struct();
centralityNames = {'degree', 'betweenness', 'closeness'};
nC = length(centralityNames);

%% Calculate centralities
for i = 1:nC
    cName = centralityNames{i};
    c = calculate_normalized_centrality(G, cName);
    Centralities.(well).(cName) = c;
end

%% Plot cell measures versus centralities
figure()
count = 0;
nM = length(measurementNames);
for i = 1:nM
    
    mName = measurementNames{i};
    m = allData.(experiment).(magnification).(well).(mName);
    
    for j = 1:nC
        count = count + 1;
        
        cName = centralityNames{j};
        c = Centralities.(well).(cName);

        % fit linear model c = b1 + b2*m
        mdl = fitlm(m, c);
        b1 = mdl.Coefficients.Estimate(1);
        b2 = mdl.Coefficients.Estimate(2);
        R2 = mdl.Rsquared.Ordinary;

        subplot(nM, nC, count)
        plot(m, c, '.', m, b1 + b2*m, '-r')
        xlabel(mName)
        ylabel(cName)
        title(sprintf('R^2 = %.2f', R2))
    end
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 12 6])
figName = fullfile('Figures/Correlations',[experiment, '_', magnification, '_', well,'_measurementCentralityCorrelations.png']);
saveas(gcf, figName)

%% Plot all measurements
figure()
count = 0;
for i = 1:nM
    
    mName1 = measurementNames{i};
    m1 = allData.(experiment).(magnification).(well).(mName1);
    
    for j = 1:nM
        count = count + 1;
        
        mName2 = measurementNames{j};
        m2 = allData.(experiment).(magnification).(well).(mName2);

        % fit linear model c = b1 + b2*m
        mdl = fitlm(m1, m2);
        b1 = mdl.Coefficients.Estimate(1);
        b2 = mdl.Coefficients.Estimate(2);
        R2 = mdl.Rsquared.Ordinary;

        subplot(nM, nM, count)
        plot(m1, m2, '.', m1, b1 + b2*m1, '-r')
        xlabel(mName1)
        ylabel(mName2)
        title(sprintf('R^2 = %.2f', R2))
    end
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 6])
figName = fullfile('Figures/Correlations',[experiment, '_', magnification, '_', well,'_measurementCorrelations.png']);
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
        title(sprintf('R^2 = %.2f', R2))
        
    end
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 12 9])
figName = fullfile('Figures/Correlations',[experiment, '_', magnification, '_', well,'_centralityCorrelations.png']);
saveas(gcf, figName)
