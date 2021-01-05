close all
clear all
allData = struct;

green = [108,186,116] / 255;
pink = [144, 108, 186] / 255;
yellow = [250,255,164] / 255;

colorHela = green;
colorCos = pink;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data 10x

experiment = 'WKS024';
disp('Starting 10x analysis...')
network_specifier = '';
root = '../../Experiments/WKS024/M10';

groups(1).('Wells') = {'B02', 'C02', 'D02'};
groups(2).('Wells') = {'B03', 'C03', 'D03'};
groups(3).('Wells') = {'B04', 'C04', 'D04'};
groups(4).('Wells') = {'B05', 'C05', 'D05'};
groups(5).('Wells') = {'B06', 'C06', 'D06'};
groups(6).('Wells') = {'B07', 'C07', 'D07'};

magnification = 'M10';
fieldSize = 1104;
scale = calculate_scale(magnification, fieldSize);

xlsfileName = fullfile(root, 'Well locations.xlsx');
T = readtable(xlsfileName);

nGroups = length(groups);
descriptionArray = cell(1,nGroups);

for i = 1:nGroups
    disp(['Group ',num2str(i)])
    wellArray = groups(i).('Wells');
    
    numCells = zeros(1, length(wellArray));
    confluency = zeros(1, length(wellArray));
    meanDegree = zeros(1, length(wellArray));
    meanCloseness = zeros(1, length(wellArray));
    meanBetweenness = zeros(1, length(wellArray));
    pearson = zeros(1, length(wellArray));
    
%     meanR = zeros(1, length(wellArray));
%     density = zeros(1, length(wellArray));
%     meanD = zeros(1, length(wellArray));
%     connected = zeros(1, length(wellArray));

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
        
        %conComp = conncomp(G);
        %sizeLargestGraph = sum( conComp == mode(conComp) );
        %connected(j) = sizeLargestGraph / numnodes(G);

        %r = sqrt( (xNodes).^2 + (yNodes).^2 );

        % calculate density
        %density(j) = calculate_density(G);
        % shortest paths
        %distances = shortest_path_lengths(G);
        %meanD(j) = sum(distances,'all') * 2 / (numNodes*(numNodes-1));

        % Calculate centralities
        dc = calculate_normalized_centrality(G, 'degree');
        bc = calculate_normalized_centrality(G, 'betweenness');
        cc = calculate_normalized_centrality(G, 'closeness');
        r = sqrt(xNodes.^2 + yNodes.^2);
        pearsonCoef = corrcoef(r, cc);

        % Store parameters
        numCells(j) = numNodes;
        confluency(j) = 4*sum(area) / (pi * diameter^2);
        meanDegree(j) = mean(dc);
        meanCloseness(j) = mean(cc);
        meanBetweenness(j) = mean(bc);
        pearson(j) = pearsonCoef(1,2);
        
    end
    
    groups(i).('NumCells') = numCells;
    groups(i).('Confluency') = confluency;
    groups(i).('MeanDegree') = meanDegree;
    groups(i).('MeanCloseness') = meanCloseness;
    groups(i).('MeanBetweenness') = meanBetweenness;
    groups(i).('Pearson') = pearson;
end

%% Load data 20x
magnification = 'M20';
disp('Starting 20x analysis...')

experimentList = {'WKS024'};
fieldSize = 1104;

AV = struct();
Graphs = struct();

for e = 1:length(experimentList)
    
    experiment = experimentList{e};
    if strcmp(experiment, 'WKS024')
        network_specifier = '_ml';
        wells = {'B02', 'C02', 'D02'};
        scale = calculate_scale(magnification, fieldSize);
        color = colorHela;
                
    elseif strcmp(experiment, 'JJ005')
        network_specifier = '';
        wells = {'C02', 'D02'};
        scale = calculate_scale(magnification, fieldSize);
        color = colorCos;
    end
    
    % Load xlsx file
    root = fullfile('../../Experiments', experiment, magnification);
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);

    for w = 1:length(wells)
        well = wells{w};
        well_folder = fullfile(root, well);

        allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier);
  
        % well locations
        xc = allData.(experiment).(magnification).(well).xc;
        yc = allData.(experiment).(magnification).(well).yc;
        diameter = allData.(experiment).(magnification).(well).diameter;
        area = allData.(experiment).(magnification).(well).area;

        % Get data of current well
        G = allData.(experiment).(magnification).(well).G;
        xNodes = allData.(experiment).(magnification).(well).xNodes - xc;  % set center of well to x=0
        yNodes = yc - allData.(experiment).(magnification).(well).yNodes;  % set center of well to y=0

        d = degree(G);
        allData.(experiment).(magnification).(well).degree = d;
        allData.(experiment).(magnification).(well).closeness = calculate_normalized_centrality(G, 'Closeness');
        allData.(experiment).(magnification).(well).betweenness = calculate_normalized_centrality(G, 'Betweenness');

        % init cell array with empty arrays
        m_d = cell(1, max(d) + 1);
        for i = 0:max(d)
            m_d{i+1} = [];
        end

        % fill array with neighbours of cells with degree n
        for n = 1:numnodes(G)
            deg = degree(G, n);
            nbh = neighbors(G, n);
            m_d{deg+1} = [m_d{deg+1}; degree(G, nbh)];
        end

        % Take the mean of the arrays
        mean_m = zeros(1, max(d) + 1);
        var_m = zeros(1, max(d) + 1);
        num_cells = zeros(1, max(d) + 1);
        for i = 0:max(d)
            num_cells(i+1) = sum(d == i);
            mean_m(i+1) = mean(m_d{i+1});
            var_m(i+1) = std(m_d{i+1});
        end
        
        % Lewis law
        meanA = NaN(1, max(d)+1);
        for i = 0:max(d)
            if sum(d==i) > 3
                meanA(i+1) = mean(area(d == i));
            end
        end

        % filter mean_m with less than 3 cells
        mean_m(num_cells < 3) = NaN;
        var_m(num_cells < 3) = NaN;
        allData.(experiment).(magnification).(well).m_d = m_d;
        allData.(experiment).(magnification).(well).mean_m = mean_m;
        allData.(experiment).(magnification).(well).var_m = var_m;
        allData.(experiment).(magnification).(well).max_d = max(d);
        allData.(experiment).(magnification).(well).mean_A = meanA / mean(area);
    end
end

%% Plot
close all

% Data 10x

for i = 1:nGroups
    
    wellArray = groups(i).('Wells');
    confluency = groups(i).Confluency;
    
    meanDegree = groups(i).MeanDegree;
    meanCloseness = groups(i).MeanCloseness;
    meanBetweenness = groups(i).MeanBetweenness;
    pearson = groups(i).Pearson;
    
    figure(1)
    subplot(2,2,4) % Degree vs confluency
    errorbar(mean(confluency), mean(meanDegree), std(meanDegree), std(meanDegree),...
             std(confluency), std(confluency), 's','Color',colorHela, 'MarkerFaceColor',colorHela);
    hold on
  	
     
    figure(3)
    subplot(2,2,4)
    errorbar(mean(confluency), mean(meanCloseness), std(meanCloseness), std(meanCloseness),...
             std(confluency), std(confluency), 's','Color',colorHela, 'MarkerFaceColor',colorHela);
    hold on
    
    figure(4)
    subplot(2,1,2)
    errorbar(mean(confluency), mean(pearson), std(pearson), std(pearson), ...
             std(confluency), std(confluency), 's','Color',colorHela, 'MarkerFaceColor',colorHela);
    hold on
    
    figure(5)
    subplot(2,2,4)
    errorbar(mean(confluency), mean(meanBetweenness), std(meanBetweenness), std(meanBetweenness),...
             std(confluency), std(confluency), 's','Color',colorHela, 'MarkerFaceColor',colorHela);
    hold on
    
end

figure(1)
subplot(2,2,4) % Degree vs confluency
xlabel('Confluency')
ylabel('Mean degree')
set(gca, 'XScale', 'log')

figure(3)
subplot(2,2,4) % Closeness vs confluency
xlabel('Confluency')
ylabel('Mean closeness')
set(gca, 'XScale', 'log')

figure(4)
subplot(2,1,2)
plot([1e-2,1],[0,0], '--k', 'LineWidth', 1)
xlabel('Confluency')
ylabel('Mean correlation')
set(gca, 'XScale', 'log')

figure(5) % Betweenness vs confluency
subplot(2,2,4)
xlabel('Confluency')
ylabel('Mean betweenness')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% Data 20x
experiment = 'WKS024';

if strcmp(experiment, 'WKS024')
    color = colorHela;
    wells = {'B02', 'C02', 'D02'};
    magnification = 'M20';
elseif strcmp(experiment, 'JJ005')
    color = colorCos;
    wells = {'C02', 'D02'};
    magnification = 'M20';
end

allMaxD = 0;

for w = 1:length(wells)
    well = wells{w};
    current_max_d = allData.(experiment).(magnification).(well).max_d;
    if current_max_d > allMaxD
        allMaxD = current_max_d;
    end
end
allMeanM = NaN(length(wells), allMaxD+1);
allVarM = NaN(length(wells), allMaxD+1);
allMeanA = NaN(length(wells), allMaxD+1);

for w = 1:length(wells)
    well = wells{w};
    mean_m = allData.(experiment).(magnification).(well).mean_m;
    var_m = allData.(experiment).(magnification).(well).var_m;
    mean_A = allData.(experiment).(magnification).(well).mean_A;
    allMeanM(w,:) = [mean_m, NaN(1,allMaxD+1 - length(mean_m))];
    allVarM(w,:) = [var_m, NaN(1,allMaxD+1 - length(var_m))];
    allMeanA(w,:) = [mean_A, NaN(1,allMaxD+1 - length(mean_A))];
end

m = nanmean(allMeanM,1);
vm = nanmean(allVarM, 1);
mA = nanmean(allMeanA,1);
std_m = nanstd(allMeanM,1);
std_A = nanstd(allMeanA,1);

% fit on AV law
dAxis = 0:allMaxD;
mdl = fitlm(dAxis(5:end-1), m(5:end-1) .* dAxis(5:end-1));
c1 = mdl.Coefficients.Estimate(1);
c2 = mdl.Coefficients.Estimate(2);
R2_AW = mdl.Rsquared.Ordinary;

% fit on lewis law
yi = mA(5:end-1);             % data
fi = (dAxis(5:end-1) - 2)/4;  % model prediction
SStot = sum( (yi - mean(yi)).^2 );
SSres = sum( (yi - fi).^2 );
R2_lewis = 1 - SSres/SStot;
R_lewis = corrcoef(yi, fi);
% F = @(a,xdata)a(1) + a(2)*xdata.^a(3);
% a0 = [1,1,2];
%[a,resnorm,~,exitflag,output] = lsqcurvefit(F,a0,dAxis(4:end-1),mA(4:end-1));

figure(1) % Degree
d = allData.(experiment).(magnification).D02.degree;
meanC = mean(d);
medianC = median(d);

G = allData.(experiment).(magnification).D02.G;
xNodes = allData.(experiment).(magnification).D02.xNodes - xc;  % set center of well to x=0
yNodes = yc - allData.(experiment).(magnification).D02.yNodes;  % set center of well to y=0

subplot(2,2,[1,3]) % Graph of well D02
p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',1);
p1.NodeCData = d;
colormap jet
h = colorbar('location','NorthOutside');
h.Label.String = 'Degree';
xlabel('x (mm)')
ylabel('y (mm)')
xticks([-4000,0,4000])
xticklabels({'-4','0','4'})
yticks([-4000,0,4000])
yticklabels({'-4','0','4'})

subplot(2,2,2) % Histogram of well D02
histogram(d, 'Normalization','probability','FaceColor',color, 'LineStyle', 'none')
xlabel('Degree, n')
hold on
xline(mean(d), 'r', 'LineWidth', 1.5);
ylabel('Frequency')
leg = legend('HeLa', 'Mean', 'Location', 'NorthEast');
%leg.ItemTokenSize = [15, 1];

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 4])
figName = fullfile('Figures/Centralities/', 'degree.png');
saveas(gcf, figName)

figure(2) % AV law and Lewis's law

subplot(2,1,1) % AV law
plot(dAxis(4:end), c1 + c2*dAxis(4:end), 'Color', 'k', 'LineWidth', 1)
hold on
errorbar(dAxis(4:end), m(4:end) .* dAxis(4:end), std_m(4:end), 's', 'MarkerFaceColor', color,...
         'MarkerEdgeColor', color, 'Color', color, 'LineWidth', 1.5)

ylabel(sprintf('$m_n \\cdot n$'), 'Interpreter', 'latex')
xlim([3.7,11.3])
ylim([15,65])
leg = legend('Linear fit', 'HeLa', 'Location', 'NorthWest');
leg.ItemTokenSize = [15, 1];

subplot(2,1,2) % Lewis' law
plot(dAxis(4:end), (dAxis(4:end) - 2)/4, '-k', 'LineWidth', 1)
hold on
errorbar(dAxis(4:end), mA(4:end), std_A(4:end), 's', 'MarkerFaceColor', color,...
         'MarkerEdgeColor',  color, 'Color', color)

leg = legend("Lewis's law",'HeLa', 'Location', 'NorthWest');
leg.ItemTokenSize = [15, 1];
xlim([3.7,11.3])
ylim([0.8, 2.6])
xlabel('Degree, n')
ylabel(sprintf('$\\bar{A_n} / \\bar{A}$'), 'interpreter', 'latex')
%plot(dAxis(4:end-1), F(a,dAxis(4:end-1)), '-b')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 5 4])
figName = fullfile('Figures/Centralities/', 'AV_lewis.pdf');
saveas(gcf, figName)

figure(3) % Closeness

cc = allData.(experiment).(magnification).D02.closeness;
meanCC = mean(cc);

G = allData.(experiment).(magnification).D02.G;
xNodes = allData.(experiment).(magnification).D02.xNodes - xc;  % set center of well to x=0
yNodes = yc - allData.(experiment).(magnification).D02.yNodes;  % set center of well to y=0
r = sqrt(xNodes.^2 + yNodes.^2);

% fit line on cc versus r
mdl = fitlm(r, cc);
c1 = mdl.Coefficients.Estimate(1);
c2 = mdl.Coefficients.Estimate(2);
R2 = mdl.Rsquared.Ordinary;
disp(['hellooo R2=',num2str(R2)])

subplot(2,2,[1,3])
p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',1);
p1.NodeCData = cc;
colormap jet
h = colorbar('location','NorthOutside');
h.Label.String = 'closeness';
xlabel('x (mm)')
ylabel('y (mm)')
xticks([-4000,0,4000])
xticklabels({'-4','0','4'})
yticks([-4000,0,4000])
yticklabels({'-4','0','4'})
%ylabel(h, 'closeness')

subplot(2,2,2) % histogram closeness
histogram(cc, 'Normalization','probability','FaceColor',colorHela, 'LineStyle', 'none')
xlabel('Closeness')
hold on
xline(mean(cc), 'r', 'LineWidth', 1.5);
ylabel('Frequency')
leg = legend('HeLa', 'Mean', 'Location', 'NorthWest');
ylim([0,0.075])
%leg.ItemTokenSize = [15, 1];

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 4])
figName = fullfile('Figures/Centralities/','closeness.png');
saveas(gcf, figName)

figure(4) % closeness versus radial distance
subplot(2,1,1)
plot(r,cc,'.', 'Color', colorHela)
hold on
plot(r, c1 + c2*r, '-k', 'LineWidth', 1.5)
xlabel('Radial distance to well center (\mu m)')
ylabel('Closeness')
xlim([0,max(r)])
%ylim([-0.1e-6, 3.5e-6])
%leg = legend('HeLa', 'Linear fit', 'Location', 'SouthWest');
%leg.ItemTokenSize = [15, 1];

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 5 4])
figName = fullfile('Figures/Centralities/','closeness_vs_r.png');
saveas(gcf, figName)

figure(5) % betweenness

bc = allData.(experiment).(magnification).D02.betweenness;

subplot(2,2,[1,3])
p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
p1.NodeCData = bc;
colormap jet
h = colorbar('location','NorthOutside');
h.Label.String = 'betweenness';
xlabel('x (mm)')
ylabel('y (mm)')
xticks([-4000,0,4000])
xticklabels({'-4','0','4'})
yticks([-4000,0,4000])
yticklabels({'-4','0','4'})
%ylabel(h, 'closeness')

%

subplot(2,2,2) % histogram betweenness

% Fit on betweenness centrality
[hc, binEdges] = histcounts(bc, 'NumBins', 74, 'Normalization', 'probability');
binCenters = (binEdges(2:end) + binEdges(1:end-1) ) / 2;

F = @(a,xdata)a(1)*xdata.^(-a(2));
a0 = [1,1];
[a,resnorm,~,exitflag,output] = lsqcurvefit(F,a0,binCenters,hc);

histogram(bc, 74,'Normalization','probability','FaceColor',colorHela, 'LineStyle', 'none')
xlabel('Betweenness')
hold on
%xline(mean(bc), 'r', 'LineWidth', 1.5);
ylabel('Frequency')
ylim([1e-5 ,5])
xmin = 0;
xmax = 0.18;
xlim([xmin, xmax])
set(gca, 'YScale', 'log')
yticks([0, 1e-4, 1e-2, 1e0])
% fit exponential distribution
xaxis = linspace(xmin, xmax, 1000);
% pd = fitdist(bc, 'exponential');
% y = pdf(pd, xaxis);
plot(xaxis, F(a, xaxis), '-k', 'LineWidth', 1)
leg = legend('HeLa', 'Power law fit', 'Location', 'NorthEast');
%leg.ItemTokenSize = [15, 1];
pearson_bc = corrcoef(hc, F(a,binCenters));

%
[h,p,st] = chi2gof(binCenters,'Ctrs',binCenters,...
                              'Frequency',hc*length(bc), ...
                              'Expected',F(a,binCenters)*length(bc),...
                              'NParams',1);
%

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 4])
figName = fullfile('Figures/Centralities/','betweenness.png');
saveas(gcf, figName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate scale in um/pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale = calculate_scale(magnification, fieldSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load xlsx file with well locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% root = fullfile('../../Experiments', experiment, magnification);
% well_folder = fullfile(root, well);
% 
% % Load image and graph if this wasn't done already
% if ~exist('T','var')
%     xlsfileName = fullfile(root, 'Well locations.xlsx');
%     T = readtable(xlsfileName);
% end
% 
% if ~exist('allData','var')
%     allData = struct;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis for 1 well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Load data for this well (if it wasn't done already)
% if ~isfield(allData, well)
%     allData = update_all_data(allData, experiment, magnification,...
%                               well, well_folder, T, scale, network_specifier);
% end
% disp('All data loaded.')
% 
% % well locations
% xc = allData.(well).xc;
% yc = allData.(well).yc;
% diameter = allData.(well).diameter;
% area = allData.(well).area;
% circularity = allData.(well).circularity;
% 
% % Get data of current well
% G = allData.(well).G;
% xNodes = allData.(well).xNodes - xc; % set center of well to x=0
% yNodes = yc - allData.(well).yNodes; % set center of well to y=0
% 
% 
% %% Centrality measures
% 
% % Calculate centralities
% NetworkProperties = struct();
% NetworkCentralities = struct();
% propertyNames = {'rDensity', 'theta'};
% centralityNames = {'degree', 'betweenness', 'closeness'};
% 
% % Calculate the cell density as a function of radius
% nBins = 20;
% r = sqrt( xNodes.^2 + yNodes.^2 );
% binEdges = linspace(0, max(r), nBins+1);
% binArea = zeros(1, nBins);
% rBinned = zeros(1, nBins);
% for i = 1:nBins
%     rBinned(i) = sum(r >= binEdges(i) & r<binEdges(i+1));
%     binArea(i) = pi * ( binEdges(i+1)^2 - binEdges(i)^2 ); 
% end
% 
% binCenters = (binEdges(2:end) + binEdges(1:end-1)) / 2;
% NetworkProperties.(well).('rDensity') = rBinned ./ binArea;  
% NetworkProperties.(well).('theta') = atan2(yNodes, xNodes);
% 
% % Calculate network centralities
% nC = length(centralityNames);
% for i = 1:nC
%     cName = centralityNames{i};
%     c = calculate_normalized_centrality(G, cName);
%     NetworkCentralities.(well).(cName) = c;
% end
% 
% % Plot biophysical properties (R, theta)
% figure()
% subplot(1,2,1)
% bar(binCenters, NetworkProperties.(well).('rDensity'),'FaceColor','k');
% xlabel('Radial distance to well center (um)')
% ylabel('Density (# cells / um^2)')
% 
% subplot(1,2,2)
% histogram(NetworkProperties.(well).('theta'), nBins, 'Normalization', 'probability','FaceColor','k');
% xlabel('Angular position (radians)')
% ylabel('Frequency')
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'Color','w','Units','inches','Position',[1 1 8 4])
% figName = fullfile('Figures/Centralities/',[experiment, '_',magnification, '_', well,'_properties.png']);
% saveas(gcf, figName)
% 
% % Plot centrality measures
% figure()
% for i = 1:nC
% 
%     cName = centralityNames{i};
%     c = NetworkCentralities.(well).(cName);
%     meanC = mean(c);
%     medianC = median(c);
% 
%     subplot(2,nC,i)
%     p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
%     p1.NodeCData = c;
%     colormap jet
%     title(cName)
%     xlabel('x (\mum)')
%     ylabel('y (\mum)')
% 
%     subplot(2,nC,i+nC)
%     histogram(c, 'Normalization','probability','FaceColor','k', 'LineStyle', 'none')
%     xlabel(cName)
%     hold on
%     xline(meanC, 'r', 'LineWidth', 1.5);
%     %xline(medianC, 'g', 'LineWidth', 2);
%     if i==1
%         ylabel('Relative frequency')
%     end
%     if i==nC
%         legend('Data', 'Mean', 'Location', 'NorthWest')
%     end
% 
% end
% 
% % save
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'Color','w','Units','inches','Position',[1 1 12 7.5])
% figName = fullfile('Figures/Centralities/',[experiment, '_',magnification, '_', well,'_centralities.png']);
% saveas(gcf, figName)
% %disp('Centralities saved')
% 
% % Calculate distributions of distances
% [X,Y] = meshgrid(xNodes, yNodes);
% distances = sqrt( (X-X').^2 + (Y-Y').^2);
% edgeLength = distances .* adjacency(G);
% existingEdgeLength = edgeLength(edgeLength > 0);
% 
% % Plot distributions of distances
% figure()
% subplot(1,2,1)
% histogram(existingEdgeLength,'Normalization','probability','FaceColor','k')
% xlabel('Edge length (um)')
% ylabel('Probability')
% 
% subplot(1,2,2)
% histogram(distances(distances>0),'Normalization','probability','FaceColor','k')
% xlabel('Distances (um)')
% ylabel('Probability')
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'Color','w','Units','inches','Position',[1 1 8 4])
% figName = fullfile('Figures/Centralities/',[experiment, '_',magnification, '_', well,'_distanceDistributions.png']);
% saveas(gcf, figName)
% 
% % Calculate confluency
% confluency = 4*sum(area) / (pi * diameter^2);
% disp(confluency)
% 
% %% Closeness centrality
% figure()
% cc = NetworkCentralities.(well).(cName);
% % fit line on cc versus r
% mdl = fitlm(r, cc);
% c1 = mdl.Coefficients.Estimate(1);
% c2 = mdl.Coefficients.Estimate(2);
% R2 = mdl.Rsquared.Ordinary;
% disp(['hellooo R2=',num2str(R2)])
% 
% subplot(3,2,[1,3])
% p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
% p1.NodeCData = cc;
% colormap jet
% h = colorbar('location','NorthOutside');
% h.Label.String = 'closeness';
% xlabel('x (\mum)')
% ylabel('y (\mum)')
% %ylabel(h, 'closeness')
% 
% subplot(3,2,[2,4])
% histogram(cc, 'Normalization','probability','FaceColor','k', 'LineStyle', 'none')
% xlabel('Closeness')
% hold on
% xline(mean(cc), 'r', 'LineWidth', 1.5);
% ylabel('Relative frequency')
% leg = legend('Data', 'Mean', 'Location', 'NorthWest');
% 
% subplot(3,2,5:6)
% plot(r,cc,'.k')
% hold on
% plot(r, c1 + c2*r, '-r', 'LineWidth', 1.5)
% xlabel('Radial distance to well center (\mu m)')
% ylabel('Closeness')
% xlim([0,3300])
% leg = legend('Data', 'Linear fit');
% leg.ItemTokenSize = [15, 1];
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'Color','w','Units','inches','Position',[1 1 10 8])
% figName = fullfile('Figures/Centralities/',[experiment, '_',magnification, '_', well,'_closeness.png']);
% saveas(gcf, figName)
% 
% %% Betweenness centrality
% 
% bc = NetworkCentralities.(well).('betweenness');
% 
% figure()
% subplot(1,2,1)
% p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',4);
% p1.NodeCData = bc;
% colormap jet
% h = colorbar('location','NorthOutside');
% h.Label.String = 'Betweenness';
% xlabel('x (\mum)')
% ylabel('y (\mum)')
% 
% subplot(1,2,2)
% histogram(bc, 'Normalization','probability','FaceColor','k', 'LineStyle', 'none')
% xlabel('Betweenness')
% hold on
% xline(mean(bc), 'r', 'LineWidth', 1.5);
% ylabel('Relative frequency')
% set(gca, 'YScale', 'log')
% leg = legend('Data', 'Mean', 'Location', 'NorthEast');
% 
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'Color','w','Units','inches','Position',[1 1 10 5])
% figName = fullfile('Figures/Centralities/',[experiment, '_',magnification, '_', well,'_betweenness.png']);
% saveas(gcf, figName)
% %%
% meanB = zeros(1,numnodes(G));
% for i = 1:numnodes(G)
%     nbh = neighbors(G, i);
%     meanB(i) = mean( bc(nbh) );
% end
% 
% figure()
% plot(bc(bc>0.1), cc(bc>0.1), '.k')

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

function meanA = lewis_law(degreeAxis)
    meanA = (degreeAxis - 2) / 4;
    %meanA = (degreeAxis / 6).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

