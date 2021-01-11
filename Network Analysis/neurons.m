clear all
close all

magnification = 'M20';
fieldsize = 1104;
scale = calculate_scale(magnification, fieldsize);

experimentList = {'WKS028',  'WKS030', 'WKS029', 'WKS031'};
Description = struct;
Description.WKS028.neuronType = 'Hippocampal';
Description.WKS028.time = 'DIV14';
Description.WKS028.wellList = {'B03','C03','D03'};
Description.WKS029.neuronType = 'Cortical';
Description.WKS029.time = 'DIV14';
Description.WKS029.wellList = {'B03','C03','D03'};
Description.WKS030.neuronType = 'Hippocampal';
Description.WKS030.time = 'DIV21';
Description.WKS030.wellList = {'B03','C03','D03'};
Description.WKS031.neuronType = 'Cortical';
Description.WKS031.time = 'DIV21';
Description.WKS031.wellList = {'C03','D03'};


nExp = length(experimentList);
%% Load data

wellList = {'B03','C03','D03'};
nWells = length(wellList);

allData = struct;
for e = 1:nExp
    experiment = experimentList{e};
    root = fullfile('../../Experiments', experiment, magnification);
    % Read well locations
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    wellList = Description.(experiment).wellList;
    nWells = length(wellList);
    T = readtable(xlsfileName);
    for w = 1:nWells 
        well = wellList{w};
        wellFolder = fullfile(root, well);
        allData = load_neuronal_measurements(allData, root, experiment, magnification, T, scale, well);
    end
end

%% positional measures

experiment = 'WKS028';
% Calculate the cell density as a function of radius
nBins = 15;

% Bin edges (r and theta)
rBinEdges = linspace(0, 3150, nBins+1);
thBinEdges = linspace(-pi, pi, nBins+1);
rBinCenters = (rBinEdges(2:end) + rBinEdges(1:end-1)) / 2;
thBinCenters = (thBinEdges(2:end) + thBinEdges(1:end-1)) / 2;

rBinArea = zeros(nWells, nBins);
thBinArea = zeros(nWells, nBins);
rBinned = zeros(nWells, nBins);
thBinned = zeros(nWells, nBins);

for w = 1:nWells

    well = wellList{w};
    
    xc = allData.(experiment).(magnification).(well).xc;
    yc = allData.(experiment).(magnification).(well).yc;
    diameter = allData.(experiment).(magnification).(well).diameter;
    xNodes = allData.(experiment).(magnification).(well).xNodes - xc; % set center of well to x=0
    yNodes = yc - allData.(experiment).(magnification).(well).yNodes; % set center of well to y=0
    
    figure()
    plot(xNodes, yNodes, 'k.')
    
    r = sqrt(xNodes.^2 + yNodes.^2);
    A = pi * diameter^2 / 4;
    th = atan2(yNodes, xNodes);

    for i = 1:nBins
        rBinned(w,i) = sum(r >= rBinEdges(i) & r < rBinEdges(i+1));
        thBinned(w,i) = sum(th >= thBinEdges(i) & th < thBinEdges(i+1));
        rBinArea(w,i) = pi * ( rBinEdges(i+1)^2 - rBinEdges(i)^2 );
        thBinArea(w,i) = ( (thBinEdges(i+1) - thBinEdges(i)) / (2*pi) ) * A;
    end
end

rDensityMean = mean(rBinned ./ rBinArea, 1);
rDensityStd = std(rBinned ./ rBinArea, 0, 1);

thDensityMean = mean(thBinned ./ thBinArea, 1);
thDensityStd = std(thBinned ./ thBinArea, 0, 1);
   
figure()
subplot(2,1,1)
yscale = 'default';
color = [108,186,116] / 255;
colorCell = {color,color};
plot_shady_error(rBinCenters', rDensityMean', rDensityStd', colorCell, yscale);

subplot(2,1,2)
plot_shady_error(thBinCenters', thDensityMean', thDensityStd', colorCell, yscale)

%% Num neurons, axon density

numNeurons = NaN(3, nExp);
areaMap2 = NaN(3, nExp);
areaAtub = NaN(3, nExp);
labels = cell(1,nExp);

for e = 1:nExp
    experiment = experimentList{e};
    wellList = Description.(experiment).wellList;
    nWells = length(wellList);
    labels{e} = [Description.(experiment).neuronType, ', ', Description.(experiment).time];
    for w = 1:nWells
        well = wellList{w};
        numNeurons(w,e) = allData.(experiment).(magnification).(well).numNeurons;
        areaMap2(w,e) = allData.(experiment).(magnification).(well).areaMap2;
        areaAtub(w,e) = allData.(experiment).(magnification).(well).areaAtub;
    end
end

map2Density = areaMap2 ./ numNeurons;
atubDensity = areaAtub ./ numNeurons;

% Calculate mean of data and order in matrix
meanNumNeurons = [nanmean(numNeurons(:,1:2),1); nanmean(numNeurons(:,3:4),1)];
stdNumNeurons =  [nanstd(numNeurons(:,1:2),[],1); nanstd(numNeurons(:,3:4),[],1)];

meanMap2Density = [nanmean(map2Density(:,1:2),1); nanmean(map2Density(:,3:4),1)];
stdMap2Density =  [nanstd(map2Density(:,1:2),[],1); nanstd(map2Density(:,3:4),[],1)];

meanAtubDensity = [nanmean(atubDensity(:,1:2),1); nanmean(atubDensity(:,3:4),1)];
stdAtubDensity =  [nanstd(atubDensity(:,1:2),[],1); nanstd(atubDensity(:,3:4),[],1)];

%% ANOVA tests (DIV14 versus DIV21)
[p, tbl, stats] = anova1(map2Density(:,3:4));
[c,~,~,gnames] = multcompare(stats);

%% Plot results in barplot


orange = [243,146,0] / 255;
blue = [30,144,255] / 255;
magentha = [204,0,204] / 255;
darkGreen = [59,122,87] / 290;
lightGreen = [176, 191, 26] / 255;
barColors = {darkGreen, lightGreen};
errorbarColor= 'k';

figure()

subplot(3,1,1)
barplot_errobar_scatter(meanNumNeurons, stdNumNeurons, numNeurons, barColors, errorbarColor)
xticklabels({'',''})
legend({'14 days', '21 days'})
ylabel('Number of neurons')
ylim([0,1800])

subplot(3,1,2)
barplot_errobar_scatter(meanAtubDensity, stdAtubDensity, atubDensity, barColors, errorbarColor)
xticklabels({'',''})
ylabel('Mean neuronal area (\mum^2)')

subplot(3,1,3)
barplot_errobar_scatter(meanMap2Density, stdMap2Density, map2Density, barColors, errorbarColor)
xticklabels({'Hippocampus', 'Cortex'})
xtickangle(0)
ylabel('Mean dendritic area (\mum^2)')
ylim([0,3500])

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 6.5 6.5])
figName = fullfile('Figures/Neurons/','barplots.pdf');
saveas(gcf, figName)  


%% Functions

function barplot_errobar_scatter(y, sigma, scatterdata, barColors, errorbarColor)

    % Plot barplot with errorbars and individual data points.
    
    % INPUTS
    % -------
    % y: ngroups x nbars matrix with mean data (=height of barplots).
    %
    % sigma: ngroups x nbars array with standard deviation (height of
    % errobars).
    %
    % scatterdata: ndata x (ngroups*nbars) data with individual data points.
    %
    % barColors: 1 x nbars cell array with bar colors.
    %
    % errorbarColor: color of error bar.
    %
    % ---------------------------------------------------------------------
    

    b = bar(y, 'grouped', 'FaceColor', 'flat', 'LineStyle', 'none', 'FaceAlpha', 0.5);
    for k = 1:size(y,2)
        b(k).CData = barColors{k};
    end

    hold on
    % Find the number of groups and the number of bars in each group
    ngroups = size(y, 1);
    nbars = size(y, 2);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));

    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    xBars = [];
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        xBars = [xBars, x];
    end

    xBars = sort(xBars);
    scattercolors = repmat(barColors, 1,ngroups);
    for e = 1:ngroups*nbars
        xaxis = xBars(e) + 0.01*randn(size(scatterdata,1),1);
        scatter(xaxis, scatterdata(:,e), 'MarkerEdgeColor', scattercolors{e}, 'MarkerFaceColor', scattercolors{e})
    end
    
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, y(:,i), sigma(:,i), errorbarColor, 'linestyle', 'none');
    end
    
    hold off

end

%%

function allData = load_neuronal_measurements(allData, root, experiment, magnification, T, scale, well)

    wellFolder = fullfile(root, well);
    row = find( strcmp(T.well, well) );
    
    if ~isfloat(T.diameter(row))
        diameter = str2double( T.diameter(row) );
    else
        diameter = T.diameter(row);
    end 
    if ~isfloat(T.xc(row))
        xc = str2double( T.xc(row) );
    else
        xc = T.xc(row);
    end
    if ~isfloat(T.yc(row))
        yc = str2double( T.yc(row) );
    else
        yc = T.yc(row);
    end
    
    allData.(experiment).(magnification).(well).diameter = diameter * scale;
    allData.(experiment).(magnification).(well).xc = xc * scale;
    allData.(experiment).(magnification).(well).yc = yc * scale;
    
    fNameCOM = fullfile(wellFolder, 'nuclei_com.csv');
    nuclei_com = csvread(fNameCOM);
    
    allData.(experiment).(magnification).(well).xNodes = nuclei_com(:,2) * scale;
    allData.(experiment).(magnification).(well).yNodes = nuclei_com(:,3) * scale;
    
    % read neuron measurements
    neuronMeasurements = fullfile(wellFolder, [well,'_neuronMeasurements.txt']);
    file = fileread(neuronMeasurements);
    disp(experiment)
    disp(file)
    measurements = split(file, ',');
    
    spl1 = split(measurements{1}, '= ');
    allData.(experiment).(magnification).(well).numNuclei = str2double(spl1{2});
    spl2 = split(measurements{2}, '= ');
    allData.(experiment).(magnification).(well).numNeurons = str2double(spl2{2});
    spl3 = split(measurements{3}, '= ');
    allData.(experiment).(magnification).(well).areaMap2 = str2double(spl3{2}) * scale^2;
    spl4 = split(measurements{4}, '= ');
    allData.(experiment).(magnification).(well).areaAtub = str2double(spl4{2}) * scale^2;
end