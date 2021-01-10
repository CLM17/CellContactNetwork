clear all
close all

experiment = 'WKS028';
magnification = 'M20';
well = 'B03';               % well name
fieldsize = 1104;
scale = calculate_scale(magnification, fieldsize);

experimentList = {'WKS028',  'WKS030', 'WKS029', 'WKS031'};
Description = struct;
Description.WKS028.neuronType = 'Hippocampal';
Description.WKS028.time = 'DIV14';
Description.WKS029.neuronType = 'Cortical';
Description.WKS029.time = 'DIV14';
Description.WKS030.neuronType = 'Hippocampal';
Description.WKS030.time = 'DIV21';
Description.WKS031.neuronType = 'Cortical';
Description.WKS031.time = 'DIV21';


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
    T = readtable(xlsfileName);
    
    for w = 1:nWells 
        well = wellList{w};
        wellFolder = fullfile(root, well);
        allData = load_neuronal_measurements(allData, root, experiment, magnification, T, scale, well);
    end
end

%% positional measures

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

numNeurons = zeros(nWells, nExp);
areaMap2 = zeros(nWells, nExp);
areaAtub = zeros(nWells, nExp);
labels = cell(1,nExp);

for e = 1:nExp
    experiment = experimentList{e};
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

figure()
b = bar(mean(map2Density,1), 'grouped', 'FaceColor', 'flat', 'LineStyle', 'none');
hold on
errorbar(mean(map2Density,1), std(map2Density,[],1), 'LineStyle', 'none')

xticklabels(labels)
xtickangle(45)

%% Functions

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