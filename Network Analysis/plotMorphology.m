close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters

%   // experiment: string
%   Experiment name.

%   // well: string
%   Name of well (example: 'B02')

%   // fieldSize: double
%   Size of 1 field (is used to calculate the scale in um/pixel).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experiment = 'WKS024';
magnification = '20x';
well = 'D02';
fieldSize = 1104;
network_specifier = '_ml';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate scale in um/pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale = calculate_scale(magnification, fieldSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load xlsx file with well locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root = fullfile('../../Experiments', experiment, magnification);
well_folder = fullfile(root, well);

% Load image and graph if this wasn't done already
if ~exist('T','var')
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
end

if ~exist('allData','var')
    allData = struct;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis for 1 well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data for this well (if it wasn't done already)
if ~isfield(allData, well)
    allData = update_all_data(allData, well, well_folder, T, scale, network_specifier);
end
disp('All data loaded.')

% well locations
xc = allData.(well).xc;
yc = allData.(well).yc;
diameter = allData.(well).diameter;
area = allData.(well).area;
circularity = allData.(well).circularity;
longness = allData.(well).longness;

% Get data of current well
G = allData.(well).G;
xNodes = allData.(well).xNodes - xc; % set center of well to x=0
yNodes = yc - allData.(well).yNodes; % set center of well to y=0

%% 
% Calculate the cell density as a function of radius
nBins = 40;
A = pi * diameter^2 / 4;
r = sqrt( xNodes.^2 + yNodes.^2 );
th = atan2(yNodes, xNodes);

% Bin edges (r and theta)
rBinEdges = linspace(0, max(r), nBins+1);
thBinEdges = linspace(-pi, pi, nBins+1);

rBinArea = zeros(1, nBins);
thBinArea = zeros(1, nBins);
rBinned = zeros(1, nBins);
thBinned = zeros(1, nBins);
for i = 1:nBins
    rBinned(i) = sum(r >= rBinEdges(i) & r < rBinEdges(i+1));
    thBinned(i) = sum(th >= thBinEdges(i) & th < thBinEdges(i+1));
    rBinArea(i) = pi * ( rBinEdges(i+1)^2 - rBinEdges(i)^2 );
    thBinArea(i) = ( (thBinEdges(i+1) - thBinEdges(i)) / (2*pi) ) * A;
end

rBinCenters = (rBinEdges(2:end) + rBinEdges(1:end-1)) / 2;
thBinCenters = (thBinEdges(2:end) + thBinEdges(1:end-1)) / 2;
allData.(well).('rDensity') = rBinned ./ rBinArea;  
allData.(well).('thDensity') = thBinned ./ thBinArea;

measurementNames = {'rDensity', 'thDensity', 'area', 'circularity'};
numM = length(measurementNames);
titles = {'Radial position', 'Angular position', 'Cell area', 'Circularity'};
xlabels = {'R (\mum)', '\theta (rad)', 'A (\mum^2)', '4\pi\cdotA\cdotperimeter^{-2}'};
ylabels = {'Density (cells \mum^{-2})', 'Density (cells \mum^{-2})', 'Frequency', 'Frequency'};

figure()
% plot r density
subplot(2,2,1)
m = allData.(well).('rDensity');
bar(rBinCenters, m,'FaceColor','k', 'BarWidth', 1, 'FaceAlpha', 0.6, 'LineStyle', 'none')
xlabel(xlabels{1})
ylabel(ylabels{1})
title(titles{1})

% plot theta density
subplot(2,2,2)
m = allData.(well).('thDensity');
bar(thBinCenters, m,'FaceColor','k', 'BarWidth', 1, 'FaceAlpha', 0.6, 'LineStyle', 'none')
xlabel(xlabels{2})
ylabel(ylabels{2})
title(titles{2})
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})

% Plot area
subplot(2,2,3)
m = allData.(well).('area');
histogram(m, nBins, 'Normalization','probability','FaceColor','k', 'LineStyle', 'none')
xlabel(xlabels{3})
ylabel(ylabels{3})
title(titles{3})

% Plot circularity
subplot(2,2,4)
m = allData.(well).('circularity');
histogram(m, nBins, 'Normalization','probability','FaceColor','k', 'LineStyle', 'none')
xlabel(xlabels{4})
ylabel(ylabels{4})
title(titles{4})

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 6])
figName = fullfile('Figures/Morphologies/',[experiment, '_',magnification, '_', well,'_morphologies.png']);
saveas(gcf, figName)

%%
distName = 'gamma';
figure()
weibull = fit_dist(allData.(well).('circularity'), distName);
histfit(allData.(well).('area'), nBins, distName)

%% functions
function dist = fit_dist(data, distname)
% Fits a Weibull distribution to crystal data an performs a chi-squared goodness of fit test. 
% If h=0, the distribution fits the data, and if h=1 the distribution does not fit.
% i is an integer (1 or 2) that denotes the dataset number. 

    pd = fitdist(data, distname);
    [h, p] = chi2gof(data,'CDF',pd);
    
    dist = struct('h', h, 'p', p);
    
     if h==1
        text = sprintf('Does NOT follow a %s distribution, p=%.3f',distname, p);
    else
        text = sprintf('Follows a %s distribution, p=%.3f',distname, p);
     end
    disp(text)
end
