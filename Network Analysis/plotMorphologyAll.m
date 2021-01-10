close all
clear all

magnification = 'M20';
fieldSize = 1104;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate scale in um/pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale = calculate_scale(magnification, fieldSize);

%% Plot positions  (mean for all wells)
green = [108,186,116] / 255;
pink = [144, 108, 186] / 255;
yellow = [250,255,164] / 255;
allData = struct;

colorHela = green;
colorCos = pink;

experiments = {'JJ005', 'WKS024'};

figure()
for e = 1:2

    experiment = experiments{e};
    %allData = struct();
    disp(experiment);
    root = fullfile('../../Experiments', experiment, magnification);
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
    
    if strcmp(experiment, 'JJ005')
        wellList = {'C02', 'D02'};
        network_specifier = '';
        color = colorCos;
        leg = 'COS-7';
    elseif strcmp(experiment, 'WKS024')
        wellList = {'B02', 'C02', 'D02'};
        network_specifier = '_ml';
        color = colorHela;
        leg = 'HeLa';
    end
    
    % Calculate the cell density as a function of radius
    nBins = 15;
    nWells = length(wellList);

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
        well_folder = fullfile(root, well);
        allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier);
        
        % well locations
        xc = allData.(experiment).(magnification).(well).xc;
        yc = allData.(experiment).(magnification).(well).yc;
        xNodes = allData.(experiment).(magnification).(well).xNodes - xc; % set center of well to x=0
        yNodes = yc - allData.(experiment).(magnification).(well).yNodes; % set center of well to y=0
        diameter = allData.(experiment).(magnification).(well).diameter;

        A = pi * diameter^2 / 4;
        r = sqrt( xNodes.^2 + yNodes.^2 );
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
    
    % statistical tests
    disp(['r ', experiment])
    alfa = 0.05;
    frequency = sum(rBinned,1)';
    custom_kstest(rBinEdges', rBinCenters', frequency, alfa, 'linear')
    
    disp(['th ', experiment])
    alfa = 0.05;
    frequency = sum(thBinned,1)';
    custom_kstest(thBinEdges', thBinCenters', frequency, alfa, 'uniform')
    
    colorCell = {color,color};
    figure(1)
    subplot(2,1,e)
    yscale = 'default';
    plot_shady_error(rBinCenters', rDensityMean', rDensityStd', colorCell, yscale)
    ylim([0,5.5e-4])
    xlim([min(rBinCenters), max(rBinCenters)])
    ylabel('Density (cells \mum^{-2})')
    if e==2
        xlabel('Radial distance to well center (\mum)')
    end
    legend(leg)

    figure(2)
    subplot(2,1,e)
    orange = [243,146,0] / 255;
    color = {orange, 'r'};
    yscale = 'default';
    plot_shady_error(thBinCenters', thDensityMean', thDensityStd', colorCell, yscale)
    xlim([min(thBinCenters), max(thBinCenters)])
    ylim([0,6e-4])
    xticks([-pi/2,0,pi/2])
    xticklabels({'-\pi/2','0','\pi/2'})
    ylabel('Density (cells \mum^{-2})')
    legend(leg)
    if e==2
        xlabel('Angular position (radians)')
    end
end

figure(1)
set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 6])
figName = fullfile('Figures/Morphologies/',[magnification,'_meanRadialPos.png']);
saveas(gcf, figName)

figure(2)
set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 6])
figName = fullfile('Figures/Morphologies/',[magnification,'_meanAngularPos.png']);
saveas(gcf, figName)

%% Plot morphology distribitions
allAreas = struct;
allCircularities = struct;
allEdgeLengths = struct;
meanEdgeLengths = struct;
%allData = struct;

for e = 1:2
    %allData = struct();
    experiment = experiments{e};
    allAreas.(experiment) = [];
    allCircularities.(experiment) = [];
    allEdgeLengths.(experiment) = [];
    meanEdgeLengths.(experiment) = [];
    root = fullfile('../../Experiments', experiment, magnification);
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
    
    if strcmp(experiment, 'JJ005')
        wellList = {'C02', 'D02'};
        network_specifier = '';
        nBins = 750;
        color = green;
    elseif strcmp(experiment, 'WKS024')
        wellList = {'B02', 'C02', 'D02'};
        network_specifier = '_ml';
        nBins = 50;
        color = pink;
    end
    nWells = length(wellList);
    for w = 1:nWells
        well = wellList{w};
        %well_folder = fullfile(root, well);
        %allData = update_all_data(allData, experiment, magnification,...
        %                      well, well_folder, T, scale, network_specifier);
        allAreas.(experiment) = [allAreas.(experiment); allData.(experiment).(magnification).(well).area];
        allCircularities.(experiment) = [allCircularities.(experiment); allData.(experiment).(magnification).(well).circularity];
        
        % find edge lengths
        G = allData.(experiment).(magnification).(well).G;
        xc = allData.(experiment).(magnification).(well).xc;
        yc = allData.(experiment).(magnification).(well).yc;
        xNodes = allData.(experiment).(magnification).(well).xNodes - xc;
        yNodes = yc - allData.(experiment).(magnification).(well).yNodes;
        
        D = find_euclidian_distances(xNodes, yNodes);
        E = adjacency(G) .* D;
        allEdgeLengths.(experiment) = [allEdgeLengths.(experiment); E(E>0)];
        
        meanEdgelengthArray = NaN(numnodes(G), 1);
        % correlate mean edge length with area
        for n = 1:numnodes(G)
            nbh = neighbors(G,n);
            if any(nbh)
                meanEdgelengthArray(n) = mean( D(n,nbh) );
            end
        end
        meanEdgeLengths.(experiment) = [meanEdgeLengths.(experiment); meanEdgelengthArray];
    end
end

%%
figure(1)
figure(2)

for e = 1:2
    experiment = experiments{e};
    if strcmp(experiment, 'JJ005')
        color = colorCos;
    elseif strcmp(experiment, 'WKS024')
        color = colorHela;
    end
    [fArea, edgesArea] = histcounts(allAreas.(experiment), 'BinWidth', 200, 'normalization', 'probability');
    centersArea = (edgesArea(2:end) + edgesArea(1:end-1)) / 2;
    
    [fCirc, edgeCirc] = histcounts(allCircularities.(experiment), 'BinWidth', 0.05, 'normalization', 'probability');
    centersCirc = (edgeCirc(2:end) + edgeCirc(1:end-1)) / 2;
    
    [fEdge, edgeEdgelength] = histcounts(allEdgeLengths.(experiment), 'BinWidth',2.5, 'Normalization', 'probability');
    centersEdgelength = ( edgeEdgelength(2:end) + edgeEdgelength(1:end-1) ) / 2;

    figure(1)
    subplot(2,2,1)
    fill([centersArea,fliplr(centersArea)], [fArea,zeros(1,length(fArea))], ...
         color, 'FaceAlpha', 0.5, 'edgeColor', color)
    hold on
    
 
    ylabel('Normalised frequency')

    subplot(2,2,2)
    fill([centersCirc,fliplr(centersCirc)], [fCirc,zeros(1,length(fCirc))], ...
         color, 'FaceAlpha', 0.5,  'edgeColor', color)
    hold on
    
    figure(2)
    fill([centersEdgelength,fliplr(centersEdgelength)], [fEdge,zeros(1,length(fEdge))], ...
         color, 'FaceAlpha', 0.5,  'edgeColor', color)
    hold on
    xlabel('Edge length (\mum)')
    ylabel('Normalised frequency')
    

end
figure(1)
subplot(2,2,1)
legend('COS-7', 'HeLa')
xlim([0,1.6e4])
subplot(2,2,2)
ylim([0,0.175])

figure(2)
legend('COS-7', 'HeLa')
set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 3])
figName = fullfile('Figures/Morphologies/',[magnification,'_edgeLengthDistribution.png']);
saveas(gcf, figName)
% subplot(3,1,3)
% xlim([0,200])

figure(1)
%
[c1_AH, c2_AH, R2_AH] = fit_linear_model(sqrt(allAreas.WKS024), meanEdgeLengths.WKS024);
[c1_AC, c2_AC, R2_AC] = fit_linear_model(sqrt(allAreas.JJ005), meanEdgeLengths.JJ005);
[c1_CH, c2_CH, R2_CH] = fit_linear_model(allCircularities.WKS024, meanEdgeLengths.WKS024);
[c1_CC, c2_CC, R2_CC] = fit_linear_model(allCircularities.JJ005, meanEdgeLengths.JJ005);

areaAxis = 0:100:max(allAreas.WKS024);
subplot(2,2,3)
scatter(allAreas.WKS024, meanEdgeLengths.WKS024, 7, 'MarkerFaceColor', colorHela, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
hold on
plot(areaAxis, c1_AH + c2_AH*sqrt(areaAxis), '-', 'Color', 'k')
%scatter(allAreas.JJ005, meanEdgeLengths.JJ005, 7, 'MarkerFaceColor', colorCos, 'MarkerFaceAlpha',0.3, 'MarkerEdgeColor', 'none');
%plot(allAreas.JJ005, c1_AC + c2_AC*allAreas.JJ005, '-', 'Color', colorCos)

%plot([allAreas.WKS024; allAreas.JJ005], c1 + c2*[allAreas.WKS024; allAreas.JJ005], '-');
hold off
%set(gca, 'XScale', 'log')
%legend('Hela', 'COS-7')
xlim([0,1.6e4])
xlabel('Area (\mum^2)')
ylabel('Mean edge length (\mum)')

%xlim([0,16000])

subplot(2,2,4)
scatter(allCircularities.WKS024, meanEdgeLengths.WKS024, 7, 'MarkerFaceColor', colorHela, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
hold on
plot(allCircularities.WKS024, c1_CH + c2_CH*allCircularities.WKS024, '-', 'Color', 'k')
%scatter(allCircularities.JJ005, meanEdgeLengths.JJ005, 7, 'MarkerFaceColor', colorCos, 'MarkerFaceAlpha',0.3, 'MarkerEdgeColor', 'none');
%set(gca, 'YScale', 'log')
xlabel('Circularity')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 10 6])
figName = fullfile('Figures/Morphologies/',[magnification,'_allMorphologies.png']);
saveas(gcf, figName)
%plot(allCircularities.WKS024, c1 + c2*allCircularities.WKS024, '-');
%title(sqrt(R2))



figure()

areaAxis = 0:100:max(allAreas.JJ005);

subplot(1,2,1)
scatter(allAreas.JJ005, meanEdgeLengths.JJ005, 7, 'MarkerFaceColor', colorCos, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
hold on
plot(areaAxis, c1_AC + c2_AC*sqrt(areaAxis), '-', 'Color', 'k')
xlabel('Area (\mum)')
ylim([0,350])
ylabel('Mean edge length (\mum)')

subplot(1,2,2)
scatter(allCircularities.JJ005, meanEdgeLengths.JJ005, 7, 'MarkerFaceColor', colorCos, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
hold on
plot(allCircularities.JJ005, c1_CC + c2_CC*allCircularities.JJ005, '-', 'Color', 'k')
xlabel('Circularity')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 10 3])
figName = fullfile('Figures/Morphologies/',[magnification,'_allMorphologiesCOS.png']);
saveas(gcf, figName)

%%
for e = 1:2
    experiment = experiments{e};
    wellList = fields(allData.(experiment).M20);
    for w = 1:length(wellList)
        well = wellList{w};
        diameter = allData.(experiment).M20.(well).diameter;
        area = allData.(experiment).M20.(well).area;
        
        confluency = 4*sum(area) / (pi * diameter^2);
        fprintf("Exp %s, well %s: confluency = %4.2f, mean area = %4.2f, max area = %4.2f \n", experiment, well, confluency, mean(area), max(area))
    end
end

[hArea, pArea] = kstest2(allAreas.WKS024, allAreas.JJ005);
[hCirc, pCirc] = kstest2(allCircularities.WKS024, allCircularities.JJ005);


%%
function [c1, c2, R2] = fit_linear_model(x, y)

    mdl = fitlm(x,y);
    c1 = mdl.Coefficients.Estimate(1);
    c2 = mdl.Coefficients.Estimate(2);
    R2 = mdl.Rsquared.Ordinary;
end

function D = find_euclidian_distances(xNodes,yNodes)

    [X,Y] = meshgrid(xNodes, yNodes);
    D = sqrt( (X - X').^2 + (Y - Y').^2);
end

function custom_kstest(binEdges, binCenters, observedFrequency, alfa, distribution)

    if size(observedFrequency,1) < size(observedFrequency,2)
        error('Input column vector!')
    end

    varnames = {'edges', 'centers', 'frequency'};
    KS = table(binEdges(2:end), binCenters, observedFrequency,'VariableNames',varnames);
    n = sum(KS.frequency);

    Sn = cumsum(KS.frequency) / n;
    dx = 1/length(KS.frequency);
    
    % linear distribution
    if strcmp(distribution, 'linear')
        F = cumsum(dx:dx:1)' / sum(dx:dx:1);
        disp('H0: Data is linearly distributed')
    elseif strcmp(distribution, 'uniform')
        F = cumsum( zeros(length(KS.frequency),1) + dx );
        disp('H0: Data is uniformly distributed')
    end
    figure()
    plot(KS.centers, Sn, KS.centers, F)
    legend('Observed', 'Expected')
    xlabel('x')
    ylabel('Cumulative distribution')

    difference = abs(Sn - F);
    Dmax = max(difference);
    disp(Dmax)

    Dna = sqrt(-log(alfa/2) / (n));
    disp(Dna)

    if Dmax < Dna
        disp('Failure to reject H0. No significant difference between data an distribution.')
    else
        disp('H0 rejected')
    end
end

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
