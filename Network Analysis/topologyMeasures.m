close all
clear all
addpath('BrainConnectivity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
magnification = 'M20';
well = 'D02';
fieldSize = 1104;
network_specifier = '_ml';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
green = [108,186,116] / 255;
pink = [144, 108, 186] / 255;
yellow = [250,255,164] / 255;

colorHela = green;
colorCos = pink;

groups = struct();
groups(1).WKS024 = {'B02', 'C02', 'D02'};
groups(2).WKS024 = {'B03', 'C03', 'D03'};
groups(3).WKS024 = {'B04', 'C04', 'D04'};
groups(4).WKS024 = {'B05', 'C05', 'D05'};
groups(5).WKS024 = {'B06', 'C06', 'D06'};
groups(6).WKS024 = {'B07', 'C07', 'D07'};

nGroups = length(groups);
nWells = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate topoology measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathLength = struct;
clusteringCoefficient = struct;
assort = struct;
largestGraph = struct;
confluency = struct;
areas = struct;
diameters = struct;
numNodes = struct;
flag = 0;

experimentList = {'JJ005', 'WKS024'};
for e = 2:2
    experiment = experimentList{e};
    if strcmp(experiment, 'WKS024')
        network_specifier = '';
        
        magnification = 'M10';
        scale = calculate_scale(magnification, fieldSize);
                
    elseif strcmp(experiment, 'JJ005')
        network_specifier = '';
        wellList = {'B03', 'C02', 'D02', 'D03'};
        magnification = 'M20';
        scale = calculate_scale(magnification, fieldSize);
    end
    
    root = fullfile('..','..','Experiments',experiment,magnification);
    well_folder = fullfile(root, well);
    
    allData = struct;
    
    avgPathLength.(experiment).mean = zeros(nGroups,nWells);
    avgPathLength.(experiment).stdv = std(nGroups,nWells);
    clusteringCoefficient.(experiment).mean = zeros(nGroups,nWells);
    clusteringCoefficient.(experiment).stdv = std(nGroups,nWells);
    assort.(experiment).mean = zeros(nGroups,nWells);
    assort.(experiment).stdv = std(nGroups,nWells);
    largestGraph.(experiment).mean = zeros(nGroups,nWells);
    confluency.(experiment) = zeros(nGroups,nWells);
    largestGraph.(experiment).stdv = zeros(nGroups,nWells);
    numNodes.(experiment) = zeros(nGroups,nWells);
    
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);

    for g=1:nGroups
        wellList = groups(g).WKS024;
        for w = 1:length(wellList)
            well = wellList{w};
            well_folder = fullfile(root, well);
            allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier);
            xc = allData.(experiment).(magnification).(well).xc;
            yc = allData.(experiment).(magnification).(well).yc;
            diameter = allData.(experiment).(magnification).(well).diameter;
            area = allData.(experiment).(magnification).(well).area;
            circularity = allData.(experiment).(magnification).(well).circularity;
            G = allData.(experiment).(magnification).(well).G;
            A = adjacency(G);
            xNodes = allData.(experiment).(magnification).(well).xNodes - xc; % set center of well to x=0
            yNodes = yc - allData.(experiment).(magnification).(well).yNodes; % set center of well to y=0

            N = numnodes(G);
            conComp = conncomp(G);
            sizeLargestGraph = sum( conComp == mode(conComp) );
            
            d = distances(G);
            d(d==inf) = [];
            d(d==0) = [];

            % topology measures
            avgPathLength.(experiment).mean(g,w) = mean(d);
            clusteringCoefficient.(experiment).mean(g,w) = sum( clustering_coef_bu(A) ) / N;
            assort.(experiment).mean(g,w) = assortativity(A,flag);
            largestGraph.(experiment).mean(g,w) = sizeLargestGraph / N;
            confluency.(experiment)(g,w) = 100* 4*sum(area) / (pi * diameter^2);
            numNodes.(experiment)(g,w) = N;
        end
    end
end


%% Plot result
numPlots = 4;

figure()
for e = 2:2
    experiment = experimentList{e};
    if strcmp(experiment, 'WKS024')
        network_specifier = '';
        wellList = {'B02', 'B03', 'B04', 'B05', 'B06', 'B07', ...
                    'C02', 'C03', 'C04', 'C05', 'C06', 'C07', ...
                    'D02', 'D03', 'D04', 'D05', 'D06', 'D07'};
        magnification = 'M10';
        scale = calculate_scale(magnification, fieldSize);
        color = colorHela;
                
    elseif strcmp(experiment, 'JJ005')
        network_specifier = '';
        wellList = {'B03', 'C02', 'D02', 'D03'};
        magnification = 'M20';
        scale = calculate_scale(magnification, fieldSize);
        color = colorCos;
    end
    
    
    maxAvgPathLentgh = 0;
    maxLogMeanN = 0;
    h = gobjects(1,2);
    
    for g=1:nGroups
        meanAvgPathLength = mean(avgPathLength.(experiment).mean(g,:));
        meanLargestGraph = mean(largestGraph.(experiment).mean(g,:));
        meanCC = mean(clusteringCoefficient.(experiment).mean(g,:));
        meanAssortativity = mean(assort.(experiment).mean(g,:));
        
        stdAvgPathLength = std(avgPathLength.(experiment).mean(g,:));
        stdLargestGraph = std(largestGraph.(experiment).mean(g,:));
        stdCC = std(clusteringCoefficient.(experiment).mean(g,:));
        stdAssortativity = std(assort.(experiment).mean(g,:));
        
        meanN = mean(numNodes.(experiment)(g,:));
    	meanConfluency = mean(confluency.(experiment)(g,:));
        
        allmeanN(g) = meanN;
        allmeanConfluency(g) = meanConfluency;
        
        stdN = std(numNodes.(experiment)(g,:));
    	stdConfluency = std(confluency.(experiment)(g,:));

        if meanAvgPathLength > maxAvgPathLentgh 
            maxAvgPathLentgh = meanAvgPathLength;
        end
        if log(meanN) > maxLogMeanN
            maxLogMeanN = log(meanN);
        end
        
        subplot(2,2,1)
        errorbar(meanConfluency, meanN, ...
                 stdN, stdN, ... 
                 stdConfluency, stdConfluency, ...
                 's','Color',color, 'MarkerFaceColor',color,'MarkerSize', 5);
        hold on
        
        subplot(2,2,2)
        errorbar(meanConfluency, meanLargestGraph, ...
                 stdLargestGraph, stdLargestGraph, ... 
                 stdConfluency, stdConfluency, ...
                 's','Color',color, 'MarkerFaceColor',color,'MarkerSize', 5);
        hold on
            
        subplot(2,2,3)
        h(1) = plot(log(numNodes.(experiment)(g,:)), avgPathLength.(experiment).mean(g,:),...
             '.','Color',color,'MarkerFaceColor',color, 'MarkerSize', 15);
        %errorbar(log(meanN), meanAvgPathLength, ...
        %         stdAvgPathLength, stdAvgPathLength, ... 
        %         log(stdN), log(stdN), ...
        %        's','Color',color, 'MarkerFaceColor',color);
        hold on

        subplot(2,2,4)
        errorbar(meanConfluency, meanCC, ...
                 stdCC, stdCC, ... 
                 stdConfluency, stdConfluency, ...
                 's','Color',color, 'MarkerFaceColor',color,'MarkerSize', 5);
        hold on
        
        
        
        %subplot(4,1,4)
        %errorbar(meanConfluency, meanAssortativity, ...
        %         stdAssortativity, stdAssortativity, ... 
        %         stdConfluency, stdConfluency, ...
        %         's','Color',color, 'MarkerFaceColor',color);
        %hold on
    end
    minX = 5;
    maxX = 9.5;
    
    minN = 100;
    maxN = 3e4;
    nAxis = [minN, maxN];
    
    minC = 2e-2;
    maxC = 1;
    cAxis = [minC,maxC];
    
    % fit linear model
    %[c1, c2, R2] = fit_linear_model(allmeanConfluency, allmeanN);
    
    subplot(2,2,1)
    %plot(cAxis, c1 + c2*cAxis, '-k','LineWidth', 1)
    xlabel('Confluency (%)')
	ylabel('Mean number of cells')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xticks([1,10,100])
    xticklabels({'1', '10', '100'})
    ylim(nAxis)
    %xlim([])
    
    subplot(2,2,2)
    xlabel('Confluency (%)')
	ylabel('Mean connectedness')
    set(gca, 'XScale', 'log')
    xticks([1,10,100])
    xticklabels({'1', '10', '100'})
    ylim([-0.05,1.1])
    
    subplot(2,2,3)
    h(2) = plot(minX:0.1:maxX, minX:0.1:maxX, '--r', 'LineWidth', 1);
    xlabel('ln(N)')
    ylabel('Average path length')
    xlim([minX,maxX])
    ylim([-1,65])
    lgnd = legend(h, {'HeLa', 'Random graph'}, 'Location', 'NorthWest');
    lgnd.ItemTokenSize = [7, 1];
    %set(gca, 'XScale', 'log')
    
    subplot(2,2,4)
    xlabel('Confluency (%)')
    ylabel('Mean clustering coefficient')
    set(gca, 'XScale', 'log')
    xticks([1,10,100])
    xticklabels({'1', '10', '100'})
    ylim([0,0.45]);
    
    
    
    
    %subplot(4,1,4)
    %xlabel('Confluency')
    %ylabel('Mean Assortativity')
    %set(gca, 'XScale', 'log')
end

%set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 5])
figName = fullfile('Figures/Topologies/',[experiment, '_', magnification, '_topologies.pdf']);
saveas(gcf, figName)  

function [c1, c2, R2] = fit_linear_model(x, y)

    mdl = fitlm(x,y);
    c1 = mdl.Coefficients.Estimate(1);
    c2 = mdl.Coefficients.Estimate(2);
    R2 = mdl.Rsquared.Ordinary;
end
