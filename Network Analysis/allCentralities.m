magnification = 'M20';
wellRows = 2:7;
wellColumns = {'B','C', 'D'};
allData = struct;
experiment = 'WKS024';
fieldSize = 1104;
network_specifier = '';


root = fullfile('..','..','Experiments', experiment, magnification);
xlsfileName = fullfile(root, 'Well locations.xlsx');
T = readtable(xlsfileName);

scale = calculate_scale(magnification, fieldSize);

%% Load data
for r = 1:length(wellRows)
    row = [num2str(0), num2str(wellRows(r))];
    if r == 1
        magnification = 'M20';
        network_specifier = '_ml';
    else
        magnification = 'M10';
        network_specifier = '';
    end
    
    scale = calculate_scale(magnification, fieldSize);
    root = fullfile('..','..','Experiments', experiment, magnification);
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
    
    for c = 1:length(wellColumns)
        column = wellColumns{c};
        well = [column, row];
        well_folder = fullfile(root, well);
 
        allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier);
    end
end

%% Calculate centrality measures

cMeasures = {'degree', 'closeness', 'betweenness'};
minValues.degree = inf;
minValues.closeness = inf;
minValues.betweenness = inf;

maxValues.degree = 0;
maxValues.closeness = 0;
maxValues.betweenness = 0;

numNodes = zeros(length(wellRows), length(wellColumns));

for m = 1:length(cMeasures)
    cName = cMeasures{m};
    disp(cName)
    for r = 1:length(wellRows)
        if r == 1
            magnification = 'M20';
            network_specifier = '_ml';
        else
            magnification = 'M10';
            network_specifier = '';
        end
        row = [num2str(0), num2str(wellRows(r))];
        for c = 1:length(wellColumns)
            column = wellColumns{c};
            well = [column, row];

            G = allData.(experiment).(magnification).(well).G;
            numNodes(r, c) = numnodes(G);
            
            cMeasure = calculate_normalized_centrality(G, cName);
            allData.(experiment).(magnification).(well).(cName) = cMeasure;
            
            if min(cMeasure) < minValues.(cName)
                minValues.(cName) = min(cMeasure);
            end
            if max(cMeasure) > maxValues.(cName)
                maxValues.(cName) = max(cMeasure);
            end
            
        end
    end
end
        
 %% Plot
 close all
 
for m = 1:length(cMeasures)
    count = 0;
    cName = cMeasures{m};
    figure(m)
    for r = 1:length(wellRows)
        if r == 1
            magnification = 'M20';
            network_specifier = '_ml';
        else
            magnification = 'M10';
            network_specifier = '';
        end
        row = [num2str(0), num2str(wellRows(r))];
        for c = 1:length(wellColumns)
            column = wellColumns{c};
            well = [column, row];

            G = allData.(experiment).(magnification).(well).G;
            xc = allData.(experiment).(magnification).(well).xc;
            yc = allData.(experiment).(magnification).(well).yc;
            xNodes = allData.(experiment).(magnification).(well).xNodes - xc;
            yNodes = yc - allData.(experiment).(magnification).(well).yNodes;

            cMeasure = allData.(experiment).(magnification).(well).(cName);

            count = count + 1;
            subplot(6,3,count)
            p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',1);
            p1.NodeCData = cMeasure;
            
            caxis manual
            caxis([minValues.(cName), maxValues.(cName)]);
            colormap jet
            %colorbar;
%             if c==3
%                 colorbar
%             end
            %xticks([-4000,0,4000]);
            %yticks([-4000,0,4000]);
            %xticklabels({-4,0,4});
            xticklabels('');
            yticklabels('');

        end
    end
    
    hp16 = get(subplot(6,3,16),'Position');
    hp18 = get(subplot(6,3,18),'Position');
    
    w = hp18(1) + hp18(3) - hp16(1);
    h = 0.02;
    pos = [hp16(1), hp16(2)-0.05, w, h];
    
    colorbar('Location', 'southoutside','Position', pos)
    
    %set(gcf,'PaperOrientation','landscape');
    set(gcf,'Color','w','Units','inches','Position',[1 1 6 11])
    figName = fullfile('Figures/Centralities',[experiment,'_all_',cName,'.png']);
    saveas(gcf, figName)
end