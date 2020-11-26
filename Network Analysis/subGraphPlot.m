close all

experiment = 'WKS024';
magnification = '10x';
well = 'B02';
fieldSize = 1104;           % Size of 1 field.

wellArray = {'B02', 'B03', 'B04', 'B05', ...
             'C02', 'C03', 'C04', 'C05', ...
             'D02', 'D03', 'D04', 'D05'};
numRows = 3;
numCols = 4;

%% ------------------------------START CODE--------------------------------

root = fullfile('..','..','Experiments', experiment, '10x old', magnification);

% Load image and graph if this wasn't done already
if ~exist('T','var')
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
end


if ~exist('allData','var')
    allData = struct;
end

scale = 1; % do not convert sizes to um.
numCell = zeros(numRows, numCols);
numCyclic = zeros(numRows, numCols);
numLine = zeros(numRows, numCols);
p = zeros(1, numCols);

count = 0;
for r = 1:numRows
    for c = 1:numCols
        count = count + 1;
        well = wellArray{count};
        disp(well)
        wellFolder = fullfile(root, well);
        if ~isfield(allData, well)
            allData = update_all_data(allData, well, wellFolder, T, scale);
        end

        % Get data of current well
        G = allData.(well).G;
        
        subGraphFile = fullfile(wellFolder, 'subGraphs_checked.csv');
        table = readtable(subGraphFile);
        N3_values = table.total;
        
        numCell(r,c) = numnodes(G);
        numCyclic(r,c) = N3_values(1);
        numLine(r,c) = N3_values(3);
        
    end
end

for c = 1:numCols
    p(c) = anova1([numCyclic(:,c), numLine(:,c)], {'Cyclic','Line'}, 'off');
end

%%
close all
figure()
errorbar(mean(numCell), mean(numCyclic), std(numCyclic), std(numCyclic),...
    std(numCell), std(numCell), 's','Color','b', 'MarkerFaceColor','b');

hold on

errorbar(mean(numCell), mean(numLine), std(numLine), std(numLine),...
    std(numCell), std(numCell), 's','Color','r', 'MarkerFaceColor','r');

legend('Cyclic', 'Line')
xlabel('Number of cells')
ylabel('Number of subgraphs')
set(gca, 'XSCale', 'log')

for c = 1:numCols
    x = mean(numCell(:,c));
    y = max( mean(numLine(:,c)) + std(numLine(:,c)), mean(numCyclic(:,c)) + std(numCyclic(:,c))) + 1;
    text(x,y,['p = ', num2str(p(c),2)], 'HorizontalAlignment', 'center');
end

xlim([-100,11000])
set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 6 4])
figName = fullfile('SubgraphAnalysis',[experiment, '_', magnification, '_subgraphResults.png']);
saveas(gcf, figName) 
