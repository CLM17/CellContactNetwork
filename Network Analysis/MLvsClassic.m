close all

%% Read data from the table
T = readtable('MLvsClassic.xlsx');
maxImg = 5;
wellList = {'B02'};

%% Initialize structure arrays
Q = struct();
Q.ml.FPNode = zeros(1,maxImg);
Q.ml.FNNode = zeros(1,maxImg);
Q.ml.NNodes = zeros(1,maxImg);
Q.ml.FPEdge = zeros(1,maxImg);
Q.ml.FNEdge = zeros(1,maxImg);
Q.ml.NEdges = zeros(1,maxImg);

Q.classic.FPNode = zeros(1,maxImg);
Q.classic.FNNode = zeros(1,maxImg);
Q.classic.NNodes = zeros(1,maxImg);
Q.classic.FPEdge = zeros(1,maxImg);
Q.classic.FNEdge = zeros(1,maxImg);
Q.classic.NEdges = zeros(1,maxImg);

%% Fill structure arrays with data from the table
for w = 1:length(wellList)
    well = wellList{w};
    wellT = T(strcmp(T.Well, well),:);
    
    for i=1:maxImg
        imgT = wellT(wellT.Img == i,:);
        
        classicT = imgT(strcmp(imgT.Method, 'classic'),:);
        mlT = imgT(strcmp(imgT.Method, 'ML'),:);
        
        Q.classic.FPNode(i) = classicT.FPNode;
        Q.classic.FNNode(i) = classicT.FNNode;
        Q.classic.NNodes(i) = classicT.NNodes;
        Q.classic.FPEdge(i) = classicT.FPEdge;
        Q.classic.FNEdge(i) = classicT.FNEdge;
        Q.classic.NEdges(i) = classicT.NEdges;
        
        Q.ml.FPNode(i) = mlT.FPNode;
        Q.ml.FNNode(i) = mlT.FNNode;
        Q.ml.NNodes(i) = mlT.NNodes;
        Q.ml.FPEdge(i) = mlT.FPEdge;
        Q.ml.FNEdge(i) = mlT.FNEdge;
        Q.ml.NEdges(i) = mlT.NEdges;

    end
end

dataClassic = [Q.classic.FPNode ./ Q.classic.NNodes; ...
               Q.classic.FNNode ./ Q.classic.NNodes; ...
               Q.classic.FPEdge ./ Q.classic.NEdges; ...
               Q.classic.FNEdge./ Q.classic.NEdges];
           
dataML = [Q.ml.FPNode ./ Q.ml.NNodes; ...
          Q.ml.FNNode ./ Q.ml.NNodes; ...
          Q.ml.FPEdge ./ Q.ml.NEdges; ...
          Q.ml.FNEdge./ Q.ml.NEdges];

allData = [mean(dataClassic,2), mean(dataML,2)] * 100;
allStd = [std(dataClassic,0,2), std(dataML,0,2)];

figure()
orange = [243,146,0] / 255;
blue = [30,144,255] / 255;
colors = {blue, orange};

b = bar(allData, 'grouped', 'FaceColor', 'flat', 'LineStyle', 'none');
for k = 1:size(allData,2)
    b(k).CData = colors{k};
end

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(allData, 1);
nbars = size(allData, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, allData(:,i), allStd(:,i), 'k', 'linestyle', 'none');
end

hold off

xticklabels({'False Positive Nodes', 'False Negative Nodes', ...
             'False Positive Edges', 'False Negative Edges'})
legend('Classic watershed', 'Learned watershed', 'Location', 'Northwest')
xtickangle(45)
ylabel('Percentage erroneous predictions (%)')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 5 2.8])
figName = fullfile('Figures', 'Fig_MLvsClassic.pdf');
saveas(gcf, figName)
