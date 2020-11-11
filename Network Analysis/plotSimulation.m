close all

experiment = 'WKS024';
magnification = '10x';
well = 'C02';
description = 'keepDistances';

%% ------------------------------START CODE--------------------------------
root = fullfile('..','Experiments', experiment, magnification);
well_folder = fullfile(root, well);
load(fullfile(well_folder, [well,'_simulatedGraph_keepDistances.mat']));

%% Calulcate centralities
Centralities = struct();
centralityNames = {'degree', 'betweenness', 'closeness', 'pagerank', 'eigenvector'};
nC = length(centralityNames);

for i = 1:nC
    cName = centralityNames{i};
    Centralities.(well).('G').(cName) = centrality(G, cName);
    Centralities.(well).('GSim').(cName) = centrality(GSim, cName);
end

%% Plot distributions
figure()
plot(smallDistances, pDist, '.')
hold on
plot(smallDistances, pEdgelength, '.')
plot(smallDistances, pConnect, '.')
hold off
legend('P(d_{ij})', 'P(d_{ij}|E_{ij})', 'P(E_{ij}|d_{ij})')
xlabel('Distance (pixels)')
ylabel('Probability')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 8 6])
figName = fullfile('Figures/Simulations/',[experiment, '_', magnification, '_', well,'_',description,'_probDistributions.png']);
saveas(gcf, figName)

%% Plot Centralities

figure()
for i = 1:nC
    cName = centralityNames{i};
    cObserved = Centralities.(well).('G').(cName);
    cSimulated = Centralities.(well).('GSim').(cName);
    
    subplot(3,nC,i)
    p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
    p1.NodeCData = cObserved;
    colormap jet
    title(cName)
    
    subplot(3,nC,i + nC)
    p2 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
    p2.NodeCData = cSimulated;
    colormap jet
    
    subplot(3,nC,i + 2*nC)
    histogram(cObserved)
    hold on
    histogram(cSimulated)
    hold off
end

subplot(3,nC,1)
ylabel('Observed')

subplot(3,nC,nC + 1)
ylabel('Simulated')

subplot(3,nC,2*nC + 1)
ylabel('Histogram')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 16 9])
figName = fullfile('Figures/Simulations/',[experiment, '_', magnification, '_', well,'_',description,'_centralities.png']);
saveas(gcf, figName)

%%
N = numnodes(G);
bco = Centralities.(well).('G').('betweenness') * 2 / ((N-1)*(N-2));
bcs = Centralities.(well).('GSim').('betweenness') * 2 / ((N-1)*(N-2));

figure()

subplot(2,2,1)
p1 = plot(G, 'XData', xNodes, 'YData', yNodes, 'markersize',2);
p1.NodeCData = bco;
colormap jet
colorbar
title('Observed')

subplot(2,2,2)
p1 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
p1.NodeCData = bcs;
colormap jet
colorbar
title('Simulated')

subplot(2,2,3)
histogram(bco)
title('Observed')

subplot(2,2,4)
histogram(bcs)
title('Simulated')
