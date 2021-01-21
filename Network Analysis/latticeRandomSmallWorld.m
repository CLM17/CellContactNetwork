close all
clear all

% Plotting specifications
nodeSize = 7;               % node size
nodeColor = 'k';            % color of nodes ('w'=white, 'k'=black, 'g'=green, etc)
lineWidth = 1.5;              % thickness of edges (= lines)
edgeColor = 'k';            % color of edges ('w'=white, 'k'=black, 'g'=green, etc)
edgeTransparency = 1;       % transparency of edges (0=fully transparent, 1=not transparent)

% Graph specifications
%xNodes = [0,1.2,2.1,3.5,3,5.1,4.9]';
%yNodes = [0,2.1,0.8,4.1,5.8,4,6.2]';
% edges = {[1,2], ...
%          [1,3], ...
%          [2,3], ...
%          [3,4], ...
%          [4,5], ...
%          [4,6], ...
%          [5,7], ...
%          [6,7]};
% numNodes = length(xNodes);

xNodes = [0:2:8, ...
          1:2:9, ...
          0:2:8, ...
          1:2:9, ...
          0:2:8];
yNodes = [zeros(1,5), ...
          ones(1,5), ...
          ones(1,5)+1, ...
          ones(1,5)+2, ...
          ones(1,5)+3] * sqrt(3);
      
numNodes = length(xNodes);
          
[X,Y] = meshgrid(xNodes, yNodes);
d = sqrt( (X-X').^2 + (Y-Y').^2 );
A = (d <= 2.1);
A(d==0) = 0;

% Make graph
G = graph(A);
% for i = 1:length(edges)
%     edge = edges{i};
%     node1 = edge(1);
%     node2 = edge(2);
%     G = addedge(G, node1, node2);
% end

%%
addpath('BrainConnectivity')
A = adjacency(G);

alphaArray = [0, 0.15, 1];

figure()
for a = 1:length(alphaArray)
    
    % Rewire edges
    alpha = alphaArray(a);
    R = randomizer_bin_und(A,alpha);
    G_rand = graph(R);
    cc = clustering_coef_bu(R);
    disp(mean(distances(G_rand), 'all'));
    subplot(1,3,a)

    % Plot graph
    p = plot(G_rand, 'XData', xNodes, 'YData', yNodes, ...
            'NodeCData', centrality(G_rand, 'betweenness'),...
            'MarkerSize',nodeSize,...
            'NodeLabel',{},...
            'LineWidth',lineWidth,...
            'EdgeColor',edgeColor,...
            'EdgeAlpha',edgeTransparency);
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 13 3])
figName = fullfile('Figures/','latticeRandomSmallWorld.pdf');
saveas(gcf, figName);


    