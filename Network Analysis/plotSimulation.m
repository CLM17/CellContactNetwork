close all
addpath('BrainConnectivity')

experiment = 'WKS024';
magnification = 'M20';
well = 'D02';
description = '20201230';
N = 2; % number of simulations

%% ------------------------------START CODE--------------------------------

root = fullfile('..','..','Experiments', experiment, magnification);
well_folder = fullfile(root, well);

% Get simulated data
simData = struct();
for n = 1:N
    nr = num2str(n-1);
    data = load(fullfile(well_folder, [well,'_simulatedGraph_',description,'_',nr,'.mat']));
    simData(n).GSim = data.GSim;
    simData(n).xSim = data.xSim;
    simData(n).ySim = data.ySim;
end

% Get observed graph
obsData = struct;
obsData.G = data.G;
obsData.xNodes = data.xNodes;
obsData.yNodes = data.yNodes;

%% Calulcate centralities
centralityNames = {'degree', 'closeness', 'betweenness'};
nC = length(centralityNames);

for i = 1:nC
    cName = centralityNames{i};
    % observed data
    obsData.(cName) = calculate_normalized_centrality(obsData.G, cName);
    % simulated data
    for n = 1:N
        GSim = simData(n).GSim;
        simData(n).(cName) = calculate_normalized_centrality(GSim, cName);
    end
end

%% Collect centralities of all simulations in 3 arrays
simDegree = [];
simCloseness = [];
simBetweenness = [];

for n = 1:N
    simDegree = [simDegree; simData(n).degree];
    simBetweenness = [simBetweenness; simData(n).betweenness];
    simCloseness = [simCloseness; simData(n).closeness];
end

%% Assortativity, connectedness & avg path length

% Observed

% Connectedness
conComp = conncomp(obsData.G);
sizeLargestGraph = sum( conComp == mode(conComp) );
% Avg path length
d = distances(obsData.G);
d(d==inf) = [];
d(d==0) = [];
% Assortativity
flag = 0;

obsData.assort = full(assortativity(adjacency(obsData.G), flag));
obsData.connectedness = sizeLargestGraph / numnodes(obsData.G);
obsData.avgPathLength = mean(d);

% simulated

for i = 1:N

    GSim = simData(n).GSim;
    % Connectedness
    conComp = conncomp(GSim);
    sizeLargestGraph = sum( conComp == mode(conComp) );
    % Avg path length
    d = distances(GSim);
    d(d==inf) = [];
    d(d==0) = [];
    % Assortativity
    flag = 0;

    simData(n).assort = full(assortativity(adjacency(GSim), flag));
    simData(n).connectedness = sizeLargestGraph / numnodes(GSim);
    simData(n).avgPathLength = mean(d);
end



%% Compare graph measures
figure()

 GSim = simData(1).GSim;
 xSim = simData(1).xSim;
 ySim = simData(1).ySim;
for i = 1:nC 
    cName = centralityNames{i};
   
    subplot(2,3,i)
    p1 = plot(GSim, 'XData', xSim, 'YData', ySim, 'markersize',2);
    p1.NodeCData = simData(1).(cName);
    h = colorbar('location','NorthOutside');
    h.Label.String = ['Simulated ',cName];

    colormap jet
    xticks([-4e3,0,4e3])
    yticks([-4e3,0,4e3])
    xticklabels({'-4','0','4'})
    yticklabels({'-4','0','4'})
    xlabel('x (mm)')
    ylabel('y (mm)')
end

subplot(2,3,4)
plotDegree.observed = obsData.degree;
plotDegree.simulated = simDegree;
v = violinplot(plotDegree,1);
for i = 1:2
    v(i).EdgeColor = [1,1,1];
    v(i).ShowData = false;
end
ylabel('Degree centrality')

subplot(2,3,5)
plotCloseness.observed = obsData.closeness;
plotCloseness.simulated = simCloseness;
v = violinplot(plotCloseness, 0.0007);
for i = 1:2
    v(i).EdgeColor = [1,1,1];
    v(i).ShowData = false;
end
ylabel('Closeness centrality')

subplot(2,3,6)
plotBetweenness.observed = obsData.betweenness;
plotBetweenness.simulated = simBetweenness;
v = violinplot(plotBetweenness, 0.001);
for i = 1:2
    v(i).EdgeColor = [1,1,1];
    v(i).ShowData = false;
end
ylim([0,0.05])
ylabel('Betweenness centrality')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 10 7.5])
figName = fullfile('Figures/Simulations/',[experiment, '_', well, '_violinplots.png']);
saveas(gcf, figName)

%% Fit power law on simulated betweenness centrality

nBins = 60;
bcSim = simData(1).betweenness;

% Fit on betweenness centrality
[hc, binEdges] = histcounts(bcSim, 'NumBins', nBins, 'Normalization', 'probability');
binCenters = (binEdges(2:end) + binEdges(1:end-1) ) / 2;

F = @(a,xdata)a(1)*xdata.^(-a(2));
a0 = [1,1];
[a,resnorm,~,exitflag,output] = lsqcurvefit(F,a0,binCenters,hc);

xmin = 0;
xmax = max(bcSim);
xaxis = linspace(xmin, xmax, 1000);
% pd = fitdist(bc, 'exponential');
% y = pdf(pd, xaxis);
figure()
histogram(bcSim, nBins,'Normalization','probability', 'LineStyle', 'none')
hold on
plot(xaxis, F(a, xaxis), '-k', 'LineWidth', 1)
set(gca, 'YScale', 'log','XScale', 'log')

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
    cObserved = simData.(well).('G').(cName);
    cSimulated = simData.(well).('GSim').(cName);
    
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
bco = simData.(well).('G').('betweenness') * 2 / ((N-1)*(N-2));
bcs = simData.(well).('GSim').('betweenness') * 2 / ((N-1)*(N-2));

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

subplot(2,2,3:4)
histogram(bco)
hold on
histogram(bcs)
legend('Observed', 'Simulated')

%% Subgraphs
subgraphs = conncomp(G);
%%
for s = unique(subgraphs)
    GSub = subgraph(G, subgraphs == s);
end
%%
subgraphs = conncomp(GSim);
rgb = vals2colormap(subgraphs, 'hsv');
figure()

subplot(4,5,[1:4, 6:9, 11:14, 16:19])

plot(GSim, 'XData', xSim, 'YData', ySim, 'NodeColor', rgb, 'markersize',3);
colormap hsv

c = 1;
while c < 5
    nodes = (subgraphs == c);
    if sum(nodes)==0
        continue
    end
    subplot(4,5,5*c)
    GSub = subgraph(GSim, nodes);
    N = numnodes(GSub);
    colour = rgb(subgraphs==c,:);
    histogram( centrality(GSub, 'betweenness') / ((N-1)*(N-2)), 'FaceColor', colour(1,:))
    c = c + 1;
end

%% Check if distances are in the same distribution
[X,Y] = meshgrid(xNodes, yNodes);
A = adjacency(G);
distE = sqrt( (X-X').^2 + (Y-Y').^2 ) .* A;

[Xsim, Ysim] = meshgrid(xSim, ySim);
Asim = adjacency(GSim);
distSimE = sqrt( (Xsim-Xsim').^2 + (Ysim-Ysim').^2 ) .* Asim;

%%
[h,p] = kstest2(full(distE(distE>0)), full(distSimE(distSimE>0)));

%%
xSim = data.xSim;
ySim = data.ySim;
xNodes = data.xNodes;
yNodes = data.yNodes;

rObs = sqrt(xNodes.^2 + yNodes.^2);
rSim = sqrt(xSim.^2 + ySim.^2);

[h,p] = kstest2(rObs, rSim);
%%
nNodes = numnodes(data.G);
pdf_r = fitdist(rObs,'kernel','Width',5);
pdf_r = truncate(pdf_r, 0, max(rObs));

rSim = random(pdf_r, nNodes, 1);
[h,p] = kstest2(rObs, rSim);
disp(h)

figure()
histogram(rObs)
hold on
histogram(rSim)
%% Plot subgraphs

%% 
% subgraphs = conncomp(GSim);
% nSubGraphs = length(subgraphs);
% count = zeros(1, length(subgraphs));    
% 
% % count(i) number of nodes in subgraph to which node i belongs.
% for i = 1:nSubGraphs
%     count(i) = sum(subgraphs == subgraphs(i));
% end
% 
% % S(n) is the amount of subgraphs with size n.
% subgraphSize = cell(1, length(unique(count)));
% S = zeros(1, length(unique(count)));
% c = 0;
% for n = unique(count)
%     c = c + 1;
%     subgraphSize{c} = num2str(n);
%     S(c) = sum(count == n) / n;
% end
% 
% c = 0;
% for n = 3:5
%     c = c + 1;
%     subplot(2,3,3+c)
%     GSub = subgraph(GSim, count == n);
%     plot(GSub, 'NodeLabel',{})
%     title(['n = ', num2str(n)])
% end



%%
function rgb = vals2colormap(vals, colormap, crange)

% Take in a vector of N values and return and return a Nx3 matrix of RGB

% values associated with a given colormap

%

% rgb = AFQ_vals2colormap(vals, [colormap = 'jet'], [crange])

%

% Inputs:

% vals     = A vector of values to map to a colormap or a cell array of

%            vectors of values

% colormap = A matlab colormap. Examples: colormap = 'autumn';

%            colormap = 'jet'; colormap = 'hot';

% crange   = The values to map to the minimum and maximum of the colormap.

%            Defualts to the full range of values in vals.

%

% Outputs:

% rgb      = Nx3 matrix of rgb values mapping each value in vals to the

%            corresponding rgb colors.  If vals is a cell array then rgb

%            will be a cell array of the same length

%

% Example:

% vals = rand(1,100);

% rgb = AFQ_vals2colormap(vals, 'hot');

%

% Copyright Jason D. Yeatman, June 2012



if ~exist('colormap','var') || isempty(colormap)

    colormap = 'jet';

end



%

if ~iscell(vals)

    if ~exist('crange','var') || isempty(crange)

        crange = [min(vals) max(vals)];

    end

    % Generate the colormap

    cmap = eval([colormap '(256)']);

    % Normalize the values to be between 1 and 256

    vals(vals < crange(1)) = crange(1);

    vals(vals > crange(2)) = crange(2);

    valsN = round(((vals - crange(1)) ./ diff(crange)) .* 255)+1;

    % Convert any nans to ones

    valsN(isnan(valsN)) = 1;

    % Convert the normalized values to the RGB values of the colormap

    rgb = cmap(valsN, :);

elseif iscell(vals)

    if ~exist('crange','var') || isempty(crange)

        crange = [min(vertcat(vals{:})) max(vertcat(vals{:}))];

    end

    % Generate the colormap

    cmap = eval([colormap '(256)']);

    for ii = 1:length(vals)

        % Normalize the values to be between 1 and 256 for cell ii

        valsN = vals{ii};

        valsN(valsN < crange(1)) = crange(1);

        valsN(valsN > crange(2)) = crange(2);

        valsN = round(((valsN - crange(1)) ./ diff(crange)) .* 255)+1;

        % Convert any nans to ones

        valsN(isnan(valsN)) = 1;

        % Convert the normalized values to the RGB values of the colormap

        rgb{ii} = cmap(valsN, :);

    end

end

return
end
