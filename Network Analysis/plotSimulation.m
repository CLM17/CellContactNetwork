close all

experiment = 'WKS024';
magnification = '10x';
well = 'B02';
description = 'simulateDistances';

%% ------------------------------START CODE--------------------------------
root = fullfile('..','..','Experiments', experiment, magnification);
well_folder = fullfile(root, well);
load(fullfile(well_folder, [well,'_simulatedGraph_',description,'.mat']));

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
