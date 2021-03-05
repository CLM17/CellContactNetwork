clear all
close all

disp('Choose the fused_RGB image you want to show')
[FusedImage,RawPath] = uigetfile('.tif'); %Choose the fused image you want
cd(RawPath)

disp('Choose the corresponding fused_segmentation image')
SegmentedImage = uigetfile('.tif'); %Choose the segmented image you want

fused = imread(FusedImage);
[seg, cmap] = imread(SegmentedImage);

rgb = label2rgb(seg,'jet',[0,0,0],'shuffle');

% Check the fields of 'Network' datastructure
% to see all cell measurements
Network = load('Network','-mat');
G = graph(Network.contact_matrix);

%% Physical parameters
Origin = round([length(fused) length(fused)]/2);
xNodes = Network.centroid1 - Origin(1);
yNodes = Network.centroid0 - Origin(2);
pixtomm = 6.3/length(fused); 
pixtoum = 6300/length(fused);

%% Define patch of interest
xpatch = 4000;
ypatch = 4000;
patchsize = 1104;

%%

%Show fused image
figure(1)
imshow(fused)

%Show segmented fused image
figure(2)
imshow(rgb)
imwrite(rgb,'SegmentedImage_rgb.tif')

%Show patch of interest 
figure()
fused_patch = fused(ypatch:(ypatch+patchsize),xpatch:(xpatch+patchsize),:);
imshow(fused_patch)
imwrite(fused_patch,'fused_patch.tif')

%Show segmented patch of interest
rgb_patch = label2rgb(seg(ypatch:(ypatch+patchsize),xpatch:(xpatch+patchsize)),'jet',[0,0,0],'shuffle');
figure()
imshow(rgb_patch)
imwrite(rgb_patch,'SegmentedImage_rgb_patch.tif')

%Show network on patch of interest
figure()
imshow(fused) 
hold on
plot(G,'XData',Network.centroid1,'YData',Network.centroid0,'MarkerSize',1,'LineWidth',1,'NodeColor','w','EdgeColor','w')
hold on
%rectangle('Position',[ypatch xpatch patchsize patchsize])
%axis([ypatch ypatch+patchsize xpatch xpatch+patchsize])
axis equal

saveas(gcf,'NetworkOverlay.png')

%%
figure()
plot(G,'XData',Network.centroid1,'YData',Network.centroid0,'MarkerSize',1,'LineWidth',1,'NodeColor','k','EdgeColor','k')
set(gcf,'Color','w')
box off
axis off
axis equal
axis([0 length(fused)/3 0 length(fused)/7])
saveas(gcf,'MCNGraph.png')


%%
NumNodes = numnodes(G);
Confluency = sum(seg>0,'all')/(pi*(length(fused)/2)^2);
disp(strcat("There are ",num2str(NumNodes)," cells in this well, ",'and the confluency is'," ",num2str(Confluency*100,3),'%'))

% Calculate distributions of distances
[X,Y] = meshgrid(xNodes, yNodes);
distances = sqrt( (X-X').^2 + (Y-Y').^2)*pixtoum;
edgeLength = distances .* adjacency(G);
existingEdgeLength = edgeLength(edgeLength > 0);

figure()
histogram(existingEdgeLength,'Normalization','probability','FaceColor','k')
xlabel('Edge length (\mum)')
ylabel('Probability')

% Calculate the cell density as a function of radius
nBins = 20;
WellDiameter = 6.3; %mm

Area = pi * WellDiameter^2 / 4;
r = sqrt( xNodes.^2 + yNodes.^2 )*pixtomm;
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
    thBinArea(i) = ( (thBinEdges(i+1) - thBinEdges(i)) / (2*pi) ) * Area;
end

rBinCenters = (rBinEdges(2:end) + rBinEdges(1:end-1)) / 2;
thBinCenters = (thBinEdges(2:end) + thBinEdges(1:end-1)) / 2;
rDensity = rBinned ./ rBinArea;  
thDensity = thBinned ./ thBinArea;

% Plot biophysical properties (R, theta)
figure()

bar(rBinCenters, rDensity,'FaceColor','k');
xlabel('Radial distance to well center (mm)')
ylabel('Density (cells / mm^2)')
ylim([0 1250])
xticks(0:0.5:3.5)
yticks(0:250:1250)

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[9 1 3.5 2])
saveas(gcf,'rDistributionOfCells.png')

figure()
bar(thBinCenters,thDensity,'FaceColor','k');
xlabel('Angular position (radians)')
ylabel('Density (cells/mm^2)')
ylim([0 1250])
xticks(-4:1:4)
yticks(0:250:1250)

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[9 1 3.5 2])
saveas(gcf,'thDistributionOfCells.png')

%% Graph theoretical analysis
figure()
plot(G,'MarkerSize',1,'LineWidth',1,'NodeColor','k','EdgeColor','k')
saveas(gcf,'Subgraphs.png')

Centralities = struct;
Centralities.Degree = centrality(G,'degree');
Centralities.Closeness = centrality(G,'closeness')*(NumNodes-1);
Centralities.Betweenness = 2*centrality(G,'betweenness')/((NumNodes-1)*(NumNodes-2)); %Normalized
Centralities.Pagerank = centrality(G,'pagerank');

figure()
subplot(2,2,1)
histogram(Centralities.Degree)
xlabel('Degree')
ylabel('Number of cells')
subplot(2,2,2)
histogram(Centralities.Closeness)
xlabel('Normalized closeness centrality')
ylabel('Number of cells')
subplot(2,2,3)
histogram(Centralities.Betweenness)
set(gca,'YScale','log')
xlabel('Normalized betweenness centrality')
ylabel('Number of cells')
subplot(2,2,4)
histogram(Centralities.Pagerank)
xlabel('PageRank')
ylabel('Number of cells')

%set(gcf, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[9 1 10 8])
saveas(gcf,'CentralityDistributions.png')

%% Plot the network with a color determined by a quantity (ColorQuant)
figure()
imshow(fused) 
hold on

ColorQuant = 'Closeness';
Color = vals2colormap(Centralities.(ColorQuant),'jet',[min(Centralities.(ColorQuant)) max(Centralities.(ColorQuant))]);

plot(G,'XData',object_com(:,2),'YData',object_com(:,1),'MarkerSize',4,'LineWidth',2,'NodeColor',Color,'EdgeColor','w')
axis equal
%axis([ypatch ypatch+patchsize xpatch xpatch+patchsize])
set(gcf,'Color','k')

caxis manual
caxis([min(Centralities.(ColorQuant)) max(Centralities.(ColorQuant))]);
colormap jet
colorbar('Location', 'southoutside','Color','w')
xlabel(strcat((ColorQuant)," centrality"),'Color','w')

saveas(gcf,strcat(ColorQuant,'NetworkOverlay.png'))



%% Betweenness paths
figure()
h = histogram(Centralities.betweenness,'Normalization','probability','FaceColor','k','FaceAlpha',1);
hold on

% xdata = (h.BinEdges(1:end-1)+diff(h.BinEdges))/2;
% ydata = h.Values;
% 
% fun = @(param,xdata)param(1)*xdata.^(-param(2));
% param0 = [1 1];
% lsqcurvefit(fun,param0,xdata,ydata)
% plot(xdata,fun(ydata,xdata),'b-')

set(gca,'YScale','log')
yticks(10.^(-5:1:0))
xlabel('Betweenness centrality')
ylabel('Probability')

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[9 1 3.5 2])
saveas(gcf,'BetweennessHistogram.png')

%% Functions

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