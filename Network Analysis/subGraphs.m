experiment = 'WKS024';
magnification = '10x';
well = 'B02';
unit = 'pixels';
fieldSize = 1104;           % Size of 1 field.

loadImage = true;           % Set to true if you want to display an image
N = 3;                      % Number of nodes in subgraph you want to display on the image.

nodeSize = 3;               % node size
nodeColor = 'w';            % color of nodes ('w'=white, 'k'=black, 'g'=green, etc)
lineWidth = 2;              % thickness of edges (= lines)
edgeColor = 'w';            % color of edges ('w'=white, 'k'=black, 'g'=green, etc)
edgeTransparency = 0.7;     % transparency of edges (0=fully transparent, 1=not transparent)

%% ------------------------------START CODE--------------------------------

outputFolder = fullfile('SubgraphAnalysis', experiment, magnification, well);
if ~isfolder(outputFolder)
    mkdir(outputFolder);
    disp('Created new output folder for this well.')
end

root = fullfile('Experiments', experiment, magnification);
well_folder = fullfile(root, well);

% Load image and graph if this wasn't done already
if ~exist('T','var')
    well_folder = fullfile(root, well);
    xlsfileName = fullfile(root, 'Well locations.xlsx');
    T = readtable(xlsfileName);
end


if ~exist('allData','var')
    allData = struct;
end

scale = 1; % do not convert sizes to um.
if ~isfield(allData, well)
    allData = update_all_data(allData, well, well_folder, T, scale);
    
    if (loadImage)
        disp(['Loading the image of well ',well,'...'])
        [fused, cmap] = imread(fullfile(well_folder,[well,'_fused_RGB.tif']));
        allData.(well).fused = fused;
        allData.(well).cmap = cmap;
    end
end
disp('All data loaded.')

% Get data of current well
G = allData.(well).G;
xNodes = allData.(well).xNodes;
yNodes = allData.(well).yNodes;

% well locations
xc = allData.(well).xc;
yc = allData.(well).yc;
diameter = allData.(well).diameter;

% image
if loadImage
    fused = allData.(well).fused;
    cmap = allData.(well).cmap;
end

%% 
subgraphs = conncomp(G);
nSubGraphs = length(subgraphs);
count = zeros(1, length(subgraphs));    

% count(i) number of nodes in subgraph to which node i belongs.
for i = 1:nSubGraphs
    count(i) = sum(subgraphs == subgraphs(i));
end

% S(n) is the amount of subgraphs with size n.
subgraphSize = cell(1, length(unique(count)));
S = zeros(1, length(unique(count)));
c = 0;
for n = unique(count)
    c = c + 1;
    subgraphSize{c} = num2str(n);
    S(c) = sum(count == n) / n;
end

figure()

subplot(2,3,1:3)
xaxis = 1:length(unique(count));
bar(xaxis, S);
xticks(xaxis)
xticklabels(subgraphSize);
ylabel('Number of subgraphs')
xlabel('Size of subraph')

c = 0;
for n = 3:5
    c = c + 1;
    subplot(2,3,3+c)
    GSub = subgraph(G, count == n);
    plot(GSub, 'NodeLabel',{})
    title(['n = ', num2str(n)])
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 10 6])
figName = fullfile(outputFolder,[experiment, '_', magnification '_subGraphs.png']);
%saveas(gcf, figName)

%% Plot
% Draw graph as network & save the file

% initialise
families = struct();
nArray = [3,4,5];
for i = 1:length(nArray)
    families(i).complete = 0;
    families(i).cyclic = 0;
    families(i).line = 0;
    families(i).star = 0;
    families(i).other = 0;
end

for i = 1:length(nArray)
    n = nArray(i);
    GSub = subgraph(G, count == n);
    subsubgraphs = conncomp(GSub);
    for j = 1:max(subsubgraphs)
        GSubSub = subgraph(GSub, subsubgraphs == j);
        
        complete = iscomplete(GSubSub);
        cyclic = iscyclic(GSubSub);
        line = isline(GSubSub);
        star = isstar(GSubSub);
        
        families(i).complete = families(i).complete + complete;
        families(i).cyclic = families(i).cyclic + cyclic;
        families(i).line = families(i).line + line;
        families(i).star = families(i).star + star;
        
        if (~complete && ~cyclic && ~line && ~star)
            families(i).other = families(i).other + 1;
        end
    end
end
%% Plot subgraph
if loadImage

    GSub = subgraph(G, count == N);

    figure()
    imshow(fused,cmap)

    hold on

    plot(GSub, 'XData', xNodes(count == N), 'YData', yNodes(count == N),...
        'MarkerSize',nodeSize,...
        'NodeLabel',{},...
        'NodeColor',nodeColor,...
        'LineWidth',lineWidth,...
        'EdgeColor',edgeColor,...
        'EdgeAlpha',edgeTransparency);
end

%% ------------------------------FUNCTIONS---------------------------------

function b = iscomplete(G)
    % number of edges is N(N-1)/2
    N = numnodes(G);
    E = numedges(G);
    b = (E == N*(N-1)/2);
end

function b = isline(G)
    % 2 nodes have degree 1, N-2 nodes have degree 2.
    N = numnodes(G);
    d = degree(G);
    b = false;
    if (sum(d==1) == 2) && (sum(d==2) == N-2)
        b = true;
    end
end

function b = iscyclic(G)
    % All nodes have degree 2.
    N = numnodes(G);
    d = degree(G);
    b = sum(d==2) == N;
end

function b = isstar(G)
    % The degree of one node equals the number of edges. 
    b = false;
    N = numnodes(G);
    if N <= 3
        return
    end
    E = numedges(G);
    d = degree(G);
    for i = 1:N
        if d(i) == E
            b = true;
            return
        end
    end
end
