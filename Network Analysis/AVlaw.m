close all

experiment = 'JJ005';
magnification = '20x';
well = 'C02';
fieldSize = 1104;
network_specifier = '';

%% ------------------------------START CODE--------------------------------

root = fullfile('..','..','Experiments', experiment, magnification);
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

if ~isfield(allData, well)
    scale = calculate_scale(magnification, fieldSize);
    allData = update_all_data(allData, well, well_folder, T, scale, network_specifier);
end
disp('All data loaded.')

% Get data of current well
measurementNames = allData.measurementNames;
G = allData.(well).G;
xNodes = allData.(well).xNodes;
yNodes = allData.(well).yNodes;

% well locations
xc = allData.(well).xc;
yc = allData.(well).yc;
diameter = allData.(well).diameter;
area = allData.(well).area;

%% Lewis' law: 
d = degree(G);

meanA = zeros(1, max(d)+1);
for i = 0:max(d)
    meanA(i+1) = mean(area(d == i));
end
figure()
plot(0:max(d), meanA / mean(area), '.')
hold on
plot(0:max(d), lewis_law(max(d)), '-r')

%% AV law:
figure()
mean_d = zeros(1, numnodes(G));
for n = 1:numnodes(G)
    nbh = neighbors(G, n);
    mean_d(n) = mean(degree(G, nbh));
end

% fit linear model m_d*d = c1 + c2*d
mdl = fitlm(d, mean_d .* d');
c1 = mdl.Coefficients.Estimate(1);
c2 = mdl.Coefficients.Estimate(2);
R2 = mdl.Rsquared.Ordinary;

figure()
plot(d, mean_d .* d', '.', d, c1 + c2*d, '-r');
ylabel('(degree) * (mean degree of neighbors)')
xlabel('degree')
title(['m_d*d = ', num2str(c1), ' + ', num2str(c2), ' * d, R2 = ',num2str(R2)])

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 5.5 4])
figName = fullfile('Figures/',[experiment, '_', magnification, '_', well,'AVLaw.png']);
saveas(gcf, figName)


% figure()
% plot(d, mean_d, '.')

%%
% init
figure()
m_d = cell(1, max(d) + 1);
for i = 0:max(d)
    m_d{i+1} = [];
end

for n = 1:numnodes(G)
	deg = degree(G, n);
    nbh = neighbors(G, n);
    m_d{deg+1} = [m_d{deg+1}; degree(G, nbh)];
end

mean_m = zeros(1, max(d) + 1);
for i = 0:max(d)
    mean_m(i+1) = mean(m_d{i+1});
end

dAxis = 0:max(d);
plot(dAxis, mean_m, '*')
hold on
plot(dAxis, 5 + 8./dAxis, '-r')
legend('Data', '5 + 8/n')
xlabel('n')
ylabel('m')

mdl = fitlm(dAxis, mean_m .* dAxis);
c1 = mdl.Coefficients.Estimate(1);
c2 = mdl.Coefficients.Estimate(2);
R2 = mdl.Rsquared.Ordinary;

%% confluency
confluency = 4*sum(area) / (pi * diameter^2);

figure()
plot(dAxis, mean_m .* dAxis, '*')
xlabel('d')
ylabel('m*d')

%%
function meanA = lewis_law(maxD)
    degrees = 0:maxD;
    meanA = (degrees - 2) / 4;
    %meanA = (degrees / 6).^2;
end
