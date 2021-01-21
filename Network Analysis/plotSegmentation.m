experiment = 'WKS024';
magnification = 'M20';
well = 'D02';               % well name
network_specifier = '_classic';

root = fullfile('..','..','Experiments',experiment,magnification);
well_folder = fullfile(root, well);
xlsfileName = fullfile(root, 'Well locations.xlsx');
T = readtable(xlsfileName);
scale = 1;

%allData = update_all_data(allData, experiment, magnification,...
%                              well, well_folder, T, scale, network_specifier);

[seg, cmap] = imread(fullfile(well_folder,[well,'_wts',network_specifier,'.tif']));

%%
rgb = label2rgb(seg,'jet',[0,0,0],'shuffle');
figure()
imshow(rgb)

saveas(gcf, fullfile('Figures/FullNetworks',[experiment,'_',magnification,'_',well,'_segmentation_.tif']))