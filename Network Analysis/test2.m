
allData = struct;
experiment = 'WKS024';
magnification = 'M20';
well = 'D02';
root = fullfile('../../Experiments', experiment, magnification);
well_folder = fullfile(root, well);
fieldSize = 1104;
network_specifier = '_ml';

xlsfileName = fullfile(root, 'Well locations.xlsx');
T = readtable(xlsfileName);

scale = calculate_scale(magnification, fieldSize);
allData = update_all_data(allData, experiment, magnification,...
                              well, well_folder, T, scale, network_specifier);