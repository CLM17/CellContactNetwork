function [Measurements, measurementNames, cellValues] = read_measurements_table(well_folder, scale, network_specifier)

    fName = ['cell_measurements', network_specifier, '.csv'];
    cellMeasurementsTable = readtable(fullfile(well_folder, fName));

    measurementNames = {'area', 'circularity', 'longness'};

    cellValues = cellMeasurementsTable.Mode;
    [cellValues,ind] = sort(cellValues);

    Measurements = struct();
    Measurements.area = cellMeasurementsTable.Area(ind) * scale^2;
    Measurements.circularity = cellMeasurementsTable.Circ_(ind);
    Measurements.longness = cellMeasurementsTable.Major(ind) ./ cellMeasurementsTable.Minor(ind);
    
end