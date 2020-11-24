function [Measurements, measurementNames, cellValues] = read_measurements_table(well_folder, scale)

    cellMeasurementsTable = readtable(fullfile(well_folder,'cell_measurements.csv'));

    measurementNames = {'area', 'circularity', 'longness'};

    cellValues = cellMeasurementsTable.Mode;
    [cellValues,ind] = sort(cellValues);

    Measurements = struct();
    Measurements.area = cellMeasurementsTable.Area(ind) * scale^2;
    Measurements.circularity = cellMeasurementsTable.Circ_(ind);
    Measurements.longness = cellMeasurementsTable.Major(ind) ./ cellMeasurementsTable.Minor(ind);
    
end