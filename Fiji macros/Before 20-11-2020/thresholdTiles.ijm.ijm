#@ File (label="root", style="directory") root
#@ String (label="wells you want to process (separated by commas)") well_string
#@ int w
#@ String (label="names of the channels you want to stitch the thresholded image for") ch_names_string

// Convert user inputs separated by commas into lists
well_list = split(well_string, ",");
ch_names_list = split(ch_names_string, ",");

for (l = 0; l < ch_names_list.length; t++) {

	close("*");
	
	// Define paths
	well = well_list[l];
	well_folder = root + "/" + well;
	tile_folder = well_folder + "/" + "tiles";
	threshold_folder = tile_folder + "/" + "thresholds";
	
	// Open fused image to get its size
	print("Opening fused image from well "+well);
	open(well_folder + "/" + well+"_fused.tif");
	run("8-bit");
	width = getWidth();
	height = getHeight();
	close();
	
	// Read configuration of tiles from fused image
	text = File.openAsString(tile_folder + "/TileConfiguration.txt");
	rows = split(text, "\n");
	setBatchMode(true);
	
	// Loop through thresholded channels
	for (c = 0; c < ch_names_list.length; c++) {
	
		// Make new fused image with the same size as the fused image
		fused_filename = well + "_th_" ch_names_list[c];
		newImage(fused_filename, "8-bit black", width, height, 1);
	
		// Loop through tiles
		for (i = 0; i < w*w; i++) {
			
			print("Processing well "+well+", channel "+ch_name+", tile "+d2s(i,0)+".");
	
			// Get tile name
			rows_split = split(rows[i+4], ";");
			file_name = rows_split[0];
			path = tile_folder + "/" + file_name;
			file_name_split = split(file_name, ".");
			tile_name = file_name_split[0];
	
			// Get tile coordinates
			tile_coordinate = rows_split[2];
			tile_coordinate_split = split(tile_coordinate, ",");
			x = tile_coordinate_split[0];
			x = parseInt( substring(x, 2, x.length) );	// Remove '(' at start
			y = tile_coordinate_split[1];
			y = parseInt( substring(y, 1, y.length-1) );	// Remove ')' at end
		
			// Open thresholded tile image
			tile_filename = tile_name + "_th_" + ch_name;
			open(threshold_folder + "/" + tile_filename + ".tif");
			
			wt = getWidth();
			ht = getHeight();
	
			// Loop through tile pixels
			for (xi = 0; xi < wt; xi++) {
				
				for (yi = 0; yi < ht; yi++) {
	
					// Get pixel value in thresholded tile
					selectWindow(tile_filename + ".tif");
					p = getPixel(xi, yi);
	
					// Put the pixel in the fused image
					selectWindow(fused_filename);
					setPixel(x+xi, y+yi, p);
		
				}
			}
			
			close(tile_filename + ".tif");
		}
		
		selectWindow(fused_filename);
		save(well_folder + "/" + fused_filename + ".tif")
		close();
	}
}
