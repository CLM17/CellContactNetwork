//-----------------------------------------------------------------------
// STITCH A GRID OF IMAGES

// Macro for stitching tiles (created by high content microscope) into a fused image.
// Created 16-10-2020 by Lukas van den Heuvel.
//
// There is also the possibility to do brute-force stitching (= pasting tiles against each other
// without calculating any correlations).
//
// What this macro does (in chronological order):
// (1) Ask the user for experiment directory (=root), well, and number of stitched images on one axis (w).
// (2) Perform grid stitching and save the result.;
// (4) Stitch the thesholded images, if the user wants to. 
//-----------------------------------------------------------------------

function number_to_string(nr){
	//------------------------------------------------------
	// This function converts an integer nr into a string.
	// Examples: 0 --> "00", 3 --> "03", 11 --> "11", 102 --> "102".
	//------------------------------------------------------
	
	if (nr < 10){
		nr_string = "0" + d2s(nr, 0);
	}
	else{
		nr_string = d2s(nr, 0);
	}
	return nr_string;
}

// User input
#@ File (label="Root", style="directory") root
#@ String (label="Well name") well 
#@ String (label="Width/height of fused image") w
#@ String (label="Overlap (%)") overlap
#@ String (label="How many bits?") bits
#@ boolean (label="Stitch the thresholded tiles?") stitch_th
#@ boolean (label="Brute-force stitch? (not recommended)") brute_force
//#@ String (label="Number of color channels") c

//close("*");
setBatchMode(true);

wellFolder = root + "/" + well;
output_file = wellFolder + "/"+well+"_fused.tif";

// If the fused image already exists, ask the user if they want to overwrite:
if (File.exists(output_file)){
	showMessageWithCancel("Fused image already exists!","A fused image of well "+well+" already exists.\nDo you want to continue and overwrite the old fused image?");
}

// Get dimensions of the first tile:
fname = root+"/"+well+"/tiles/tile_00.tif";
tileFolder = root + "/" + well + "/tiles";
thresholdTileFolder = tileFolder + "/thresholds";
open(fname);
getDimensions(width, height, channels, slices, frames);
close();

// If the images were saved as z-stack instead of color channels,
// convert them into color channels.
if (slices > channels){

	for (i = 0; i < w*w; i++) {
		if (i < 10) {
			n = "0"+d2s(i,0);
		}
		else{
			n = d2s(i,0);
		}
		fname = root+"/"+well+"/tiles/tile_"+n+".tif";
		open(fname);
		print("Converting z-stack to color channel in tile "+n);
		run("Properties...", "channels="+d2s(slices,0)+" slices=1 frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
		save(fname);
		close();
		}
	channels = slices;
}

// Do grid stitching
print("Starting the stitching...");

// Normal stitching
if (!(brute_force)){
	print("stitching");
	run("Grid/Collection stitching", "type=[Grid: column-by-column] order=[Down & Right                ] grid_size_x="+w+" grid_size_y="+w+" tile_overlap="+overlap+" first_file_index_i=0 directory="+tileFolder+" file_names=tile_{ii}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");
}

// Brute-force stitching, only if normal stitching does not work
else if (brute_force){

	// Initialize a string to store the tile configurations
	outputString = "";
	outputString = outputString + "# Define the number of dimensions we are working on\n";
	outputString = outputString + "dim = 2\n";
	outputString = outputString + "\n";
	outputString = outputString + "# Define the image coordinates\n";

	// Initialize tile positions on the grid
	x = 0;
	y = 0;

	// Loop over tiles
	for (i = 0; i < w*w; i++) {
		
		nr = number_to_string(i);
		outputString = outputString + "tile_" + nr + ".tif; ; ";
		outputString = outputString + "(" + d2s(x,1) + ", " + d2s(y,1) + ")\n";

		y = y + height;
		if (y == w*height){
			y = 0;
			x = x + width;
		}
	}
	textFileConfigurations = tileFolder + "/TileConfiguration.registered.txt";
	File.saveString(outputString, textFileConfigurations);
	fusedWidth = d2s(w*width, 0);
	fusedHeight = d2s(w*height, 0);
	run("Tiles to Fused", "tiledirectory="+tileFolder+" welllocationstextfile="+textFileConfigurations+" fusedwidth="+fusedWidth+" fusedheight="+fusedHeight);
}

// Copy registered tile configuration to thresholds folder
File.copy(tileFolder + "/TileConfiguration.registered.txt", tileFolder + "/thresholds/TileConfiguration.registered.txt")

// Convert to lower bount of bits
print("Converting to "+bits+"-bit. Please wait for a message box.");
run(bits+"-bit");

// Let the user change the LUT
setBatchMode("show");
title = "Set the right colors.";
message = "Change the LUT of the channels to set the right colors.\nPress OK when you are done.";
waitForUser(title, message);

// Save result
print("Saving the result...");
saveAs("Tiff", output_file);
close();

// Stitch the thresholded images if the user wants to
//run("Tiles to Fused", "tiledirectory="+tileFolder+" welllocationstextfile="+textFileConfigurations+" fusedwidth="+width+" fusedheight="+height);
if(stitch_th){
	print("Starting the stitch of thresholded images...");
	run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by TileConfiguration] directory="+thresholdTileFolder+" layout_file=TileConfiguration.registered.txt fusion_method=[Max. Intensity] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");
	getDimensions(width, height, channels, slices, frames);
	rename("threshold");
	run("Split Channels");

	// Run over the channels
	chNames = newArray(channels);
	for (c = 0; c < channels; c++) {
		selectWindow("C"+d2s(c+1,0)+"-threshold");
		setBatchMode("show");
		run("Grays");
		
		// Let the user choose a name for this channel
		// The computer suggests a name based on the channel nr.
		if(c%2 == 0){
			suggestedName = "dapi";
		}
		else{
			suggestedName = "phalloidin";
		}
		chNames[c] = getString("Name of this channel:", suggestedName);
		print("Saving the thresholded result...");
		save(wellFolder + "/" + well+"_th_"+ chNames[c] +".tif");
		close();
	}
}

setBatchMode("exit and display");
print("Successfully stitched well "+well+".");
