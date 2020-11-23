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
//
//---------------------------START FUNCTIONS-----------------------------

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

function brute_force_stitch(tileFolder, w, width, height){
	//------------------------------------------------------
	// This function creates a file called textFileConfigurations.txt
	// inside the tileFolder.
	// The file contains the tile locations in a brute-force stitched fused image.
	// Brute force stitching means just pasting the tiles next to each other,
	// without computing overlap.
	//
	// When the file is created, the function makes the fused image
	// by running the plugin "Tiles to fused".
	//
	// INPUT PARAMETERS
	//-----------------
	//		(1) tileFolder: string
	//		Path to tiles.
	//
	//		(2) w: int
	//		Number of rows/columns in fused image.
	//
	//		(3,4) width, height: int
	//		Width and height (in pixels) of one tile.
	//------------------------------------------------------
	
	
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

//-----------------------------START SCRIPT------------------------------

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
	print("stitching...");
	run("Grid/Collection stitching", "type=[Grid: column-by-column] order=[Down & Right                ] grid_size_x="+w+" grid_size_y="+w+" tile_overlap="+overlap+" first_file_index_i=0 directory="+tileFolder+" file_names=tile_{ii}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");
}

// Brute-force stitching, only if normal stitching does not work
else if (brute_force){
	print("Brute-force stitching...");
	brute_force_stitch(tileFolder, w, width, height);
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

// Save raw fused result
print("Saving the result...");
saveAs("Tiff", output_file);

// Adjust brightness and contrast
run("Brightness/Contrast...");
title0 = "STEP0: Adjust brightness / contrast";
message0 = "Adjust brightness/contrast.\nPress OK when you are done.";
waitForUser(title0, message0);

// Let the user make an oval selection
setTool("oval");
title1 = "STEP1: Make the selection";
message1 = "\nFit a circular selection to the well borders while holding shift.\nPress OK when you are done.";
waitForUser(title1, message1);

if (roiManager("count") > 0){
	roiManager("Deselect");
	roiManager("Delete");
}
roiManager("Add");

// Measure the oval selection
print("Measuring well size and copying the results to your clipboard...");
String.resetBuffer; 
run("Set Measurements...", "centroid bounding fit redirect=None decimal=3");
run("Measure");

// Copy the center coordinates and diameter to clipboard
String.append(well);
String.append("\t");
String.append(getResult("X", nResults-1));
String.append("\t");
String.append(getResult("Y", nResults-1)); 
String.append("\t");
String.append(getResult("Width", nResults-1));
String.copy(String.buffer);

// Let the user paste the measurements to the xlsx file
title2 = "STEP2: Copy results";
message2 = "Measurement results are copied to your clipboard.\nPaste them into the .xlsx file 'Well locations'.\nPress OK when you are done.";
waitForUser(title2, message2);

// Convert to 8-bit RGB
setBatchMode(true);
selectWindow(fname);
print("Converting to 8-bit RGB and saving it as seperate tif tile...");
run("RGB Color");

// Clear everything outside the oval selection
print("Clearing the outside of the well...");
roiManager("Select", 0);
run("Clear Outside");

// Save the RGB color. Will later be used by Matlab to draw the network on.
saveAs(root+"/"+well+"/"+well+"_fused_RGB.tif");
close();
setBatchMode("exit and display");

// Stitch the thresholded images if the user wants to
//run("Tiles to Fused", "tiledirectory="+tileFolder+" welllocationstextfile="+textFileConfigurations+" fusedwidth="+width+" fusedheight="+height);
if(stitch_th){
	print("Starting the stitch of thresholded images...");
	run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by TileConfiguration] directory="+thresholdTileFolder+" layout_file=TileConfiguration.registered.txt fusion_method=[Max. Intensity] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");
	getDimensions(width, height, channels, slices, frames);
	rename("threshold");
	
	// Clear outside of well of thresholded image
	roiManager("Select", 0);
	run("Clear Outside");
	run("Select None");
	roiManager("Deselect");
	roiManager("Delete");

	// Pre-define messages to user
	title3 = "Make adjustments?";
	message3 = "If you want to, make adjustments. Use color picker and brush tool to add foreground.\n Select regions that you want to remove from all thresholded images and add them to the ROI manager ('t').\n Press OK when you are done.";

	// plit channels
	run("Split Channels");
	
	// Loop over thresholded images and ask the user the name of this channel.
	// Also, ask the user to select noise on the image and add it to ROi manager.

	roiManager("show all with labels");
	chNames = newArray(channels);
	for (c = 0; c < channels; c++) {
		selectWindow("C"+d2s(c+1,0)+"-threshold");
		//setBatchMode("show");
		run("Grays");
		run("Select None");
		
		// Let the user choose a name for this channel
		// The computer suggests a name based on the channel nr.
		if(c%2 == 0){
			suggestedName = "dapi";
		}
		else{
			suggestedName = "phalloidin";
		}
		chNames[c] = getString("Name of this channel:", suggestedName);
		rename(well+"_th_"+ chNames[c] +".tif");
		waitForUser(title3, message3);
	}

	// Get the number of regions to clear
	roiManager("deselect");
	nROI = roiManager("count");
	
	// Loop over thresholded images and clear the regions selected 
	for (c = 0; c < channels; c++) {
		selectWindow(well+"_th_"+ chNames[c] +".tif");
		
		for (r = 0; r < nROI; r++) {
			roiManager("select", r);
			run("Clear", "slice");
			roiManager("deselect");
		}
		print("Saving the thresholded result...");
		save(wellFolder + "/" + well+"_th_"+ chNames[c] +".tif");
		close();
	}
}

close("*");
print("Successfully stitched and processed well "+well+".");