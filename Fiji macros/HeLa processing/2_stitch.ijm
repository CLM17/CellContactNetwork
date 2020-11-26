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

function initialize_position(width, height){

	//-----------------------------------------------------
	// Finds the center of a spiral grid, i.e. the 
	// (x,y) location in the spiral grid with value 0.

	// Inputs
	// width & height: the dimensions of the spiral grid.

	// Output
	// init_position: array of length 2. 
	// Index 0 is the x position, index 1 is the y position.
	//-----------------------------------------------------
	
	init_position = newArray(2);
	if (width%2 == 0){
		x = width / 2 - 1;
	}
	else{
		x = floor(width / 2);
	}
	if (height%2 == 0){
		y = height / 2 - 1;
	}
	else{
		y = floor(height / 2);
	}
	init_position[0] = x;
	init_position[1] = y;

	return init_position;

}

function compare_arrays(arr1, arr2){

	//----------------------------------------------------------
	// Function compares array arr1 and array arr2 and
	// returns true if they are equal, and false if they are not.
	//----------------------------------------------------------
	
	same = true;
	l = arr1.length;
	for (i = 0; i < l; i++){
		if (arr1[i] != arr2[i]){
			same = false;
		}
	}
	return same;
}

function turn_right(current_direction){

	//-----------------------------------------------------
	// This function is called by make_spiral_grid().
	
	// It takes as input the current direction (north, south, west and east),
	// and outputs the direction after turning right:
	// north->east, south->west, east->south, west->north.
	
	// All directions are arrays of length 2:
	// index 0 is dx, index 1 is dy.
	//-----------------------------------------------------
	
	NORTH = newArray(0,-1);
	S = newArray(0,1);
	W = newArray(-1,0);
	E = newArray(1,0);

	if (compare_arrays(current_direction, NORTH)){
		new_direction = E;
	}
	else if (compare_arrays(current_direction, S)){
		new_direction = W;
	}
	else if (compare_arrays(current_direction, E)){
		new_direction = S;
	}
	else if (compare_arrays(current_direction, W)){
		new_direction = NORTH;
	}
	return new_direction;
	
}

function make_spiral_grid(width,height){
	//-----------------------------------------------------
	// This function makes a spiral grid array.
	// Note that 2D arrays are not supported by Fiji macro language.
	// The output will be a 1D array with length width*height.
	
	// If you index it with matrix[x + y * width], 
	// the matrix is effectively 2D.
	//-----------------------------------------------------
	
	// Check if dimenstions are at least 1:
	if (width < 1 || height < 1){
		return;
	}

	// Directions to walk (dx, dy) are arrays:
	NORTH = newArray(0,-1);
	S = newArray(0,1);
	W = newArray(-1,0);
	E = newArray(1,0);

	// Initial position:
	init_position = initialize_position(width, height);
	x = init_position[0];
	y = init_position[1];
	// We want to start walking to the west. 
	// This means our initial direction is north:
	// we then turn right immediately, and end up going west.
	direction = NORTH;
	dx = direction[0];
	dy = direction[1];

	// Initialize matrix:
	matrix = newArray(w*w);
	Array.fill(matrix, NaN);
	count = 0;

	while (true){
		// Fill matrix at current position with value count:
		matrix[x + y*width] = count;
		
		count = count + 1;
		
		// Try to turn right:
		new_direction = turn_right(direction);
		new_dx = new_direction[0];
		new_dy = new_direction[1];
		new_x = x + new_dx;
		new_y = y + new_dy;

		// Turn right if you are not yet at the boundaries,
		// and if the position at your right-hand side is not yet visited:
		if ((0<=new_x) & (new_x<width) & (0<=new_y) & (new_y<height) & (isNaN(matrix[new_x + new_y*width]))){
			x = new_x;
			y = new_y;
			dx = new_dx;
			dy = new_dy;
			direction = new_direction;
		}
		// If not, go straight:
		else{
			x = x + dx;
			y = y + dy;
			// If you are at the boundary, stop the process and return the matrix:
			if (!((0 <= x) & (x < width)) & (0 <= y) & (y < height)){
				return matrix;
			}
		}
	}
}

function make_column_grid(width, height){

	//-----------------------------------------------------
	// This function makes a column grid array.
	// Note that 2D arrays are not supported by Fiji macro language.
	// The output will be a 1D array with length width*height.
	
	// If you index it with matrix[x + y * width], 
	// the matrix is effectively 2D.
	//-----------------------------------------------------
	
	matrix = newArray(width*height);
	Array.fill(matrix, NaN);
	count = 0;
	x = 0;
	y = 0;
	for (i = 0; i < width*height; i++) {
		matrix[x + y*width] = count;
		count = count + 1;
		y = y + 1;
		if (y == height){
			y = 0;
			x = x + 1;
		}
	}
	return matrix;
}

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

close("*");
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
	
	// Get dimensions of first thresholded tile
	fname = thresholdTileFolder + "/tile_00.tif";
	open(fname);
	getDimensions(width, height, channels, slices, frames);
	
	print("Starting the stitch of thresholded images...");
	run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by TileConfiguration] directory="+thresholdTileFolder+" layout_file=TileConfiguration.registered.txt fusion_method=[Max. Intensity] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");
	rename("threshold");
	run("Split Channels");

	spiral_grid = make_spiral_grid(w,w);
	setJustification("center");
	setColor("red");
	setFont("SansSerif", 25);

	// Run over the channels
	chNames = newArray(channels);
	for (c = 0; c < channels; c++) {
		selectWindow("C"+d2s(c+1,0)+"-threshold");
		setBatchMode("show");
		run("Grays");

		// Draw ROIs
		x = 0;
		y = 0;
		for (i = 0; i < w*w; i++){
			img_nr = spiral_grid[x + y*w];
			drawString(number_to_string(img_nr), (x+0.5)*width, (y+0.5)*height);
			
			// Update position
			y = y + 1;
			if (y == w){
				y = 0;
				x = x + 1;
			}
		}
	
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
