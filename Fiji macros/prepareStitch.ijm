//-----------------------------------------------------------------------
// PREPARE IMAGES FOR STITCHING
// Created 30-10-2020 by Lukas van den Heuvel.

// The CX7 high-content microscope arranges its images in a spiral,
// and stores each color channel in a separate image.

// This macro converts the filenames from a spiral grid into a column grid,
// so that Fiji can use it for column-grid stitching.
// Furthermore, it performs operations on seperate color channels 
// (e.g. background subtraction, contrast stretching, histogram equalization),
// and then merges all color channels into one image per tile.
 
// The output is a folder inside the root folder of the experiment
// with the name of the well.
// Inside this well_folder a folder "tiles" is created. It contains one image per tile.
// The names of the tiles are ordered in a column grid.
//-----------------------------------------------------------------------

//---------------------------START FUNCTIONS-----------------------------

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

function find_ch_operation_boolean(nr_channels, ch_operation_string){
	
	//----------------------------------------------------------------------------------
	// This function is used to make a user choose on which color 
	// channels they want to do an operation.

	// Examples of operations are: background subtraction, 
	// contrast enhancement, histogram equalization.

	// The user enters the channels into the dialog window,
	// seperated by commas.
	// This results in a string, example: ch_enhance_string = "1,2".
	// This function converts the string in a boolean with length
	// nr_channels (=number of channels).

	// So if nr_channels=3, the function outputs ch_enhance_boolean = [false, true, true].
	//----------------------------------------------------------------------------------

	// Split the string on commas:
	ch_operation_split = split(ch_operation_string, ",");
	// Initialize a boolean array with length nr_channels:
	ch_operation_boolean = newArray(nr_channels);

	// Loop over all channels
	for (nr = 0; nr < nr_channels; nr++){

		do_operation = false;
		
		// if this channel appears in ch_operation_string, make do_operation = true.
		for (i = 0; i < ch_operation_split.length; i++) {
			ch_nr = ch_operation_split[i];
			if (d2s(nr,0) == ch_nr){
				do_operation = true;
			}
		}
		ch_operation_boolean[nr] = do_operation;
	}
	return ch_operation_boolean;

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

function get_img_file_name(file0, well, nr, ch){
	
	//------------------------------------------------------
	// This function returns the filename of a raw tile
	// generated by the CX7 high content microscope.

	// Inputs
	// ------
	// file0: string
	// Name of the first file in the image folder. This serves as a template.

	// well: string
	// Name of the well, e.g. "B02".

	// nr: int
	// Number of the tile in the spiral grid.

	// ch: int
	// Channel number.

	// Output
	// ------
	// file_name: string
	// Name of the file specified by well, nr and ch.
	// Example: file_name = "MFGTMP_201029210001_B03f00d0.TIF"
	//------------------------------------------------------
	
	parts = split(file0, "_");
	nr_string = number_to_string(nr);
	file_name = parts[0] + "_" + parts[1] + "_" + well + "f" + nr_string + "d" + d2s(ch,0) + ".TIF";
	return file_name;
}
//---------------------------END FUNCTIONS-----------------------------

//---------------------------START SCRIPT------------------------------

// User input
#@ File (label="raw", style="directory") raw
#@ File (label="root", style="directory") root
#@ String well
#@ int w
#@ int (label="number of channels") nr_channels
#@ String(label="enhance contrast on channels") ch_enhance_string
#@ String(label="equalize histogram on channels") ch_equalize_string
#@ String(label="subtract background on channels") ch_subtract_backround_string
#@ int(label="rolling ball radius") rolling_ball_radius
#@ boolean (label="Do you want to display the images?") see

if (!(see)){
	setBatchMode(true);
}

// Close all open images
close("*");

input_folder = raw + "/" + well;
well_folder = root + "/" + well;
tile_folder = well_folder + "/tiles";

// Create well_folder if it did not exist already.
if (!(File.isDirectory(well_folder))){
	File.makeDirectory(well_folder);
	print("Created a new folder for this well.");
}
// Create tile_folder if it did not exist already.
// If it did exist, ask the user if they want to overwrite it.
if (!(File.isDirectory(tile_folder))){
	File.makeDirectory(tile_folder);
}
else{
	showMessageWithCancel("Tiles already created!","Tiles for well "+well+" were already created.\nDo you want to continue and overwrite the old tiles?");
}

// Convert user inputs about channel operations (contrast enhancement, etc.)
// into booleans.
ch_enhance = find_ch_operation_boolean(nr_channels, ch_enhance_string);
ch_equalize = find_ch_operation_boolean(nr_channels, ch_equalize_string);
ch_subtract_background = find_ch_operation_boolean(nr_channels, ch_subtract_backround_string);

// Make spiral and column grid arrays.
// Note that 2D arrays are not supported by Fiji macro language, so the arrays are 1D.
// If you index it them with matrix[x + y * width], they are effectively 2D.
column_grid = make_column_grid(w,w);
spiral_grid = make_spiral_grid(w,w);

// List all images in the raw input folder:
file_list = getFileList(input_folder);
// Get the name of the first file (it serves as a template).
file0 = file_list[0];

// Initialize positions on grid
x = 0;
y = 0;

// Loop over tiles in the grid
for (i = 0; i < w*w; i++) {
	print("Processing image "+d2s(i+1,0)+" out of "+d2s(w*w,0)+".");
	
	tile_nr = column_grid[x + y*w]; // column index corresponding with (x,y) position
	img_nr = spiral_grid[x + y*w];  // spiral index corresponding with (x,y) position
	
	// Initialize arrays with filenames and channel numbers.
	// They will later be used to merge the channels.
	file_names = newArray(nr_channels);
	ch_counts = newArray(nr_channels);

	// Loop over the channels in the image
	for (ch = 0; ch < nr_channels; ch++){
		file_name = get_img_file_name(file0, well, img_nr, ch);
		file_names[ch] = file_name;
		ch_counts[ch] = "c" + d2s(ch+1,0);
		
		// Open image
		open(input_folder + "/" + file_name);

		// Perform the operations if the user asked for it:
		if ((ch_enhance[ch] == true) && (ch_equalize[ch] == false)){
			run("Enhance Contrast...", "saturated=0.1 normalize");
		}
		if ((ch_enhance[ch] == true) && (ch_equalize[ch] == true)){
			run("Enhance Contrast...", "saturated=0.1 normalize equalize");
		}
		if (ch_subtract_background[ch] == true){
			run("Subtract Background...", "rolling="+d2s(rolling_ball_radius,0)+" disable");
		}
	}

	// Find the right color for right channel (reverse order)
	channel_string = "";
	for (ch = 0; ch < nr_channels; ch++){
		 channel_string = channel_string + ch_counts[ch] + "=" + file_names[nr_channels - ch - 1] + " ";
	}
	// Merge the channels
	run("Merge Channels...", channel_string + "create");

	// Give the image the correct tile number.
	// This number is based on the column grid.
	tile_nr_string = number_to_string(tile_nr);
	tile_name = "tile_" + tile_nr_string;

	// Save the image and close it.
	save(tile_folder + "/" + tile_name);
	close();

	// Update position
	y = y + 1;
	if (y == w){
		y = 0;
		x = x + 1;
	}
}
print("Finished the preparation, you can now start stitching this well.");