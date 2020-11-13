//-----------------------------------------------------------------------
// STITCH A GRID OF IMAGES

// Macro for stitching tiles (created by high content microscope) into a fused image.
// Created 16-10-2020 by Lukas van den Heuvel.

// What this macro does (in chronological order):
// (1) Ask the user for experiment directory (=root), well, and number of stitched images on one axis (w).
// (2) Perform grid stitching;
// (3) Save the fused result as 16-bit tif.
//-----------------------------------------------------------------------

// User input
#@ File (label="Root", style="directory") root
#@ String (label="Well name") well 
#@ String (label="Width/height of fused image") w
#@ String (label="Overlap (%)") overlap
#@ String (label="How many bits?") bits
//#@ String (label="Number of color channels") c

close("*");
setBatchMode(true);

output_file = root+"/"+ well+"/"+well+"_fused.tif";

// If the fused image already exists, ask the user if they want to overwrite:
if (File.exists(output_file)){
	showMessageWithCancel("Fused image already exists!","A fused image of well "+well+" already exists.\nDo you want to continue and overwrite the old fused image?");
}

// Get dimensions of the first tile:
fname = root+"/"+well+"/tiles/tile_00.tif";
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
}

// Do grid stitching
print("Starting the stitching...");
setBatchMode(false);
run("Grid/Collection stitching", "type=[Grid: column-by-column] order=[Down & Right                ] grid_size_x="+w+" grid_size_y="+w+" tile_overlap="+overlap+" first_file_index_i=0 directory="+root+"/"+well+"/tiles file_names=tile_{ii}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");

// Convert to lower bount of bits
print("Converting to "+bits+"-bit. Please wait for a message box.");
run(bits+"-bit");

// Let the user change the LUT
title = "Set the right colors.";
message = "Change the LUT of the channels to set the right colors.\nPress OK when you are done.";
waitForUser(title, message);

// Save result
print("Saving the result...");
saveAs("Tiff", output_file);
print("Fused image is saved, you can now close it.");

