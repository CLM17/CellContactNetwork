//-----------------------------------------------------------------------
// MAKE FUSED THRESHOLDS
// Created 20-11-2020 by Lukas van den Heuvel.
//
// This macro can be run after tiles are stitched into a fused image with stitch.ijm.
//
// It makes fused thresholded images, based on the thresholded tiles that were
// made with the macro prepareTiles.ijm and the tile positions that were made
// with stitch.ijm.
// It asks the user to make adjustments to the thresholded images (remove noise) and saves them.
//
// OUTPUTS
//--------- 

//		1) One RGB image called <well>_fused_RGB.tif.
// 		This is the fused image converted to RGB. It is used by Matlab to draw the network on.
//
// 		2) One thresholded image for each thresholded channel.
// 		The name of each thresholded image is <well>_th_<ch_name>.tif.
// 		Both <well> and <ch_name> are specified by the user when running this macro.
//
// All output images are saved automatically into the well_folder.
//-----------------------------------------------------------------------
//
//-----------------------------START SCRIPT------------------------------

// close all images
close("*") 

// User input
#@ File (label="root", style="directory") root
#@ String (label="Well") well 
#@ String (label="Names of the thresholded channels") chNamesString

// Specify the thresholds tile folder as tilefolder.
wellFolder = root + "\\" + well;
tileFolder = wellFolder + "\\tiles\\thresholds";
textFileConfigurations = wellFolder + "\\tiles\\TileConfiguration.txt";
print(tileFolder);

// Check if the thresholded imgages exist.
chNames = split(chNamesString, ",");
for (c = 0; c < chNames.length; c++) {
	thFile = wellFolder + "/" + well+"_th_"+ chNames[c] +".tif";
	if(!(File.exists(thFile))){
		exit("No image called "+well+"_th_"+ chNames[c] +".tif was found.\nCheck whether all thresholded files exist and whether you entered the names of the thresholded images correctly.");
	}
}

// Open the fused image
print("Reading the fused RGB image...");
open(wellFolder + "\\" + well+"_fused_RGB.tif");
getDimensions(width, height, channels, slices, frames);
fname =  getTitle();

// Adjust brightness and contrast
run("Brightness/Contrast...");
title0 = "STEP0: Adjust brightness / contrast";
message0 = "Adjust brightness/contrast.\nPress OK when you are done.";
waitForUser(title0, message0);

// Pre-set an oval selection
w =  getWidth();
h =  getHeight();
xc = w/2;
yc = h/2;
//makeOval(w/2 - w/4, h/2 - h/4, w/2, h/2);
setTool("oval");

// Let the user adjust the oval selection
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
selectWindow(fname);
//print("Converting to 8-bit RGB and saving it as seperate tif tile...");
//run("RGB Color");

// Clear everything outside the oval selection
print("Clearing the outside of the well...");
roiManager("Select", 0);
run("Clear Outside");

// Save the RGB color. Will later be used by Matlab to draw the network on.
saveAs(root+"\\"+well+"\\"+well+"_fused_RGB.tif");
//close();

// Clear the outside of the well for boundary image, if it exists
boundaryName = wellFolder + "/" + well +"_boundaries.tif";
if (File.exists(boundaryName)){
	open(boundaryName);
	roiManager("Select", 0);
	run("Clear Outside");
	waitForUser("Check threshold");
	save(boundaryName);
	close();
}

// Clear the outside of the well in the thresholded images
channels = chNames.length;
for (c = 0; c < channels; c++) {
	open(wellFolder + "/" + well+"_th_"+ chNames[c] +".tif");
	run("Grays");
	roiManager("Select", 0);
	run("Clear Outside");
}

run("Select None");
roiManager("Deselect");
roiManager("Delete");

title3 = "Make adjustments?";
message3 = "If you want to, make adjustments. Use color picker and brush tool to add foreground.\n Select regions that you want to remove from all thresholded images and add them to the ROI manager ('t').\n Press OK when you are done.";

// Loop over thresholded images and ask the user the name of this channel.
// Also, ask the user to select noise on the image and add it to ROi manager.

for (c = 0; c < channels; c++) {
	selectWindow(well+"_th_"+ chNames[c] +".tif");
	run("Select None");
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
	save(wellFolder + "\\" + well+"_th_"+ chNames[c] +".tif");
	close();
}

close("*");
print("All thresholded images are saved.");