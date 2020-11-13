// Macro for preprocessing images of HeLa cell cultures.
// Created 12/10/2020 by Lukas.
// To run this macro, make sure you have the MorphoLibJ plugin installed:
// https://github.com/ijpb/MorphoLibJ/

// What this macro does (in chronological order):
// (1) it lets the user select the well;
// (2) it crops the outside of the well and measures its center and diameter;
// (3) applies a Gaussian blur to the 3D image;
// (4) splits the channels;
// (5) lets the user apply a threshold on the individual channel;
// (6) saves the thresholded images.

close("*") // close all images

// User input
#@ File (label="root", style="directory") root
#@ String (label="Well") well 
#@ Character (label="DAPI color channel") chDAPI
#@ Character (label="Marker color channel") chMARKER
#@ int (label="Magnification.") M

if (M==10){
	sigma = "2";
	print("Magnification 10x");
}
else if (M==20){
	sigma = "4";
	print("Magnification 20x");
}
else{
	exit("Invalid magnification.");
}

folder = root + "/" + well;

// Check if thresholded images are already creted.
// If it is, ask the user if they want to continue.
if ((File.exists(folder + "/" + well + "_th_dapi.tif")) | (File.exists(folder + "/" + well + "_th_marker.tif"))){
	showMessageWithCancel("Already processed!","Thresholded DAPI and/or marker images for well "+well+" already exist.\nDo you want to continue and overwrite the old image(s)?");
}

// Open the image
print("Reading the fused image...");
open(folder + "/" + well+"_fused.tif");
run("8-bit");

getDimensions(width, height, channels, slices, frames);
if (slices > channels){
	print("Converting z-stack to color channels.");
	run("Properties...", "channels="+d2s(slices,0)+" slices=1 frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
}

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
message1 = "Select the Oval tool (tip left corner).\nFit the circular selection to the well borders while holding ctrl+shift.\nPress OK when you are done.";
waitForUser(title1, message1);

if (roiManager("count") > 0){
	roiManager("Deselect");
	roiManager("Delete");
}
roiManager("Add");

// Make a duplication
print("Duplicating the image...");
run("Select None");
run("Duplicate...", "duplicate");
fname =  getTitle();

roiManager("Select", 0);

// Clear everything outside the oval selection
print("Clearing the outside of the well...");
selectWindow(fname);
run("Clear Outside");

// Measure the oval selection
print("Measuring well size and copying the results to your clipboard...");
String.resetBuffer; 
run("Set Measurements...", "centroid bounding fit redirect=None decimal=3");
run("Measure");

// Get the bounding box results (for later)
bx = getResult("BX", nResults-1);
by = getResult("BY", nResults-1);
d = getResult("Width", nResults-1);

// Copy the center coordinates and diameter to clipboard
String.append(well);
String.append("\t");
String.append(getResult("X", nResults-1));
String.append("\t");
String.append(getResult("Y", nResults-1)); 
String.append("\t");
String.append(d);
String.copy(String.buffer);

// Let the user paste the measurements to the xlsx file
title2 = "STEP2: Copy results";
message2 = "Measurement results are copied to your clipboard.\nPaste them into the .xlsx file 'Well locations'.\nPress OK when you are done.";
waitForUser(title2, message2);

// Convert to RGB color and save it. 
// Will later be used by Matlab to draw the network on.
setBatchMode(true);
selectWindow(fname);
print("Converting to 8-bit RGB and saving it as seperate tif tile...");
run("RGB Color");
saveAs(root+"/"+well+"/"+well+"_fused_RGB.tif");
close();
setBatchMode(false);if (M==10){
	selem = 4;
	dynamic = 2;
	print("Magnification 10x");
}
else if (M==20){
	selem = 8;
	dynamic = 4;
	print("Magnification 20x");
}
else{
	exit("Invalid magnification.");
}

// Apply a Gaussian blur to the image
print("Applying a Gaussian blur and splitting the channels...");
selectWindow(fname);
run("Select None");
run("Gaussian Blur...", "sigma="+sigma);

// Split the color channels and select the DAPI channel
run("Split Channels");
selectWindow("C"+chDAPI+"-"+fname);

// Let the user threshold the DAPI channel (and remove contamination)
//makeOval(bx, by, d, d);
run("Threshold...");
setOption("BlackBackground", true);
title3 = "STEP3: Threshold DAPI";
message3 = "Set the right thresholding method for DAPI.\nIn case of contamination, select it with the Polygon selection tool\n(top left) and remove it (edit > clear).\nReport the thresholding method in the .xlsx file 'Well locations'.\n Press OK when you are done.";
waitForUser(title3, message3);
selectWindow("C"+chDAPI+"-"+fname);
run("Convert to Mask");
title3a = "Make adjustments?";
message3a = "If you want to, make adjustments. Use color picker and brush tool.\n Press OK when you are done.";
waitForUser(title3a, message3a);

selectWindow("C"+chDAPI+"-"+fname);

// Select the marker channel image and make sure there are no selections on it
selectWindow("C"+chMARKER+"-"+fname);
run("Select None");

// Let the user threshold the marker channel (and remove contamination)
//makeOval(bx, by, d, d);
run("Threshold...");
setOption("BlackBackground", true);
title4 = "STEP4: Threshold marker";
message4 = "Set the right thresholding method for the marker.\nIn case of contamination, select it with the Polygon selection tool\n(top left) and remove it (edit > clear).\nReport the thresholding method in the .xlsx file 'Well locations'.\n Press OK when you are done.";
waitForUser(title4, message4);
selectWindow("C"+chMARKER+"-"+fname);
run("Convert to Mask");
title4a = "Make adjustments?";
message4a = "If you want to, make adjustments. Use color picker and brush tool.\n Press OK when you are done.";
waitForUser(title4a, message4a);

// Save the thresholded marker image
selectWindow("C"+chMARKER+"-"+fname);
print("Saving the thresholded marker image...");
save(folder + "/" + well+"_th_marker.tif");
close();

// Save the thresholded dapi image
selectWindow("C"+chDAPI+"-"+fname);
print("Saving the segmented DAPI image...");
save(folder + "/" + well+"_th_dapi.tif");

// Close all images except the original fused image
selectWindow(well+"_fused.tif");
close("\\Others");
print("Finished processing well "+well+"! You can now find the network edges.");

