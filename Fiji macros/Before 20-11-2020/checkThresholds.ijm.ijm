close("*") // close all images

// User input
#@ File (label="root", style="directory") root
#@ String (label="Well") well 

wellFolder = root + "/" + well;
tileFolder = wellFolder + "/tiles/thresholds";
textFileConfigurations = wellFolder + "/tiles/TileConfiguration.txt";

// Open the fused image
print("Reading the fused image...");
open(wellFolder + "/" + well+"_fused.tif");
run("8-bit");

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

// Clear everything outside the oval selection
print("Clearing the outside of the well...");
selectWindow(fname);
roiManager("Select", 0);
run("Clear Outside");

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

// Convert to RGB color and save it. 
// Will later be used by Matlab to draw the network on.
selectWindow(fname);
print("Converting to 8-bit RGB and saving it as seperate tif tile...");
run("RGB Color");
saveAs(root+"/"+well+"/"+well+"_fused_RGB.tif");
close();

// Make thresholded fused image
setBatchMode(true);
run("Tiles to Fused", "tiledirectory="+tileFolder+" welllocationstextfile="+textFileConfigurations);
setBatchMode("exit and display");
getDimensions(width, height, channels, slices, frames);
rename("threshold");
run("Split Channels");


// Clear the outside of the well in the thresholded images
for (i = 0; i < channels; i++) {
	selectWindow("C"+d2s(i+1,0)+"-threshold");
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
chNames = newArray(channels);
for (c = 0; c < channels; c++) {
	selectWindow("C"+d2s(c+1,0)+"-threshold");
	chNames[c] = getString("Name of this channel:", "dapi");
	waitForUser(title3, message3);
}

// Get the number of regions to clear
roiManager("deselect");
nROI = roiManager("count");

// Loop over thresholded images and clear the regions selected 
for (c = 0; c < channels; c++) {
	selectWindow("C"+d2s(c+1,0)+"-threshold");
	
	for (r = 0; r < nROI; r++) {
		roiManager("select", r);
		run("Clear", "slice");
		roiManager("deselect");
	}
	save(wellFolder + "/" + well+"_th_"+ chNames[c] +".tif");
	close();
}

print("All thresholded images are saved.")

