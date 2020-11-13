// Use the 'Oval' tool to draw an oval on the image.
// To make the oval a circle, hold Ctrl + Shift while drawing.
// If your oval marks the well edges perfectly, run this macro and ...
// copy the results in the excel file 'Well locations.xlsx'.
// Note: X,Y are the locations of the center, Major, Minor is the diameter.

close("*")
#@ String well 
#@ File (label="folder", style="directory") folder
#@ boolean (label="Do you want to save the thresholded DAPI image?") sDAPI
#@ boolean (label="Do you want to save the thresholded marker image?") sMARKER

open(folder + "/" + well+"_fused.tif");

run("Duplicate...", "duplicate");
fname =  getTitle();
w =  getWidth();
h =  getHeight();
xc = w/2;
yc = h/2;
makeOval(w/2 - w/4, h/2 - h/4, w/2, h/2);

title1 = "STEP1: Make the selection";
message1 = "Select the Oval tool (tip left corner).\nFit the circular selection to the well borders while holding ctrl+shift.\nPress OK when you are done.";
waitForUser(title1, message1);
run("Clear Outside");

run("Set Measurements...", "centroid fit redirect=None decimal=3");
run("Measure");

title2 = "STEP2: Copy results";
message2 = "Copy the Measurement results to the .xlsx file 'Well locations'.\nPress OK when you are done.";
waitForUser(title2, message2);

run("Select None");
run("Gaussian Blur...", "sigma=2");
run("Split Channels");
close("C2-"+fname);

run("Threshold...");
setOption("BlackBackground", true);
title3 = "STEP3: Threshold DAPI";
message3 = "Set the right thresholding method for DAPI.\nIn case of contamination, select it with the Polygon selection tool\n(top left) and remove it (edit > clear).\nReport the thresholding method in the .xlsx file 'Well locations'.\n Press OK when you are done.";
waitForUser(title3, message3);
//run("Threshold...");
run("Convert to Mask");

if (sDAPI) {
	save(folder + "/" + well+"_th_dapi.tif");
}
close();
//close();

run("Threshold...");
setOption("BlackBackground", true);

title4 = "STEP4: Threshold marker";
message4 = "Set the right thresholding method for the marker.\nnReport the thresholding method in the .xlsx file 'Well locations'.";
waitForUser(title4, message4);
run("Convert to Mask");

if (sMARKER) {
	save(folder + "/" + well+"_th_maker.tif");
}
close();
print("You finished the pre-processing of well "+well+".")


 