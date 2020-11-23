close("*") // close all images

#@ File (label="root", style="directory") root
#@ String well

folder = root + "/" + well
open(folder + "/" + well+"_th_DAPI.tif");
run("Duplicate...", " ");

// Open the image
run("Erode");
run("Erode");
run("Erode");
run("Dilate");
run("Dilate");
run("Dilate");

// Separate touching nuclei with a distance transform watershed
run("Distance Transform Watershed", "distances=[Borgefors (3,4)] output=[32 bits] normalize dynamic=1 connectivity=8");
setThreshold(1, 1E30);
setOption("BlackBackground", true);
run("Convert to Mask");
save(folder + "/" + well+"_th_dapi.tif");

close("\\Others");// close all images