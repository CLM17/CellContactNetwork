#@ File (label="root", style="directory") root
#@ String (label="Wells to process, separated by commas.") wells
#@ int (label="Magnification.") M 

close("*");
well_list = split(wells, ",");
Array.print(well_list);

if (M==10){
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

well = well_list[0];
folder = root + "/" + well;
open(folder + "/" + well+"_th_marker.tif");
open(folder + "/" + well+"_th_dapi.tif");

// Only take into account nuclei that overlap with marker
imageCalculator("Multiply", well+"_th_dapi.tif",well+"_th_marker.tif");

// Open and close the thresholded image
print("Opening and closing the DAPI image...");
run("Morphological Filters", "operation=Opening element=Disk radius="+selem);
run("Morphological Filters", "operation=Closing element=Disk radius="+selem);

// Separate touching nuclei with a distance transform watershed
print("Segmenting the nuclei with a distance transform watershed...");
run("Distance Transform Watershed", "distances=[Borgefors (3,4)] output=[16 bits] normalize dynamic="+dynamic+" connectivity=8");

setThreshold(1, 1E30);
setOption("BlackBackground", true);
run("Convert to Mask");

fName = getTitle();
run("Make Markers","outputdir="+folder+" currentData="+fName);

marker_name = getTitle();
print(marker_name);
	
// Marker controlled watershed
run("Marker-controlled Watershed", "input="+well+"_th_marker.tif marker="+marker_name+" mask="+well+"_th_marker.tif");

//run("Ultimate Points");
//setThreshold(1, 1E30);
//setOption("BlackBackground", true);
//run("Convert to Mask");

// Measure positions of center pixels and save the result
//run("Set Measurements...", "center redirect=None decimal=3");
//run("Analyze Particles...", "display clear add");
//saveAs("Results", folder + "/" + "nuclei_com.csv");
//run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
