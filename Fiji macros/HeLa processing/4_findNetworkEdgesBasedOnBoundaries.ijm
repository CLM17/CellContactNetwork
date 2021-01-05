//-----------------------------------------------------------------------
// FIND NETWORK EDGES
// Created 20-11-2020 by Lukas van den Heuvel.
//
// This macro processes the fused thresholded images created in makeFusedThresholds.ijm.
// It then creates the cell network by finding which cells are neighbours.
//
// The macro requires no user intervention and can be run for multiple wells.
//
// These are the processing operations in chronological order:
// 1) Multiply th_DAPI and th_SOMA <-- new th_DAPI;
// 2) Open and close th_DAPI <-- new th_DAPI;
// 3) Seperate touching nuclei on th_DAPI with a distance transform watershed <-- wts_DAPI;
// 4) Put pixels at the centers of mass of the nuclei <-- markers_DAPI;
// 5) Marker-controlled watershed with th_SOMA as mask and markers_DAPI as markers.
// 6) Measures cell properties (area, circularity, ...) and stores the result in cell_measurements.csv.
// 7) Finds the network edges with the plugin "Network Creator" and stores the result in edges.csv.
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

#@ File (label="root", style="directory") root
#@ String (label="Wells to process, separated by commas.") wells
#@ String (label="Dapi channel name.") nDAPI
#@ String (label="Soma channel name.") nSOMA
#@ int (label="Magnification.") M 
#@ boolean (label="Do you want to save the watershed result?") sWTS
#@ boolean (label="Do you want to remove the thresholded dapi image?") rDAPI
#@ boolean (label="Do you want to remove the thresholded soma image?") rSOMA

close("*");
//setBatchMode(true);
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

// Check if thresholded files exist for all wells.
// Also, check if the network is already created.
// If it is, ask the user if they want to continue.
for (i = 0; i < well_list.length; i++) {
	well = well_list[i];
	wellFolder = root + "/" + well;

	if(!(File.exists( wellFolder + "/" + well+"_boundaries"+".tif" ))){
		exit(well+"_boundaries"+".tif not found. Make a prediction with the Python code 'make_nn_prediction.py'");
	}

	if(!(File.exists( wellFolder + "/" + well+"_th_"+nSOMA+".tif" ))){
		exit(well+"_th_"+nSOMA+".tif not found.");
	}
	if(!(File.exists( wellFolder + "/" + well+"_th_"+nDAPI+".tif" ))){
		exit(well+"_th_"+nDAPI+".tif not found.");
	}
	
	if ((File.exists(wellFolder + "/" + "edges.csv")) | (File.exists(wellFolder + "/" + "nuclei_com.csv"))){
		showMessageWithCancel("Network already created!","CSV files for the newtork in well "+well+" already exist.\nDo you want to continue and overwrite the old network?");
		if (File.exists(wellFolder + "/" + "edges.csv")){
			File.delete(wellFolder + "/" + "edges.csv");
		}
		if (File.exists(wellFolder + "/" + "nuclei_com.csv")){
			File.delete(wellFolder + "/" + "nuclei_com.csv");
		}
	}
}

// Loop over all wells
for (i = 0; i < well_list.length; i++) {
	well = well_list[i];
	print("Processing well " + well);
	
	wellFolder = root + "/" + well;
	fileNameSOMA = well+"_th_"+nSOMA;
	fileNameDAPI = well+"_th_"+nDAPI;
	fileNameBoundary = well+"_boundaries";
	open(wellFolder + "/" + fileNameSOMA+".tif");
	run("Select None");
	setThreshold(180, 255);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	
	open(wellFolder + "/" + fileNameDAPI+".tif");
	run("Select None");
	open(wellFolder + "/" + fileNameBoundary+".tif");
	run("Select None");

	// Only take into account nuclei that overlap with marker
	//imageCalculator("Multiply", fileNameDAPI+".tif", fileNameSOMA+".tif");
	
	// Open and close the thresholded image
	selectWindow(fileNameDAPI+".tif");
	print("Opening and closing the DAPI image...");
	run("Morphological Filters", "operation=Opening element=Disk radius="+selem);
	run("Morphological Filters", "operation=Closing element=Disk radius="+selem);

	close(fileNameDAPI+".tif");
	close(fileNameDAPI+"-Opening");
	selectWindow(fileNameDAPI+"-Opening-Closing");
	
	// Separate touching nuclei with a distance transform watershed
	print("Segmenting the nuclei with a distance transform watershed...");
	run("Distance Transform Watershed", "distances=[Borgefors (3,4)] output=[16 bits] normalize dynamic="+dynamic+" connectivity=8");

	close(fileNameDAPI+"-Opening-Closing");
	selectWindow(fileNameDAPI+"-Opening-Closing-dist-watershed");

	setThreshold(1, 1E30);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	
	fName = getTitle();
	run("Make Markers","outputdir="+wellFolder+" currentData="+fName);

	marker_name = getTitle();

	// Set threshold on boundary image
	selectWindow(fileNameBoundary + ".tif");
	run("Select None");
	setAutoThreshold("Yen dark");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	
	imageCalculator("Subtract create", fileNameSOMA+".tif", fileNameBoundary+".tif");
	run("Morphological Filters", "operation=Erosion element=Disk radius=2");
	rename("cell_bodies");

	print("Eroded");
	// Marker controlled watershed
	run("Marker-controlled Watershed", "input=cell_bodies marker="+marker_name+" mask=cell_bodies calculate use");
	cell_bodies_name = getTitle();
	print("Marker controlled watershed");
	//close("cell_bodies");
	run("Marker-controlled Watershed", "input="+fileNameSOMA+".tif marker="+cell_bodies_name+" mask="+fileNameSOMA+".tif  calculate use");
	wts_name = getTitle();
	close(cell_bodies_name);

	selectWindow(wts_name);
	close("\\Others");

	// Save watershed (if the user asked for it)
	if (sWTS) {
		print("Saving the watershed result...");
		save(wellFolder + "/" + well+"_wts.tif");
	}
	// Remove dapi or soma threshold if the user asked for it.
	if (rDAPI) {
		File.delete(wellFolder + "/" + fileNameDAPI +".tif");
	}
	if (rSOMA) {
		File.delete(wellFolder + "/" + fileNameSOMA +".tif");
	}

	selectWindow(wts_name);
	// Find the edges of the network
	run("Network Creator", "outputdir="+wellFolder+" currentData="+wts_name);

	// Measure cell properties
	selectWindow(wts_name);
	run("Duplicate...", " ");
	
	setThreshold(1, 1E30);
	setOption("BlackBackground", true);
	run("Convert to Mask");

	if (roiManager("count") > 0){
		roiManager("delete");
	}
	if (nResults > 0){
		run("Clear Results");
	}

	run("Analyze Particles...", "add");
	selectWindow(wts_name);
	roiManager("Show All");

	run("Set Measurements...", "area modal bounding fit shape redirect=None decimal=3");
	roiManager("Measure");
	saveAs("Results",  wellFolder + "/cell_measurements.csv"); 
	
	print("Network created for well "+well+". You can now plot it in Matlab.");
	close("*");
}

print("All wells are processed.");
