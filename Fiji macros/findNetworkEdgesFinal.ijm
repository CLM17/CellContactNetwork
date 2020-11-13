#@ File (label="root", style="directory") root
#@ String (label="Wells to process, separated by commas.") wells
#@ int (label="Magnification.") M 
#@ boolean (label="Do you want to save the watershed result?") sWTS
#@ boolean (label="Do you want to remove the thresholded dapi image?") rDAPI
#@ boolean (label="Do you want to remove the thresholded marker image?") rMARKER

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

// Check if the network is already created.
// If it is, ask the user if they want to continue.
for (i = 0; i < well_list.length; i++) {
	well = well_list[i];
	folder = root + "/" + well;
	
	if ((File.exists(folder + "/" + "edges.csv")) | (File.exists(folder + "/" + "nuclei_com.csv"))){
		showMessageWithCancel("Network already created!","CSV files for the newtork in well "+well+" already exist.\nDo you want to continue and overwrite the old network?");
		if (File.exists(folder + "/" + "edges.csv")){
			File.delete(folder + "/" + "edges.csv");
		}
		if (File.exists(folder + "/" + "nuclei_com.csv")){
			File.delete(folder + "/" + "nuclei_com.csv");
		}
	}
}

for (i = 0; i < well_list.length; i++) {
	well = well_list[i];
	print("Processing well " + well);
	
	folder = root + "/" + well;
	open(folder + "/" + well+"_th_marker.tif");
	open(folder + "/" + well+"_th_dapi.tif");

	// Only take into account nuclei that overlap with marker
	imageCalculator("Multiply", well+"_th_dapi.tif",well+"_th_marker.tif");
	
	// Open and close the thresholded image
	print("Opening and closing the DAPI image...");
	run("Morphological Filters", "operation=Opening element=Disk radius="+selem);
	run("Morphological Filters", "operation=Closing element=Disk radius="+selem);

	close(well+"_th_dapi.tif");
	close(well+"_th_dapi-Opening");
	selectWindow(well+"_th_dapi-Opening-Closing");
	
	// Separate touching nuclei with a distance transform watershed
	print("Segmenting the nuclei with a distance transform watershed...");
	run("Distance Transform Watershed", "distances=[Borgefors (3,4)] output=[16 bits] normalize dynamic="+dynamic+" connectivity=8");

	close(well+"_th_dapi-Opening-Closing");
	selectWindow(well+"_th_dapi-Opening-Closing-dist-watershed");

	setThreshold(1, 1E30);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	
	fName = getTitle();
	run("Make Markers","outputdir="+folder+" currentData="+fName);

	marker_name = getTitle();
	print(marker_name);
	// Marker controlled watershed
	run("Marker-controlled Watershed", "input="+well+"_th_marker.tif marker="+marker_name+" mask="+well+"_th_marker.tif  calculate use");

	wts_name = well+"_th_marker-watershed.tif";
	selectWindow(wts_name);
	close("\\Others");
	
	if (sWTS) {
		print("Saving the watershed result...");
		save(folder + "/" + well+"_wts.tif");
	}
	if (rDAPI) {
		File.delete(folder + "/" + well + "_th_dapi.tif");
	}
	if (rMARKER) {
		File.delete(folder + "/" + well + "_th_marker.tif");
	}

	selectWindow(wts_name);
	// Find the edges of the network
	run("Network Creator", "outputdir="+folder+" currentData="+wts_name);

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
	saveAs("Results",  folder + "/cell_measurements.csv"); 
	
	print("Network created for well "+well+". You can now plot it in Matlab.");
	close("*");
}

print("All wells are processed.");
