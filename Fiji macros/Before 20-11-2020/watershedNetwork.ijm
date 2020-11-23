close("*");

#@ File (label = "root", style = "directory") root
#@ String well
folder = root + "/" + well;
#@ boolean (label="Do you want to save the watershed result?") sWTS

if ((File.exists(folder + "/" + "edges.csv")) | (File.exists(folder + "/" + "nuclei_com.csv"))){
	showMessageWithCancel("Network already created!","CSV files for the newtork in well "+well+" already exist.\nDo you want to continue and overwrite the old network?");
	if (File.exists(folder + "/" + "edges.csv")){
		File.delete(folder + "/" + "edges.csv")
	}
	if (File.exists(folder + "/" + "nuclei_com.csv")){
		File.delete(folder + "/" + "nuclei_com.csv")
	}
}

open(folder + "/" + well + "_th_dapi.tif");
open(folder + "/" + well + "_th_marker.tif");

selectWindow(well+"_th_dapi.tif");

// Find center pixels of nuclei
run("Ultimate Points");
setThreshold(1, 1E30);
setOption("BlackBackground", true);
run("Convert to Mask");

// Measure positions of center pixels and save the result
run("Set Measurements...", "center redirect=None decimal=3");
run("Analyze Particles...", "display clear add");
saveAs("Results", folder + "/" + "nuclei_com.csv");

// Marker controlled watershed
run("Marker-controlled Watershed", "input="+well+"_th_marker.tif marker="+well+"_th_dapi.tif mask="+well+"_th_marker.tif binary");
if (sWTS) {
	print("Saving the watershed result...");
	save(folder + "/" + well+"_wts.tif");
}

// Find the edges of the network
run("Network Creator", "outputdir="+folder+" currentData="+well+"_th_marker-watershed.tif");
print("Network created, analysis finished.");

//run("Set Measurements...", "min redirect=None decimal=3");
//run("Measure");
//max = getResult("Max", nResults-1);

//setBatchMode(true);
//run("Clear Results");
//roiManager("reset");
//run("Set Measurements...", "area centroid center perimeter fit feret's area_fraction redirect=None decimal=3");
//for (i = 1; i < max+1; i++) {
//	selectWindow(well+"_th_marker-watershed.tif");
//	run("Duplicate...", " ");
//	setThreshold(i, i);
//	setOption("BlackBackground", true);
//	run("Convert to Mask");
//	run("Analyze Particles...", "add");
//	close();
//}
//roiManager("Show All");
//roiManager("Measure");
//saveAs("Results", folder + "/" + "cell_properties.csv");


