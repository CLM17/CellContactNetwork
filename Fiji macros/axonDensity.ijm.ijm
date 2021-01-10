// User input
#@ File (label="Root", style="directory") root
#@ String (label="Well names (separated by commas)") wellListString 
#@ String (label="Names of the thresholded channels") chNamesString

close("*");

chNames = split(chNamesString, ",");
wellList = split(wellListString, ",");

// Check if thresholded images are there
for (m = 0; m < wellList.length; m++) {
	
	well = wellList[m];
	wellFolder = root + "\\" + well; 
	
	for (c = 0; c < chNames.length; c++) {
		thFile = wellFolder + "\\" + well+"_th_"+ chNames[c] +".tif";
		if(!(File.exists(thFile))){
			exit("No image called "+well+"_th_"+ chNames[c] +".tif was found.\nCheck whether all thresholded files exist and whether you entered the names of the thresholded images correctly.");
		}
	}
}

for (m = 0; m < wellList.length; m++) {

	well = wellList[m];
	wellFolder = root + "\\" + well; 
	
	// Open thresholded images
	channels = chNames.length;
	for (c = 0; c < channels; c++) {
		open(wellFolder + "/" + well+"_th_"+ chNames[c] +".tif");
		run("Select None");
	}
	
	// Empty ROI manager
	if (roiManager("count") > 0){
		roiManager("Deselect");
		roiManager("Delete");
	}
	selectWindow(well + "_th_dapi.tif");
	run("Select None");
	run("Analyze Particles...", "size=100-Infinity add");
	numNuclei = roiManager("count");
	roiManager("Deselect");
	roiManager("Delete");
	
	imageCalculator("AND create", well+"_th_dapi.tif", well+"_th_map2.tif");
	selectWindow("Result of "+well+"_th_dapi.tif");
	run("Marker-controlled Watershed", "input="+well+"_th_map2.tif marker=[Result of "+well+"_th_dapi.tif] mask="+well+"_th_dapi.tif binary use");
	close("Result of "+well+"_th_dapi.tif");
	
	selectWindow(well+"_th_map2-watershed.tif");
	fName = getTitle();
	run("8-bit");
	setThreshold(1, 255);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	
	run("Make Markers","outputdir="+wellFolder+" currentData="+fName);
	
	print("Counting the number of neurons...");
	run("Clear Results");
	run("Set Measurements...", "min redirect=None decimal=3");
	run("Measure");
	numNeurons = getResult("Max", nResults-1);
	close();
	
	print("Measuring area of MAP2...");
	run("Clear Results");
	run("Set Measurements...", "area area_fraction redirect=None decimal=3");
	selectWindow(well + "_th_map2.tif");
	run("Select None");
	run("Measure");
	area = getResult("Area", nResults-1);
	areaMAP2 = getResult("%Area", nResults-1) * area / 100;
	close(well + "_th_map2.tif");
	
	print("Measuring area of aTUB...");
	run("Clear Results");
	selectWindow(well + "_th_atubulin.tif");
	run("Select None");
	run("Measure");
	areaTubulin = getResult("%Area", nResults-1) * area / 100;
	close(well + "_th_atubulin.tif");
	
	results = "Num nuclei = " + d2s(numNuclei,0) + "\n" + "Num neurons = " + d2s(numNeurons,0) + "\n" + "Area MAP2 = " + d2s(areaMAP2,0) + "\n" + "Area aTUB = " + d2s(areaTubulin,0);
	print(results);
	
	File.saveString(results, wellFolder + "\\" + well + "_neuronMeasurements.txt"); 
	close("*");
}

print("All wells are processed");
