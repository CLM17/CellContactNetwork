
run("Set Measurements...", "min redirect=None decimal=3");
run("Measure");
max = getResult("Max", nResults-1);

setBatchMode(true);
run("Clear Results");

run("Set Measurements...", "area centroid center perimeter fit feret's area_fraction redirect=None decimal=3");
for (i = 1; i < max+1; i++) {
	selectWindow("B02_wts.tif");
	run("Duplicate...", " ");
	setThreshold(i, i);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Analyze Particles...", "add");
	close();
}
roiManager("Show All");
roiManager("Measure");
saveAs("Results", output + "/" + "cell_properties.csv");