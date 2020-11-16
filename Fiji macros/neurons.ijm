close("*");
folder = "M:/tnw/bn/dm/Shared/Lukas/BEP/Experiments/WKS026/test/D03/";
open("M:/tnw/bn/dm/Shared/Lukas/BEP/Experiments/WKS026/test/D03/D03_fused.tif");
run("Split Channels");

//setBatchMode("hide");

selectWindow("C3-D03_fused.tif");
rename("dapi");
setAutoThreshold("Li dark");
setOption("BlackBackground", true);
run("Convert to Mask");

selectWindow("C2-D03_fused.tif");
rename("mapk2");
setAutoThreshold("Li dark");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Morphological Filters", "operation=Opening element=Disk radius=3");
run("Morphological Filters", "operation=Closing element=Disk radius=3");

 
close("mapk2");
close("mapk2-Opening");

selectWindow("C1-D03_fused.tif");
rename("atub")
setAutoThreshold("Li dark");
setOption("BlackBackground", true);
run("Convert to Mask");

imageCalculator("Multiply", "dapi","mapk2-Opening-Closing");
run("Morphological Filters", "operation=Opening element=Disk radius=5");
run("Morphological Filters", "operation=Closing element=Disk radius=5");

close("dapi");
close("dapi-Opening");

fName = "dapi-Opening-Closing";
selectWindow(fName);
run("Make Markers","outputdir="+folder+" currentData="+fName);

run("Marker-controlled Watershed", "input=mapk2-Opening-Closing marker=Output mask=mapk2-Opening-Closing calculate use");

selectWindow("Output");
setThreshold(1, 1E30);
setOption("BlackBackground", true);
run("Convert to Mask");

run("Set Measurements...", "center redirect=None decimal=3");
run("Analyze Particles...", "add");
run("Clear Results");
roiManager("Measure");

i = 1;

// Make MAKK2 image with curent MAPK2 subtracted
selectWindow("mapk2-Opening-Closing-watershed");
run("Duplicate...", " ");
setThreshold(i, i);
setOption("BlackBackground", true);
run("Convert to Mask");
run("Divide...", "value=255");
run("Multiply...", "value="+d2s(i,0));

imageCalculator("Subtract create", "mapk2-Opening-Closing-watershed","mapk2-Opening-Closing-watershed-1");
close("mapk2-Opening-Closing-watershed-1");
rename("mapk2_without_i");

xc = round(getResult("XM", i-1));
yc = round(getResult("YM", i-1));

w = getWidth();
h = getHeight();
newImage("growing_seed", "8-bit black", w, h, 1);
setPixel(xc, yc, 1);

for (t = 0; t < 500; t++) {
	selectWindow("growing_seed");
	run("Morphological Filters", "operation=Dilation element=Disk radius=20");
	close("growing_seed");
	rename("growing_seed");
	imageCalculator("AND", "growing_seed","atub");

	imageCalculator("Multiply create 32-bit", "growing_seed","mapk2_without_i");
	getStatistics(area, mean, min, max, std, histogram);
	close("Result of growing_seed");
	print(t);
	
	if (max > 0){

		print(max);
		selectWindow("mapk2_without_i");
		run("Duplicate...", " ");
		setThreshold(max, max);
		setOption("BlackBackground", true);
		run("Convert to Mask");

		imageCalculator("Subtract", "atub", "mapk2_without_i-1");

		selectWindow("mapk2_without_i-1");
		run("Divide...", "value=255");
		run("Multiply...", "value="+d2s(max,0));
		
		imageCalculator("Subtract", "mapk2_without_i", "mapk2_without_i-1");
		close("mapk2_without_i-1");

		print(d2s(i,0)+" connects to "+d2s(max,0));
	}
}