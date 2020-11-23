getDimensions(width, height, channels, slices, frames);
setThreshold(1, 1E30);
run("Convert to Mask");
run("Set Measurements...", "center redirect=None decimal=3");
run("Select None");
run("Analyze Particles...", "display clear add");
setBatchMode(true);
newImage("markers", "8-bit black", width, height, 1);

for (i = 0; i < nResults; i++) {
	x = Math.round( getResult("XM", i) );
	y =  Math.round( getResult("YM", i) );
	
	makeRectangle(x, y, 1, 1);
	run("Fill", "slice");
}
setBatchMode(false);