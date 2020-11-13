//open("M:/tnw/bn/dm/Shared/Lukas/BEP/Experiments/Test/B02/B02_wts.tif");

selectWindow("B02_wts.tif");
run("Duplicate...", " ");
setThreshold(1, 1);
setOption("BlackBackground", true);
run("Convert to Mask");

run("Clear Results");
run("Analyze Particles...", "add");