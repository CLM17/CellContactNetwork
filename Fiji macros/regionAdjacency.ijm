close("*");
open("M:/tnw/bn/dm/Shared/Lukas/BEP/Experiments/WKS024/B02/B02_th_dapi.tif");
open("M:/tnw/bn/dm/Shared/Lukas/BEP/Experiments/WKS024/B02/B02_th_marker.tif");

selectWindow("B02_th_dapi.tif");

run("Ultimate Points");
run("Marker-controlled Watershed", "input=B02_th_marker.tif marker=B02_th_dapi.tif mask=B02_th_marker.tif binary");
run("Duplicate");
run("3-3-2 RGB");

run("Region Adjacency Graph", "show image=B02_th_marker-watershed-1.tif");

