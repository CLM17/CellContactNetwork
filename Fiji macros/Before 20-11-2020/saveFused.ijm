root = "M:/tnw/bn/dm/Shared/Lukas/BEP/Experiments/WKS024";
well = "D04";
w = "8";

close("*")
run("Grid/Collection stitching", "type=[Grid: column-by-column] order=[Down & Right                ] grid_size_x="+w+" grid_size_y="+w+" tile_overlap=1 first_file_index_i=0 directory="+root+"/"+well+"/tiles file_names=tile_{ii}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap display_fusion computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");
run("16-bit");

fname = root+"/"+ well+"/"+well+"_fused.tif";
saveAs("Tiff", fname);
