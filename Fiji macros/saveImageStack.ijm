// Get specifications 
#@ File (label="Input MFGTMP folder", style="directory") input_path
#@ File (label="Output folder", style="directory") output_path
#@ String (label="Well name") well
#@ int (label="Channel of interest") ch
#@ String (label="ATP concentration in mM") atp_str
#@ String (label="Start time of measurement (in seconds after ATP addition)") startT_str

close("*");
q = File.separator;
setBatchMode(true);

// Convert strings to decimal numbers
atp = parseFloat(atp_str);
startT = parseFloat(startT_str);

// Get MFGTMP name (name of input directory)
split_path = split(input_path,q);
mfgtmp_name = split_path[split_path.length-1];  // File name is last entry of path

// Init paths to files
tif_file_path = input_path + q + well;
kinetics_path = input_path + q + mfgtmp_name + "_kineticprotocol.log";
metadata_path = output_path + q + well + "_metadata.txt";
output_file = output_path + q + well + "_stack.tif";

// check inputs
if (!File.isDirectory(tif_file_path)) {
	exit("Sorry, there is no folder with TIFF images in the input MFGTMP folder.\nCheck whether you entered the correct folder.")
}
if (!File.exists(kinetics_path)) {
	exit("Sorry, the metadata file " + mfgtmp_name + "_kineticprotocol.log is missing in " + input_path + ".\nCheck the input MFGTMP folder.")
}
if (File.exists(output_file)) {
	showMessageWithCancel("Stack already created!","WARNING: An image stack with the name "+ well + "_stack.tif already exists in the output folder.\nDo you want to continue and overwrite the existing stack and metadata file?");
}
showMessageWithCancel("Check inputs","You entered the following parameters:\n\n[ATP] = "+d2s(atp,3)+" mM\nStart time = "+d2s(startT,0)+" s.\n\nPress OK to continue.");

// Get first input TIF file
file_list = getFileList(tif_file_path);
for (i = 0; i < file_list.length; i++) {
	file_name = file_list[i];
	if (indexOf(file_name, "d"+d2s(ch,0)) > 0){
		first_file_path = tif_file_path + q + file_name;
		file_name_split1 = split(file_name, "f");
		file_name_split2 = split(file_name_split1[1], "d");
		field_nr_str = file_name_split2[0];
		break;
	}
}

// Read file with kinetics protocol
kinetics_file = File.openAsString(kinetics_path);
kinetics_lines = split(kinetics_file, "\n");

// Get metadata from kinetics protocol
metadata = "RawDataFolder = " + input_path + "\n";
metadata = metadata + "Field = " + field_nr_str + "\n";
for (i = 0; i < kinetics_lines.length; i++) {
	split_line = split(kinetics_lines[i], "=");
	if (split_line.length > 0) {
		if (split_line[0] == "Name"){
			metadata = metadata + "ProtocolName = " + split_line[1] + "\n";
		}
		else if (split_line[0] == "TotalKineticsExecutionTime") {
			metadata = metadata + "DurationOfMeasurement = " + split_line[1] + " s\n";
		}
		else if (split_line[0] == "MinimumScanInterval") {
			metadata = metadata + "MinimumScanInterval = " + split_line[1] + " s\n";
		}
	}
}

// Add ATP concentration and start time of measurement to metadata file
metadata = metadata + "[ATP] = " + d2s(atp,3) + " mM\n";
metadata = metadata + "StartTime = " + d2s(startT,0) + " s\n";

// Save metadata file
File.saveString(metadata, metadata_path)
print("\n>>>> Saved metadata file: " + metadata_path);


// Load and save image sequence
print(">>>> Loading image sequence ...");
run("Image Sequence...", "open="+first_file_path+" file=d"+d2s(ch,0)+" sort");
run("8-bit");
run("Green");
run("Enhance Contrast", "saturated=0.35");
run("Apply LUT", "stack");
print(">>>> Saving image sequence ...");
saveAs("Tiff", output_file);
print(">>>> Saved image stack in " + output_file);
setBatchMode("exit and display");
