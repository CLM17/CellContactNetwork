#@ File (label="root", style="directory") root

wells = newArray(12);
wells[0] = "C02";
wells[1] = "C03";
wells[2] = "C04";
wells[3] = "C05";
wells[4] = "C06";
wells[5] = "C07";
wells[6] = "D02";
wells[7] = "D03";
wells[8] = "D04";
wells[9] = "D05";
wells[10] = "D06";
wells[11] = "D07";

for (i = 0; i < 12; i++) {

	well = wells[i];
	open(root+"/"+well+"/"+well+"_fused.tif");

	title1 = "Adjust B/C.";
	message1 = "Adjust B/C.\nPress OK when you are done.";
	waitForUser(title1, message1);

	run("RGB Color");
	setTool("oval");
	
	title2 = "Make oval selection.";
	message2 = "Fit the circular selection to the well borders while holding ctrl+shift.\nPress OK when you are done.";
	waitForUser(title2, message2);
	
	setBackgroundColor(0, 0, 0);
	run("Clear Outside");
	saveAs(root+"/"+well+"/"+well+"_fused_RGB.tif");
	close("*");
}