/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package com.mycompany.imagej;

import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.ParticleAnalyzer;
import ij.process.ImageProcessor;
import ij.measure.ResultsTable;
import ij.measure.Measurements;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


/**
 * This example illustrates how to create an ImageJ {@link Command} plugin.
 * <p>
 * The code here is a simple Gaussian blur using ImageJ Ops.
 * </p>
 * <p>
 * You should replace the parameter fields with your own inputs and outputs,
 * and replace the {@link run} method implementation with your own logic.
 * </p>
 */
@Plugin(type = Command.class, menuPath = "Plugins>Make Markers")
public class MakeMarkers<T extends RealType<T>> implements Command 
{
    //
    // Feel free to add more parameters here...
    //
	
	@Parameter(label = "Please select directory for output file", style = "directory")
    private File outputdir;

    @Parameter
    private Dataset currentData;
	
    @Parameter
    private UIService uiService;

    @Parameter
    private OpService opService;

    public void run() {
    	
    	// Get image and make sure it's 8 bit
    	final ImagePlus imp = IJ.getImage();  
    	//final ImagePlus imp = IJ.openImage("/Users/lukasvdh/Desktop/B02_th_dapi_cropped.tif");
    	
    	IJ.run(imp, "8-bit", "");
    	
    	// Initialize particle analysis
    	ResultsTable rt = new ResultsTable();
    	rt.reset();
    	
		double minSize = 0;
		double maxSize = Double.POSITIVE_INFINITY;
		
		//int options = ParticleAnalyzer.SHOW_RESULTS;
		int options = 0;
		int measurements = Measurements.CENTER_OF_MASS;
		
		// Do particle analysis
		ParticleAnalyzer PA = new ParticleAnalyzer(options, measurements, rt, minSize, maxSize);
    	PA.analyze(imp);
    	    	
    	// Get centers of mass
    	float[] xCM = rt.getColumn(ResultsTable.X_CENTER_OF_MASS);				
    	float[] yCM = rt.getColumn(ResultsTable.Y_CENTER_OF_MASS);
    	
    	//Create empty image
    	int width = imp.getWidth();
    	int height = imp.getHeight();
    	
    	ImagePlus out = IJ.createImage("Markers", "16-bit black", width, height, 1);
    	ImageProcessor ip_out = out.getProcessor();
    	
    	// Loop through centers of mass and put a pixel on all COMs
    	for (int i=0; i < xCM.length; i++) {
    		
    		int x = (int) xCM[i];
    		int y = (int) yCM[i];
    		ip_out.putPixel(x, y, i+1);
    		
    	}
    	
    	ImagePlus output = new ImagePlus("Output", ip_out);
    	output.show();
    	
    	
    	// Write COM results to csv
    	PrintWriter pw = null;
        try {
            pw = new PrintWriter(new File(this.outputdir + "/nuclei_com.csv"));
        }
        
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        for (int i=0; i < xCM.length; i++) {
			String x = String.format("%s",xCM[i]);
			String y = String.format("%s",yCM[i]);

			pw.print(String.valueOf(i+1));
			pw.print(",");
			pw.print(x);
			pw.print(",");
			pw.print(y);
			pw.print("\n");
    	}
        pw.close();
    	
    }

    /**
     * This main function serves for development purposes.
     * It allows you to run the plugin immediately out of
     * your integrated development environment (IDE).
     *
     * @param args whatever, it's ignored
     * @throws Exception
     */
    //public static void main(final String... args) throws Exception {
    //    // create the ImageJ application context with all available services
    //    final ImageJ ij = new ImageJ();
    //    ij.ui().showUI();

        // ask the user for a file to open
     //   final File file = ij.ui().chooseFile(null, "open");

    //    if (file != null) {
            // load the dataset
    //        final Dataset dataset = ij.scifio().datasetIO().open(file.getPath());
            
            // show the image
    //        ij.ui().show(dataset);

            // invoke the plugin
    //        ij.command().run(MakeMarkers.class, true);
    //    }
    //}

}
