/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package com.mycompany.imagej;

import java.util.ArrayList;
import java.util.Iterator;
import ij.process.ImageProcessor;
import ij.ImagePlus;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.ParticleAnalyzer;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.TreeSet;
import ij.IJ;

import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import org.scijava.ui.UIService;
import net.imagej.Dataset;
import org.scijava.plugin.Parameter;
import java.io.File;
import org.scijava.plugin.Plugin;
import org.scijava.command.Command;
import net.imglib2.type.numeric.RealType;


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
@Plugin(type = Command.class, menuPath = "Plugins>Measure Cells")
public class MeasureCells<T extends RealType<T>> implements Command
{
	@Parameter(label = "Please select directory for output file", style = "directory")
    private File outputdir;
    @Parameter
    private Dataset currentData;
    @Parameter
    private UIService uiService;
    @Parameter
    private OpService opService;
    
    public void run() {
        
    	//final ImagePlus imp = IJ.getImage();
    	final ImagePlus imp = IJ.openImage("/Users/lukasvdh/Desktop/B02_wts.tif");
    	ImageProcessor ip = imp.getProcessor();
    	double max = ip.getMax();

    	// Initialize measurement arrays
    	ArrayList<Float> area = new ArrayList<Float>();
    	ArrayList<Float> major = new ArrayList<Float>();
    	ArrayList<Float> minor = new ArrayList<Float>();
    	ArrayList<Float> circularity = new ArrayList<Float>();
    	
    	ResultsTable rt = new ResultsTable();
    	double minSize = 0;
		double maxSize = Double.POSITIVE_INFINITY;
		
		//int options = ParticleAnalyzer.SHOW_RESULTS;
		int options = 0;
		
		//int measurements = Measurements.AREA+Measurements.CENTROID+Measurements.CIRCULARITY+Measurements.PERIMETER;
		int measurements = Measurements.AREA+Measurements.ELLIPSE+Measurements.CIRCULARITY+Measurements.CENTER_OF_MASS;
    	
		System.out.println("Entering for-loop");
		
		for(int i=1; i<max+1; i++) {
			
			ImageProcessor th = imp.getProcessor();
			th.setThreshold(i, i, ImageProcessor.RED_LUT);
			
			ImagePlus output = new ImagePlus("Output", th);
			//IJ.run(output,"Convert to Mask","");
			
			rt.reset();
			ParticleAnalyzer PA = new ParticleAnalyzer(options, measurements, rt, minSize, maxSize);
			PA.analyze(output);
			
			float[] a = rt.getColumn(ResultsTable.AREA);
			float[] mi = rt.getColumn(ResultsTable.MINOR);
			float[] ma = rt.getColumn(ResultsTable.MAJOR);
			float[] c = rt.getColumn(ResultsTable.CIRCULARITY);
			
			area.add(a[0]);
			minor.add(mi[0]);
			major.add(ma[0]);
			circularity.add(c[0]);
			
			if (i % 100 == 0) {
				System.out.println(String.valueOf(i+1));
			}
		}
		
		
		// Write measurement results to csv
    	PrintWriter pw = null;
        try {
            pw = new PrintWriter(new File(this.outputdir + "/cell_measurements.csv"));
        }
        
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        
        pw.print("Nr");
		pw.print(",");
		pw.print("Area");
		pw.print(",");
		pw.print("Minor");
		pw.print(",");
		pw.print("Major");
		pw.print(",");
		pw.print("Circularity");
		pw.print("\n");
        
        for (int i=0; i < area.size(); i++) {
			String a = String.format("%s",area.get(i));
			String mi = String.format("%s",minor.get(i));
			String ma = String.format("%s",major.get(i));
			String c = String.format("%s",circularity.get(i));

			pw.print(String.valueOf(i+1));
			pw.print(",");
			pw.print(a);
			pw.print(",");
			pw.print(mi);
			pw.print(",");
			pw.print(ma);
			pw.print(",");
			pw.print(c);
			pw.print("\n");
    	}
        pw.close();
        System.out.println("Done");
    }
    
    
    /**
     * This main function serves for development purposes.
     * It allows you to run the plugin immediately out of
     * your integrated development environment (IDE).
     *
     * @param args whatever, it's ignored
     * @throws Exception
     */
    public static void main(final String... args) throws Exception {
        // create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();

        // ask the user for a file to open
        final File file = ij.ui().chooseFile(null, "open");

        if (file != null) {
            // load the dataset
            final Dataset dataset = ij.scifio().datasetIO().open(file.getPath());
            
            // show the image
            ij.ui().show(dataset);

            // invoke the plugin
            ij.command().run(MeasureCells.class, true);
        }
    }
}
