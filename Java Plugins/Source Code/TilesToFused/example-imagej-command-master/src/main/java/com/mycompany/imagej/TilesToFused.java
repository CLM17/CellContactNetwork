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
import ij.*;

import ij.plugin.HyperStackConverter;
import ij.plugin.HyperStackReducer;
import java.util.List;
import java.util.Vector;
import java.util.Collections;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import ij.ImagePlus;
import ij.process.ImageProcessor;


import java.io.File;
import java.io.IOException;



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
@Plugin(type = Command.class, menuPath = "Plugins>Tiles to Fused")
public class TilesToFused<T extends RealType<T>> implements Command 
{
    //
    // Feel free to add more parameters here...
    //
	
	@Parameter(label = "Please select tile folder", style = "directory", required = true)
    private File tileDirectory;
	
	@Parameter(label = "Please select well locations file", required = true)
    private File wellLocationsTextFile;

    @Parameter
    private int fusedwidth;
    
    @Parameter
    private int fusedheight;
	
    @Parameter
    private UIService uiService;

    @Parameter
    private OpService opService;

    public void run() {    	
    	
    	String tileFolder = tileDirectory.toString();
    	
    	// Get dimensions of fused image
    	//final ImagePlus imp = IJ.openImage("/Users/lukasvdh/Documents/BEP/Experiments/WKS024/10x/B02/B02_fused.tif");
    	int width = fusedwidth;
    	int height = fusedheight;
    	
    	// Get number of channel in the first tile: this will be the number of channels in the fused image
    	final ImagePlus tile0 = IJ.openImage(tileFolder+"/tile_00.tif");
    	int ch = tile0.getNChannels();
    	
    	// Create empty array for output image
    	ImagePlus [] out = new ImagePlus[ch];
    	
    	// Initiate 1 empty image for each channel and store in out
    	for (int c=0; c<ch; c++) {
    		out[c] = IJ.createImage(Integer.toString(c+1) , "8-bit black", width, height, 1);
    	}
    	
    	// Read .txt file with tile configuration
    	List<String> lines = Collections.emptyList();
    	try {
    		
    		lines = Files.readAllLines(wellLocationsTextFile.toPath(), StandardCharsets.UTF_8); 
    	}
			catch (IOException e) {
		  		e.printStackTrace(); 
		    } 
  
    	// Loop over tiles
    	for(int i=4; i<lines.size(); i++) {
      		
    		// Get tile name and coordinates from text file
    		String line = lines.get(i);
    		String [] parts = line.split(";");
    		
    		String tileName = parts[0];
    		String coordinates = parts[2];
    		String [] splitCoordinates = coordinates.split(",");
    		int x = (int) Float.parseFloat( splitCoordinates[0].substring(2) );										// Remove "("
    		int y = (int) Float.parseFloat( splitCoordinates[1].substring(0, splitCoordinates[1].length() - 1) );	// Remove ")"
    		
    		// Open and split tile
    		ImagePlus tile = IJ.openImage(tileFolder + "/" + tileName);
    		//IJ.run(tile, "8-bit", "");
    		ImagePlus [] tileChannels = split(tile);
    		
    		int ht = tile.getHeight();
    		int wt = tile.getWidth();

    		// Loop over channels
    		for (int c=0; c<ch; c++) {
    			
    			ImageProcessor tileIP = tileChannels[c].getProcessor();
    			ImageProcessor outIP = out[c].getProcessor();
    			
    			// Loop over pixels
    			for(int xi=0; xi<wt; xi++) {
    				
    				for(int yi=0; yi<ht; yi++) {
    					
    					float p = tileIP.getf(xi, yi);
    					outIP.putPixel(x + xi, y + yi, (int) p);
    					
    				}
    			}
    			
    			out[c] = new ImagePlus(Integer.toString(c+1), outIP);
    			
    		}
    	}
     	
    	// Create output image stack
    	ImageStack stack = new ImageStack(width,height);
    	
    	for(int c=0; c<ch; c++) {
    		stack.addSlice(out[c].getProcessor());
    	}

    	// Convert stack to channels
    	ImagePlus outputImage = new ImagePlus("fused", stack);
    	outputImage = HyperStackConverter.toHyperStack(outputImage, ch, 1, 1);
 
		outputImage.show();
		
		
    }
        
        
    /** Splits the specified image into separate channels. */
    public static ImagePlus[] split(ImagePlus imp) {
        
        int width = imp.getWidth();
        int height = imp.getHeight();
        int channels = imp.getNChannels();
        int slices = imp.getNSlices();
        int frames = imp.getNFrames();
        int bitDepth = imp.getBitDepth();
        int size = slices*frames;
        Vector images = new Vector();
        HyperStackReducer reducer = new HyperStackReducer(imp);
        for (int c=1; c<=channels; c++) {
            ImageStack stack2 = new ImageStack(width, height, size); // create empty stack
            stack2.setPixels(imp.getProcessor().getPixels(), 1); // can't create ImagePlus will null 1st image
            ImagePlus imp2 = new ImagePlus("C"+c+"-"+imp.getTitle(), stack2);
            stack2.setPixels(null, 1);
            imp.setPosition(c, 1, 1);
            imp2.setDimensions(1, slices, frames);
            imp2.setCalibration(imp.getCalibration());
            reducer.reduce(imp2);
            if (imp.isComposite() && ((CompositeImage)imp).getMode()==IJ.GRAYSCALE)
                IJ.run(imp2, "Grays", "");
            if (imp2.getNDimensions()>3)
                imp2.setOpenAsHyperStack(true);
            images.add(imp2);
        }
        ImagePlus[] array = new ImagePlus[images.size()];
        return (ImagePlus[])images.toArray(array);
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
        // create the ImageJ application context with all available services
    //    final ImageJ ij = new ImageJ();
    //    ij.ui().showUI();

        // ask the user for a file to open
    //    final File file = ij.ui().chooseFile(null, "open");

    //    if (file != null) {
            // load the dataset
    //        final Dataset dataset = ij.scifio().datasetIO().open(file.getPath());
            
            // show the image
    //        ij.ui().show(dataset);

            // invoke the plugin
    //        ij.command().run(TilesToFused.class, true);
    //    }
    //}

}
