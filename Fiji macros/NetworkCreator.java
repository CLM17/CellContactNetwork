package com.tudelft.dm;

import java.io.File;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.TreeSet;

import net.imagej.ImageJ;
import ij.IJ;
import ij.process.ImageProcessor;
import ij.ImagePlus;
import ij.WindowManager;
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;

public class Network_Creator implements PlugIn{
	
	public void main(String[] args) {
		
		//Path currentRelativePath = Paths.get("");
		//String s = currentRelativePath.toAbsolutePath().toString();
		//System.out.println(s);
		
		//String filenameWts = "B02_wts.tif";
		//String filenameWts = "B03_th_watershed.tif";
		
		ImagePlus imp = IJ.getImage();
		
		//ArrayList list = new ArrayList(); 
		TreeSet<LabelPair> list = new TreeSet<LabelPair>();
		
		ImageProcessor ip = imp.getChannelProcessor();
		int count = 0;

		System.out.println("Entering for-loop");
		// loop over image pixels, but exclude the edge
		for (int x=1; x<ip.getWidth()-1; x++) {
			for (int y=1; y<ip.getHeight()-1; y++) {
				
				int center = (int) ip.getf(x,y);
				
				if (center > 0) {
					int north = (int) ip.getf(x,y-1);
					int west = (int) ip.getf(x-1,y);
					int northwest = (int) ip.getf(x-1, y-1);
					count = count + 1;
					
					if (north > 0 && center != north) {
						//String toNorth = Integer.toString(center) + "->" + Integer.toString(north);
						//String fromNorth = Integer.toString(north) + "->" + Integer.toString(center);
						
						LabelPair pairN = new LabelPair(center, north);
						
						if (!list.contains(pairN)){
							list.add(pairN);
						}
						//if (!list.contains(fromNorth)){
						//	list.add(fromNorth);
						//}
					}
					
					if (west > 0 && center != west) {
						//String toWest = Integer.toString(center) + "->" + Integer.toString(west);
						//String fromWest = Integer.toString(west) + "->" + Integer.toString(center);
						
						LabelPair pairW = new LabelPair(center, west);
						
						if (!list.contains(pairW)){
							list.add(pairW);
						}
						//if (!list.contains(fromWest)){
						//	list.add(fromWest);
						//}
					}
					
					if (northwest > 0 && center != northwest) {
						//String toNorthwest = Integer.toString(center) + "->" + Integer.toString(northwest);
						//String fromNorthwest = Integer.toString(northwest) + "->" + Integer.toString(center);
						
						LabelPair pairNW = new LabelPair(center, northwest);
						
						if (!list.contains(pairNW)){
							list.add(pairNW);
						}
						//if (!list.contains(fromNorthwest)){
						//	list.add(fromNorthwest);
						//}
					}
				}
										
			}
		}
		
		System.out.println("For-loop is done");
		imp.show();
		
		 for (LabelPair pair : list) {
			 System.out.println("(" + pair.label1 + "," + pair.label2 + ")");
		 }
		
		// Write results to csv
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(new File("output.csv"));
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		for (LabelPair pair : list) {
			pw.print(pair.label1 + "->" + pair.label2); 
			pw.print("\n");
		}
		
		pw.close();
		System.out.println("Done");
	}
	
	/**

	 * Used to stores the adjacency information between two regions. In order to
	 * ensure symmetry of the relation, the value of label1 field always
	 * contains the lower label, while the value of label2 always contains the
	 * highest label.
	 */

	public static final class LabelPair implements Comparable <LabelPair>

	{
		public int label1;
		public int label2;

		public LabelPair(int label1, int label2)

		{

			if (label1 < label2) 

			{
				this.label1 = label1;
				this.label2 = label2;
			}

			else

			{
				this.label1 = label2;
				this.label2 = label1;
			}

		}
		
		@Override
		public int compareTo(LabelPair pair) {

			if (this.label1 < pair.label1)
				return -1;

			if (this.label1 > pair.label1)
				return +1;

			if (this.label2 < pair.label2)
				return -1;

			if (this.label2 > pair.label2)
				return +1;

			return 0;

		}
	}

}
